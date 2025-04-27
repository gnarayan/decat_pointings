# Stores all airmass calculations
import os
from pathlib import Path

import astropy.units as u
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from astropy.coordinates import SkyCoord
from astropy.time import Time
import sys

#from .helpers import ctio_location
from astroplan import Observer


def ctio_location(
    temperature=5.0 * u.deg_C,
    pressure=780.0 * u.mbar,
    relative_humidity=0.5,
):
    """Generate Observer (which inherits from EarthLocation)
    for CTIO.
    """
    return Observer.at_site(
        "CTIO",
        timezone="Etc/GMT+3",
        pressure=pressure,
        temperature=temperature,
        relative_humidity=relative_humidity,
    )


class AirmassCalculator:
    """Class that performs airmass calculations and caches
    for later access."""

    def __init__(self, cache_dir, location=None):
        """Initializes location airmass is calculated from."""
        if location is None:
            self._loc = ctio_location()
        else:
            self._loc = location

        pref_a = self._loc.pressure.value
        pref_b = self._loc.temperature.value
        pref_c = self._loc.relative_humidity
        file_prefix = f"{pref_a}_{pref_b}_{pref_c}"
        
        print(f'airmass grid cache dir: {cache_dir}')
        
        #cache_dir = os.path.join(Path(__file__).parent.parent.parent.absolute(), "data/airmass_grids")
        os.makedirs(cache_dir, exist_ok=True)
        self._cache_fn = os.path.join(cache_dir, f"{file_prefix}.feather")

        # hour angle / declination grid
        # here we only go to HA = 5:15:00 because that's where DECam cuts off
        ha_deg_range = np.linspace(0.0, 78.75, num=1_000)  # symmetric so only positives
        dec_range = np.linspace(-90.0, 40.0, num=1_000)  # above 40deg not visible
        _ha_grid, _dec_grid = np.meshgrid(ha_deg_range, dec_range)
        self._ha_grid = _ha_grid.ravel()
        self._dec_grid = _dec_grid.ravel()

    def _get_hour_angles(self, ra, times):
        """Convert RA value to hour angles."""
        try:
            times_valid = [isinstance(t, Time) for t in times]
            has = np.nan * np.ones(len(times))
            if np.sum(times_valid) == 0:
                return has
            has[times_valid] = (
                self._loc.local_sidereal_time(times[times_valid]).to(u.deg).value - ra.to(u.deg).value
            )
            return has
        except:
            if not isinstance(times, Time):
                return np.nan
            return self._loc.local_sidereal_time(times).to(u.deg).value - ra.to(u.deg).value

    def generate_airmass_grid(self):
        """From grid of hour angles and declinations, generate"""
        dummy_time = Time.now()  # just to convert from HA to RA and back

        ra_values = (self._loc.local_sidereal_time(dummy_time).degree - self._ha_grid) % 360
        sky_coords = SkyCoord(ra=ra_values * u.deg, dec=self._dec_grid * u.deg, frame="icrs")

        # Convert to AltAz frame (to then retrieve airmass)
        airmasses = self._loc.altaz(time=dummy_time, target=sky_coords).secz
        airmasses[airmasses < 0.0] = np.nan  # remove invalid airmasses

        self._cache = pd.DataFrame(
            np.array([self._ha_grid, self._dec_grid, airmasses]).T, columns=["ha", "dec", "airmass"]
        )
        self._ha_range = self._cache.ha.unique()[np.newaxis, :]
        self._cache.to_feather(self._cache_fn)

    def load_airmass_grid(self):
        """Load and return airmass grid."""
        if not os.path.exists(self._cache_fn):
            print('WARNING!!! the airmass grid is not in the cache dir, creating it!!!')
            self.generate_airmass_grid() 
            #raise FileNotFoundError("Cache not yet created! Please run self.generate_airmass_grid() first.")
        self._cache = pd.read_feather(self._cache_fn)
        self._ha_range = self._cache.ha.unique()[np.newaxis, :]
        return self._cache.copy()

    def _query_from_ha_dec(self, ha, dec):
        """Query cache from hour angle(s) and SINGLE declination."""
        if not np.isscalar(dec):
            raise TypeError("dec must be a scalar")

        if dec >= 40.0:  # out of bounds
            return np.nan * np.ones(len(ha))

        # take absolute values, airmasses symmetric around HA
        abs_ha = np.atleast_1d(np.abs(ha))
        oob_mask = np.isnan(abs_ha) | (abs_ha >= 78.75)  # out of bounds for DECam

        # Find closest hour angles
        min_dec_diff = (self._cache.dec - dec).abs().min()
        dec_df = self._cache.loc[(self._cache.dec - dec).abs() == min_dec_diff, :].copy()
        differences = np.abs(abs_ha[~oob_mask, np.newaxis] - self._ha_range)
        indices = np.argmin(differences, axis=1)
        nearest_hour_angles = self._ha_range[0, indices]

        ha_df = pd.DataFrame({"ha": nearest_hour_angles, "original_index": range(len(nearest_hour_angles))})
        dec_df = dec_df.sort_values("ha").reset_index(drop=True)
        merged_df = pd.merge(ha_df, dec_df, on="ha", how="left")
        merged_df = merged_df.sort_values("original_index").reset_index(drop=True)

        # Extract airmasses
        airmasses = merged_df["airmass"].to_numpy()
        airmasses_full = np.nan * np.ones(len(abs_ha))
        airmasses_full[~oob_mask] = airmasses

        return airmasses_full

    def _separation_squared(self, ra_arr, dec_arr, ra2, dec2):
        """Helper function to calculate separation
        directly from floats in degrees. Assumes
        small angle approximation.
        """
        return ((ra_arr - ra2) * np.cos((dec_arr + dec2) * np.pi / 360.0)) ** 2 + (dec_arr - dec2) ** 2

    def _query_direct(self, times, ra, dec):
        """Directly query airmass from RA/dec.
        For comparison to caching runtime.
        """
        coord = SkyCoord(ra=ra, dec=dec)
        airmasses = self._loc.altaz(time=times, target=coord).secz
        airmasses[(airmasses < 0.0) | (airmasses > 5.0)] = np.nan
        return airmasses

    def query(self, times, ra, dec):
        """Query from cache."""
        hour_angles = self._get_hour_angles(ra, times)
        airmasses = self._query_from_ha_dec(hour_angles, dec.to(u.deg).value)
        return airmasses

    def observable_range(self, ra, dec, start_time, end_time, max_airmass=1.8):
        """Find the min and max observable time within a certain
        time range such that the object is above a given airmass.
        """
        time_linspace = (
            start_time + np.linspace(0, (end_time - start_time).to(u.hour).value, num=500) * u.hour
        )
        hour_angles = self._get_hour_angles(ra, time_linspace)
        airmasses = self._query_from_ha_dec(hour_angles, dec.to(u.deg).value)
        visible_times = time_linspace[airmasses < max_airmass]
        if len(visible_times) == 0:
            return start_time, start_time

        return visible_times[0], visible_times[-1]

    def _plot_single_airmass(self, ax, ra, dec, times, **plot_kwargs):
        """Add airmass plot to ax, assuming
        target is at RA/Dec. If times is None, plot
        across entire 24 hours from now.
        """
        time_arr = np.atleast_1d(times)
        hour_angles = self._get_hour_angles(ra, time_arr)
        airmasses = self._query_from_ha_dec(hour_angles, dec.to(u.deg).value)

        ax.plot(time_arr.to_datetime(), airmasses, **plot_kwargs)

        return ax

    def plot_airmass(self, df, start_time, end_time, tz_shifts):
        """Plot airmass for entire dataframe of targets, highlighting
        current observing windows.
        """
        fig = plt.figure(figsize=(10, 4))
        ax1 = fig.add_subplot(111)
        ax2 = ax1.twiny()
        ax3 = ax1.twiny()
        ax1.invert_yaxis()
        ax1.set_ylim([2.0, 0.8])
        ax1.set_ylabel("Relative Air Mass")
        ax1.set_xlabel("Local Time")
        ax1.grid(True)
        ax1.set_axisbelow(True)

        length_of_night = 500

        tz_shift_sidereal, tz_shift_extra = tz_shifts

        full_night_linspace = (
            start_time
            + np.linspace(0, (end_time - start_time).to(u.hour).value, num=length_of_night) * u.hour
        )
        full_night_sidereal = full_night_linspace + tz_shift_sidereal
        full_night_extra = full_night_linspace + tz_shift_extra

        ax2.plot(full_night_linspace.to_datetime(), np.zeros(length_of_night))
        ax2.set_xlabel("UTC")
        ax2.get_xaxis().set_major_formatter(mdates.DateFormatter("%H:%M"))

        num_ticks = 7
        nn = round(length_of_night / num_ticks)
        ax3_ind = [i * nn for i in range(num_ticks)]

        ax3.plot(full_night_sidereal.to_datetime(), np.zeros(length_of_night))
        ax3_ind.remove(0)
        ax3.set_xlabel("LST")
        ax3.xaxis.set_ticks_position("bottom")
        ax3.xaxis.set_label_position("bottom")

        ax3.set_xticks(full_night_sidereal[ax3_ind].to_datetime())
        ax3.get_xaxis().set_major_formatter(mdates.DateFormatter("%H:%M"))

        # Offset the twin axis below the host
        ax3.spines["bottom"].set_position(("axes", -0.18))

        colormap = sns.color_palette("husl", len(df.program.unique()))
        program_colormap = {k: v for (k, v) in zip(df.program.unique(), colormap)}

        for k in program_colormap:
            #total_time = df.loc[df.program == k, "dur"].sum()
            label = f"{k}\nTotal time: X min"
            ax1.plot(
                full_night_extra.to_datetime(),
                np.zeros(length_of_night),
                color=program_colormap[k],
                linewidth=3,
                label=label,
            )

        ax1.get_xaxis().set_major_formatter(mdates.DateFormatter("%H:%M"))

        for row in df.itertuples():
            col = program_colormap[row.program]

            if not pd.isna(row.order):
                st = row.t_start_hidden
                et = row.t_end_hidden

                obs_linspace = st + np.linspace(0, (et - st).to(u.hour).value, num=500) * u.hour
                self._plot_single_airmass(ax2, row.ra*u.deg, row.dec*u.deg, times=obs_linspace, linewidth=3.0, color=col)

                self._plot_single_airmass(
                    ax2, row.ra*u.deg, row.dec*u.deg, times=full_night_linspace, linewidth=3.0, color=col, alpha=0.1
                )

            else:
                self._plot_single_airmass(
                    ax2,
                    row.ra*u.deg,
                    row.dec*u.deg,
                    times=full_night_linspace,
                    linewidth=3.0,
                    color=col,
                    linestyle="dotted",
                    alpha=0.3,
                )

        leg = ax1.legend(bbox_to_anchor=(1.01, 1.015), loc="upper left", ncol=1, prop={"size": 8})

        # set the linewidth of each legend object
        #for legobj in leg.legendHandles:
        for legobj in leg.legend_handles:
            legobj.set_linewidth(3.0)

        return fig

    def plot_airmassDELME(self, df, start_time, end_time, tz_shifts):
        """Plot airmass for entire dataframe of targets, highlighting
        current observing windows.
        """
        fig = plt.figure(figsize=(10, 4))
        ax1 = fig.add_subplot(111)
        ax2 = ax1.twiny()
        ax3 = ax1.twiny()
        ax1.invert_yaxis()
        ax1.set_ylim([2.0, 0.8])
        ax1.set_ylabel("Relative Air Mass")
        ax1.set_xlabel("Local Time")
        ax1.grid(True)
        ax1.set_axisbelow(True)

        length_of_night = 500

        tz_shift_sidereal, tz_shift_extra = tz_shifts

        full_night_linspace = (
            start_time
            + np.linspace(0, (end_time - start_time).to(u.hour).value, num=length_of_night) * u.hour
        )
        full_night_sidereal = full_night_linspace + tz_shift_sidereal
        full_night_extra = full_night_linspace + tz_shift_extra

        ax2.plot(full_night_linspace.to_datetime(), np.zeros(length_of_night))
        ax2.set_xlabel("UTC")
        ax2.get_xaxis().set_major_formatter(mdates.DateFormatter("%H:%M"))

        num_ticks = 7
        nn = round(length_of_night / num_ticks)
        ax3_ind = [i * nn for i in range(num_ticks)]

        ax3.plot(full_night_sidereal.to_datetime(), np.zeros(length_of_night))
        ax3_ind.remove(0)
        ax3.set_xlabel("LST")
        ax3.xaxis.set_ticks_position("bottom")
        ax3.xaxis.set_label_position("bottom")

        ax3.set_xticks(full_night_sidereal[ax3_ind].to_datetime())
        ax3.get_xaxis().set_major_formatter(mdates.DateFormatter("%H:%M"))

        # Offset the twin axis below the host
        ax3.spines["bottom"].set_position(("axes", -0.18))

        colormap = sns.color_palette("husl", len(df.program.unique()))
        program_colormap = {k: v for (k, v) in zip(df.program.unique(), colormap)}

        for k in program_colormap:
            #total_time = df.loc[df.program == k, "dur"].sum()
#            label = f"{k}\nTotal time: X min"
            label = f"{k}"
            ax1.plot(
                full_night_extra.to_datetime(),
                np.zeros(length_of_night),
                color=program_colormap[k],
                linewidth=3,
                label=label,
            )

        ax1.get_xaxis().set_major_formatter(mdates.DateFormatter("%H:%M"))

        for row in df.itertuples():
            col = program_colormap[row.program]

            if not pd.isna(row.order):
                st = row.t_start_hidden
                et = row.t_end_hidden

                obs_linspace = st + np.linspace(0, (et - st).to(u.hour).value, num=500) * u.hour
                self._plot_single_airmass(ax2, row.raX, row.decX, times=obs_linspace, linewidth=3.0, color=col)

                self._plot_single_airmass(
                    ax2, row.raX, row.decX, times=full_night_linspace, linewidth=3.0, color=col, alpha=0.1
                )

            else:
                self._plot_single_airmass(
                    ax2,
                    row.raX,
                    row.decX,
                    times=full_night_linspace,
                    linewidth=3.0,
                    color=col,
                    linestyle="dotted",
                    alpha=0.3,
                )

        leg = ax1.legend(bbox_to_anchor=(1.01, 1.015), loc="upper left", ncol=1, prop={"size": 8})

        # set the linewidth of each legend object
        #for legobj in leg.legendHandles:
        for legobj in leg.legend_handles:
            legobj.set_linewidth(3.0)

        return fig
