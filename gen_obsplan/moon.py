#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Kaylee Soto
"""
import astropy.units as u
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from astropy.coordinates import SkyCoord, get_body


class MoonCalculator:
    """Class that performs airmass calculations and caches
    for later access."""

    def __init__(self):
        """Initializes location airmass is calculated from."""
        pass

    def get_moon_distance(self, ra, dec, times):
        """Return distance to certain coord."""
        moon_coords = get_body("moon", times)
        coord = SkyCoord(ra=ra, dec=dec)
        return moon_coords.separation(coord).to(u.deg).value

    def _plot_single_moon_dist(self, ax, ra, dec, times, **plot_kwargs):
        """Plot curve for single target versus moon distance."""
        time_arr = np.atleast_1d(times)
        distances = self.get_moon_distance(ra, dec, time_arr)

        ax.plot(time_arr.to_datetime(), distances, **plot_kwargs)

        return ax

    def plot_moon_distances(self, df, start_time, end_time, tz_shifts):
        """Plot airmass for entire dataframe of targets, highlighting
        current observing windows.
        """
        fig = plt.figure(figsize=(10, 4))
        ax1 = fig.add_subplot(111)
        ax2 = ax1.twiny()
        ax3 = ax1.twiny()
        ax1.set_ylim([0.0, 180.0])
        ax1.set_ylabel("Moon Distance (deg)")
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

        ax2.plot(full_night_linspace.to_datetime(), np.ones(length_of_night) * -100.0)
        ax2.set_xlabel("UTC")
        ax2.get_xaxis().set_major_formatter(mdates.DateFormatter("%H:%M"))

        num_ticks = 7
        nn = round(length_of_night / num_ticks)
        ax3_ind = [i * nn for i in range(num_ticks)]

        ax3.plot(full_night_sidereal.to_datetime(), np.ones(length_of_night) * -100.0)
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
#            total_time = df.loc[df.program == k, "dur"].sum()
            label = f"{k}"
            ax1.plot(
                full_night_extra.to_datetime(),
                np.ones(length_of_night) * -100.0,
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

                obs_linspace = st + np.linspace(0, (et - st).to(u.hour).value, num=3) * u.hour
                self._plot_single_moon_dist(
                    ax2, row.ra*u.deg, row.dec*u.deg, times=obs_linspace, linewidth=3.0, color=col
                )

#                self._plot_single_moon_dist(
#                    ax2, row.ra*u.deg, row.dec*u.deg, times=full_night_linspace, linewidth=3.0, color=col, alpha=0.1
#                )

            else:
                self._plot_single_moon_dist(
                    ax2,
                    row.ra*u.deg,
                    row.dec*u.deg,
                    times=full_night_linspace,
                    linewidth=3.0,
                    color=col,
                    linestyle="dotted",
                    alpha=0.3,
                )

        ax2.axhline(30,color='red',lw=2,linestyle='--')
        ax2.axhline(15,color='red',lw=2,linestyle='-')
        leg = ax1.legend(bbox_to_anchor=(1.01, 1.015), loc="upper left", ncol=1, prop={"size": 8})

        # set the linewidth of each legend object
        for legobj in leg.legend_handles:
        #for legobj in leg.legendHandles:
            legobj.set_linewidth(3.0)

        return fig
