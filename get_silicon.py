import numpy as np
import pandas as pd
from copy import copy
from pathlib import Path

_CORNER_FILE = Path(__file__).resolve().with_name('cornerCoords_SN-C1.dat')
c1_corners = pd.read_csv(_CORNER_FILE, sep=r'\s+')
c1_field_center_ra = 54.2743
c1_field_center_dec = -27.1116

def get_ccd_deg_offsets(df):
    ccd_ra_centers = {}
    ccd_dec_centers = {}
    for ccd in range(1,63):
        ww = df['CCD'] == ccd
        ccd_ra_centers[ccd] = -1*(np.mean(df['RA'][ww].to_numpy()) - c1_field_center_ra)
        ccd_dec_centers[ccd] = -1*(np.mean(df['Dec'][ww].to_numpy()) - c1_field_center_dec)
    return ccd_ra_centers,ccd_dec_centers

raoff,decoff = get_ccd_deg_offsets(c1_corners)


def get_field_center_for_target_on_specific_ccd(
    ras,
    decs,
    ccds,
    extra_offset_in_arcsec=50,
    seed=None,
):
    """Return field centers that place each target near its requested CCD center.

    The extra displacement is sampled uniformly by area inside a circle with
    radius ``extra_offset_in_arcsec``.  The previous implementation selected
    only +/- the requested offset on each axis, which restricted targets to
    four corners of a square.

    Parameters
    ----------
    ras, decs : array-like
        Target sky coordinates in degrees.
    ccds : array-like
        Requested DECam CCD numbers (1--62).
    extra_offset_in_arcsec : float, optional
        Radius of the circular displacement region, in arcseconds.
    seed : int or None, optional
        Seed for reproducible pointings.  With ``None``, fresh random offsets
        are generated on each call.
    """
    ras = np.asarray(ras, dtype=float)
    decs = np.asarray(decs, dtype=float)
    ccds = np.asarray(ccds, dtype=int)

    if not (ras.ndim == decs.ndim == ccds.ndim == 1):
        raise ValueError('ras, decs, and ccds must be one-dimensional')
    if not (len(ras) == len(decs) == len(ccds)):
        raise ValueError('ras, decs, and ccds must have the same length')
    if not np.isfinite(extra_offset_in_arcsec) or extra_offset_in_arcsec < 0:
        raise ValueError('extra_offset_in_arcsec must be a finite non-negative value')

    invalid_ccds = sorted(set(ccds) - set(raoff))
    if invalid_ccds:
        raise ValueError(f'CCD numbers must be between 1 and 62; got {invalid_ccds}')

    nominal_field_ras = ras + np.array([raoff[ccd] for ccd in ccds])
    nominal_field_decs = decs + np.array([decoff[ccd] for ccd in ccds])

    rng = np.random.default_rng(seed)
    radii_arcsec = extra_offset_in_arcsec * np.sqrt(rng.random(len(ras)))
    position_angles = rng.uniform(0.0, 2.0 * np.pi, len(ras))
    east_offsets_arcsec = radii_arcsec * np.cos(position_angles)
    north_offsets_arcsec = radii_arcsec * np.sin(position_angles)

    cos_dec = np.cos(np.deg2rad(nominal_field_decs))
    if np.any(np.isclose(cos_dec, 0.0)):
        raise ValueError('cannot convert an east-west offset at a celestial pole')

    field_cen_ras = nominal_field_ras + east_offsets_arcsec / (3600.0 * cos_dec)
    field_cen_decs = nominal_field_decs + north_offsets_arcsec / 3600.0
    return field_cen_ras, field_cen_decs
        

def is_on_silicon(RApoints,DECpoints,RAcands,DECcands):
    ccds = []
    raminarcsecs = []
    decminarcsecs = []
    for RApoint,DECpoint,RAcand,DECcand in zip(RApoints,DECpoints,RAcands,DECcands):
        found = np.nan
        raminarcsec = np.nan
        decminarcsec = np.nan
        for ccd in range(1,63):
            ww = c1_corners['CCD'] == ccd
            maxra = np.max(c1_corners['RA'][ww].to_numpy())+RApoint-c1_field_center_ra
            minra = np.min(c1_corners['RA'][ww].to_numpy())+RApoint-c1_field_center_ra
            maxdec = np.max(c1_corners['Dec'][ww].to_numpy())+DECpoint-c1_field_center_dec
            mindec = np.min(c1_corners['Dec'][ww].to_numpy())+DECpoint-c1_field_center_dec
            if (RAcand < maxra) & (RAcand > minra) & (DECcand < maxdec) & (DECcand > mindec):
                found = copy(ccd)
                raminarcsec = min([maxra-RAcand,RAcand-minra])*3600
                decminarcsec =	min([maxdec-DECcand,DECcand-mindec])*3600

        ccds.append(found)
        raminarcsecs.append(raminarcsec)
        decminarcsecs.append(decminarcsec)
    return np.array(ccds),np.array(raminarcsecs),np.array(decminarcsecs)
