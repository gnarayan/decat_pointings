import numpy as np
import pandas as pd
from copy import copy

c1_corners = pd.read_csv('cornerCoords_SN-C1.dat',delim_whitespace=True)
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


def get_field_center_for_target_on_specific_ccd(ras,decs,ccds,extra_offset_in_arcsec=50):
    field_cen_ras = []
    field_cen_decs = []
    for ra,dec,ccd in zip(ras,decs,ccds):
        field_cen_ras.append(ra+raoff[int(ccd)]+(np.random.randint(2)-.5)*2*extra_offset_in_arcsec*0.000277778)
        field_cen_decs.append(dec+decoff[int(ccd)]+(np.random.randint(2)-.5)*2*extra_offset_in_arcsec*0.000277778)
    return np.array(field_cen_ras),np.array(field_cen_decs)
        

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

