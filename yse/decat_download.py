import sys
import pandas as pd
import numpy as np
import requests
import json
import shutil
import os
from astropy.utils.data import download_file
from astropy.table import Table, unique, vstack
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
import urllib.error
import warnings

warnings.filterwarnings('ignore')

natroot='https://astroarchive.noirlab.edu'
adsurl = f'{natroot}/api/adv_search'

def parse_coord(ra, dec):
    coord = None
    if ':' in str(ra) and ':' in str(dec):
        coord = SkyCoord(ra, dec, unit=(u.hour, u.deg), frame='icrs')
    else:
        coord = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')
    return(coord)

class decat_download(object):
    def __init__(self):

        natroot='https://astroarchive.noirlab.edu'
        self.adsurl = f'{natroot}/api/adv_search'

        self.field_table = None

        self.filters={'u':'u DECam c0006 3500.0 1000.0',
                       'g':'g DECam SDSS c0001 4720.0 1520.0',
                       'r':'r DECam SDSS c0002 6415.0 1480.0',
                       'i':'i DECam SDSS c0003 7835.0 1470.0',
                       'z':'z DECam SDSS c0004 9260.0 1520.0',
                       'y':'Y DECam c0005 10095.0 1130.0',
                       'VR':'VR DECam c0007 6300.0 2600.0'}

        self.columns = {
                "ra_center": float,
                "dec_center": float,
                "ra_min": float,
                "ra_max": float,
                "dec_min": float,
                "dec_max": float,
                "md5sum": str,
                "archive_filename": str,
                "instrument": str,
                "telescope": str,
                "proc_type": str,
                "prod_type": str,
                "obs_type": str,
                "release_date": str,
                "proposal": str,
                "url": str,
                "dateobs_center": str,
                "ifilter": str,
            }

    def ra_dec_box_search(self, ra, dec, size=2.5, radius=1.11*u.deg,
        instrument=['decam'], proc_type=['resampled'], prod_type=['image'],
        row_limit=10000):
        coord = parse_coord(ra, dec)
        ra = coord.ra.degree
        dec = coord.dec.degree

        apiurl = f'{self.adsurl}/fasearch/?limit={row_limit}'

        scale = 1.5

        ramin = coord.ra.degree-scale*size/2.0*1./np.cos(dec*np.pi/180.0)
        ramax = coord.ra.degree+scale*size/2.0*1./np.cos(dec*np.pi/180.0)
        demin = coord.dec.degree-scale*size/2.0
        demax = coord.dec.degree+scale*size/2.0

        if demin<-90.0: demin=-90.0
        if demax>90.0: demax=90.0
        if ramin<0: ramin+=360.0
        if ramax>360.0: ramax-=360.0

        jj_base={'outfields': list(self.columns.keys()),'search':[]}
        if instrument: jj_base['search'].append(['instrument']+instrument)
        if proc_type: jj_base['search'].append(['proc_type']+proc_type)
        if prod_type: jj_base['search'].append(['prod_type']+prod_type)

        if ramin>ramax:
            # Need to perform two searches to account for overlap
            jj=jj_base
            jj['search'].append(['ra_center',ramin,360.0])
            jj['search'].append(['dec_center',demin,demax])
            ads_df1 = pd.DataFrame(requests.post(apiurl,json=jj).json()[1:])
            ads_tb1 = Table.from_pandas(ads_df1)

            jj=jj_base
            jj['search'].append(['ra_center',0.0,ramax])
            jj['search'].append(['dec_center',demin,demax])
            ads_df2 = pd.DataFrame(requests.post(apiurl,json=jj).json()[1:])
            ads_tb2 = Table.from_pandas(ads_df2)

            ads_tb = vstack([ads_tb1,ads_tb2])

        else:
            jj=jj_base
            jj['search'].append(['ra_center',ramin,ramax])
            jj['search'].append(['dec_center',demin,demax])

            ads_df = pd.DataFrame(requests.post(apiurl,json=jj).json()[1:])
            ads_tb = Table.from_pandas(ads_df)

        if radius:
            separation = coord.separation(SkyCoord(ads_tb['ra_center'],
                ads_tb['dec_center'], unit=(u.deg, u.deg), frame='icrs'))
            mask = separation < radius
            ads_tb=ads_tb[mask]

        # Set data types
        for key in ads_tb.keys():
            ads_tb[key]=ads_tb[key].astype(self.columns[key])

        return(ads_tb)

    def parse_fieldmaps(self, fieldfile):
        if not os.path.exists(fieldfile):
            w=f'WARNING: {fieldfile} does not exist!'
            print(w)
            return(None)

        self.field_table = Table.read(fieldfile, format='ascii')
        return(self.field_table)

    def download_file(self, row, outdir='', clobber=False):
        if not os.path.exists(outdir):
            try:
                os.makedirs(outdir)
            except:
                e=f'ERROR: could not make output directory {outdir}'
                print(e)
                return(None)

        # Can directly download file from url
        if 'url' not in row.colnames or 'archive_filename' not in row.colnames:
            w='ERROR: cannot download without url and archive_filename in '+\
                'search_columns'
            print(w)
            return(None)
        else:
            ref_file = outdir+os.path.basename(row['archive_filename'])
            url = row['url']
            if os.path.exists(ref_file) and not clobber:
                print(f'WARNING: {ref_file} exists and clobber=False')
                return(ref_file)

            try:
                dat = download_file(url, cache=False, show_progress=True,
                    timeout=30)
            except urllib.error.HTTPError as e:
                msg = str(e)
                print(f'WARNING: for file {ref_file}: {msg}')
                return(None)
            except urllib.error.URLError as e:
                msg = str(e)
                print(f'ERROR: for file {ref_file}: {msg}')
                return(None)

            return(shutil.move(dat, ref_file))

    def download_all_files(self, field_table, after='', before='',
        filts=['u','g','r','i','z','y'], outdir='', by_date=False,
        verbose=False):

        all_data = None
        for row in field_table:
            m='Checking available data for field {0}'.format(row['YSEID'])
            if verbose: print(m)
            ads_tb = self.ra_dec_box_search(row['RA'], row['Dec'])
            if after:
                t = Time(after)
                mask = Time(ads_tb['dateobs_center'])>t
                ads_tb = ads_tb[mask]
            if before:
                t = Time(before)
                mask = Time(ads_tb['dateobs_center'])<t
                ads_tb = ads_tb[mask]
            if filts:
                filt_list = [self.filters[f] for f in filts
                    if f in self.filters.keys()]
                mask = np.array([r['ifilter'] in filt_list for r in ads_tb])
                ads_tb = ads_tb[mask]
            if all_data:
                all_data = vstack([all_data, ads_tb])
            else:
                all_data = ads_tb

        # Get unique rows
        all_data = unique(all_data)
        all_data.sort('dateobs_center')
        print(all_data['dateobs_center','proc_type','archive_filename'])

        if all_data and len(all_data)>0:
            all_data = unique(all_data)
            m='Need to download {0} data files'.format(len(all_data))
            if verbose: print(m)
            for row in all_data:
                self.download_file(row, outdir=outdir)



de = decat_download()
field_table = de.parse_fieldmaps('fieldfile.txt')
de.download_all_files(field_table, after='2021-02-01', outdir='test/',
    verbose=True)
