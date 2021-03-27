import pandas as pd
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from glob import glob
import get_silicon as gs
import os
import make_json as mj
import ccdmap

def get_desired_ccd(x):
    return ccdmap.ccdmap[x]

def get_desired_observations(x):
    obs = x.split(';')[1].split(';')[0]
    return  obs

def get_priorities(x):
    priority = x.split('PRIORITY')[1][0]
    return  priority

def get_propid(x):
    propid = x.split('{')[1].split('}')[0]
    return  propid

def parse_infile(df):
    df['coords'] = df['ra_h'].astype(str)+' '+df['ra_m'].astype(str)+' '+df['ra_s'].astype(str)+' '+df['dec_d'].astype(str)+' '+df['dec_m'].astype(str)+' '+df['dec_s'].astype(str)
    coords = SkyCoord(df['coords'].to_numpy(), unit=(u.hourangle, u.deg))
    df['candRA'] = coords.ra.degree
    df['candDEC'] = coords.dec.degree
    df['ccd'] = df['name'].apply(get_desired_ccd)
    df['obs'] = df['comment'].apply(get_desired_observations)
    df['priority'] = df['comment'].apply(get_priorities)
    df['propid'] = df['comment'].apply(get_propid)

    return df
    
    
reformatdir = 'reformatted_target_downloads/'
indir = 'ysepz_downloads/'
infiles = glob(indir+'/Blanco*.txt')

for infile in infiles:
    date = infile.split('/')[-1].split('_')[1].split('.txt')[0].replace('-','')
    data = open(infile,'r').readlines()
    outfilep = reformatdir+'/'+infile.split('/')[-1]
    outfile = open(outfilep,'w')

    json_outpath = 'jsons/%s/'%date
    if not os.path.exists(json_outpath):
        os.mkdir(json_outpath)

    outfile.write('name ra_h ra_m ra_s dec_d dec_m dec_s equinox a b mag c d comment\n')
    for line in data[1:]:
        if line[0] == '#': continue
        outfile.write(line)
    outfile.close()

    df = pd.read_csv(outfilep,delim_whitespace=True)
    df = parse_infile(df)
    
    fieldra,fielddec = gs.get_field_center_for_target_on_specific_ccd(df['candRA'],df['candDEC'],df['ccd'])
    
    df['pointRA'] = fieldra
    df['pointDEC'] = fielddec

    df['expTypes'] = 'object'
    df['programs'] = 'DECAT'
    
    df[['name','candRA','candDEC','ccd','obs','priority','pointRA','pointDEC']].to_csv(reformatdir+date+'.txt',sep=' ',index=False)
    df[['name','candRA','candDEC']].to_csv(reformatdir+date+'_for_iObserve.txt',index=False,header=False,sep=' ')
    for i,row in df.iterrows():
        print(row['name'],'Candidate RA',row['candRA'],'Candidate DEC',row['candDEC'],'Field RA',row['pointRA'],'Field DEC',row['pointDEC'])
    mj.individual(json_outpath,df['name']+'_'+df['candRA'].round(0).astype(str)+'_'+df['candDEC'].round(0).astype(str)+'_PRIORITY'+df['priority'].astype(str),df['pointRA'],df['pointDEC'],df['obs'],df['propid'],df['name']+'_P'+df['priority'].astype(str),df['expTypes'],df['programs'])

