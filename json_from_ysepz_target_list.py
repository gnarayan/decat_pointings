import pandas as pd
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from glob import glob
import get_silicon as gs
import os
import make_json as mj
import ccdmap
import sys
import random

def get_desired_ccd(x):
    if x in ccdmap.ccdmap.keys():
        return ccdmap.ccdmap[x]
    else:
        return get_random_ccd(x)
        
def get_random_ccd(x):
    print(x,'not found in ccdmap.py')
    ccdlist = [1,3,4,5,6,7,8,9,10,\
               11,12,13,14,15,16,17,18,19,20,\
               21,22,23,24,25,26,27,28,29,30,\
               32,33,34,35,36,37,38,39,40,\
               41,42,43,44,45,46,47,48,49,50,\
               51,52,53,54,55,56,57,58,59,60,62]
    ccd = random.choice(ccdlist)
    print('putting ',x,'on chip',ccd)
    lines = open('ccdmap.py','r').readlines()
    f = open('ccdmap.py','w')
    for l in lines:
        if not '}' in l:
            f.write(l)
        else:
            f.write("    '%s': %d,\n}"%(x,ccd))
    f.close()
    return ccd
    
        
    
def get_desired_observations(x):
    return 'g15,r15,i15,z15'
    #obs = x.split(';')[1].split(';')[0]
    #return  obs

def get_priorities(x):
    return '1'
    #priority = x.split('PRIORITY')[1][0]
    #return  priority

def get_propid(x):
    return '2020B-0053'
    #propid = x.split('{')[1].split('}')[0]
    #return  propid

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
alreadyprinted = []
for infile in infiles:
    date = infile.split('/')[-1].split('_')[1].split('.txt')[0].replace('-','')
    data = open(infile,'r').readlines()
    outfilep = reformatdir+'/'+infile.split('/')[-1]
    outfile = open(outfilep,'w')

    json_outpath = 'jsons/2020B-0053_DEBASS_Brout/TEMPLATE/'
    if not os.path.exists(json_outpath):
        os.mkdir(json_outpath)

    outfile.write('name ra_h ra_m ra_s dec_d dec_m dec_s equinox a b mag c d comment\n')
    for line in data[1:]:
        if line[0] == '#': continue
        outfile.write(line)
    outfile.close()
    df = pd.read_csv(outfilep,delim_whitespace=True, error_bad_lines=False)
    df = parse_infile(df)
    
    fieldra,fielddec = gs.get_field_center_for_target_on_specific_ccd(df['candRA'],df['candDEC'],df['ccd'])
    
    df['pointRA'] = fieldra
    df['pointDEC'] = fielddec

    df['expTypes'] = 'object'
    df['programs'] = 'DECAT'
    
    df[['name','candRA','candDEC','ccd','obs','priority','pointRA','pointDEC']].to_csv(reformatdir+date+'.txt',sep=' ',index=False)
    df[['name','candRA','candDEC']].to_csv(reformatdir+date+'_for_iObserve.txt',index=False,header=False,sep=' ')

    #for i,row in df.iterrows():
    #    print(row['name'],'Candidate RA',row['candRA'],'Candidate DEC',row['candDEC'],'Field RA',row['pointRA'],'Field DEC',row['pointDEC'])
    #for i,row in df.iterrows():
    #    print("OR (power(power(t.ra - %s,2)+power(t.dec - %s,2),.5)<2 AND t.name != '%s')"%(row['candRA'],row['candDEC'],row['name']))

    #print('-'*100)
    for i,row in df.iterrows():
        if not row['name'] in alreadyprinted:
            print("OR (power(power(t.ra - %s,2)+power(t.dec - %s,2),.5)<.1 AND t.name != '%s')"%(row['candRA'],row['candDEC'],row['name']))
            alreadyprinted.append(row['name'])
    mj.individual(json_outpath,df['name'],df['pointRA'],df['pointDEC'],df['obs'],df['propid'],df['name']+'_P'+df['priority'].astype(str),df['expTypes'],df['programs'])
    
