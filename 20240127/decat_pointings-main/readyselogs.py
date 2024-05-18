#!/usr/bin/env python
# D. Jones - 4/12/21

import numpy as np
import pylab as plt
import os
import pathlib
import astropy.table as at
import glob
import datetime
import dateutil.parser
import pandas as pd

from astropy.time import Time

def date_to_mjd(date):
    time = Time(date,scale='utc')
    return time.mjd

def mjd_to_date(mjd):
    time = Time(mjd,format='mjd',scale='utc')
    return time.isot

def parse_qcinv_file(filename):

    qcinv = at.Table(names=(
        'expid','ra','dec','ut','fil','time',
        'secz','psf','sky','cloud','teff','Object'),
                     dtype=(float,float,float,'S20','S20','S20',float,float,float,float,float,'S20'))

    with open(filename) as fin:
        for line in fin:
            if line.startswith('#'): continue
            lineparts = line.split()
            cols = np.append(lineparts[0:11],' '.join(lineparts[11:]))
            qcinv.add_row(cols)

    return qcinv

def parse_qcinv_dillon(f):

    datestr = f.split('/')[-1].split('.')[0]
    date = datestr[:4]+'-'+datestr[4:6]+'-'+datestr[6:8]
    expnums = []
    objects = []
    filts = []
    ras = []
    decs = []
    teffs = []
    print(f)
    for l in open(f,'r').readlines():
        if l[0] == '#': continue
        if len(l.split()) < 2: continue
        if l.split()[0] == 'MJD': continue
        #if l.split()[0] == 'ID': continue
        expnums.append(int(l.split()[0]))
        objects.append(str(l.split()[-1]))
        if len(l.split()) <= 10:
            teffs.append(np.nan)
        else:
            teffs.append(float(l.split()[10]))
        filts.append(str(l.split()[4]))
        ras.append(float(l.split()[1]))
        decs.append(float(l.split()[2]))

    dates = [date for e in range(len(expnums))]
    tdf = pd.DataFrame.from_dict({'expnum':expnums,'object':objects,'date':dates,'filt':filts,'ra':ras,'dec':decs,'teff':teffs})

    return tdf
    
def main():
    files =  glob.glob('2021A/*/*nv')
    ysefieldfiles = glob.glob('2021A/*/decat_YSE_list_*.txt')

    fielddates = [date_to_mjd(dateutil.parser.parse(f.split('_')[-1].split('.')[0])) for \
                  f in ysefieldfiles]
    ysefieldfile = np.array(ysefieldfiles)[np.argsort(fielddates)][::-1][0]

    data = at.Table.read(ysefieldfile,format='ascii')
    fieldids = np.loadtxt('fields.txt',unpack=True,usecols=[0],dtype=str)

    ysedict = {}
    mjddict = {}
    # let's list all the YSE obs
    for f in files:
        qc = parse_qcinv_dillon(f)

        for y in data['ID']:
            if y in qc['object'].values:
                if y in ysedict.keys():
                    ysedict[y] += [f"{qc['date'][0]}: {','.join(qc['filt'][qc['object'] == y])}"]
                    mjddict[y] += [date_to_mjd(qc['date'][0])]
                else:
                    ysedict[y] = [f"{qc['date'][0]}: {','.join(qc['filt'][qc['object'] == y])}"]
                    mjddict[y] = [date_to_mjd(qc['date'][0])]
        for y in fieldids:
            if y in qc['object'].values:
                if y in ysedict.keys():
                    ysedict[y] += [f"{qc['date'][0]}: {','.join(qc['filt'][qc['object'] == y])}"]
                    mjddict[y] += [date_to_mjd(qc['date'][0])]
                else:
                    ysedict[y] = [f"{qc['date'][0]}: {','.join(qc['filt'][qc['object'] == y])}"]
                    mjddict[y] = [date_to_mjd(qc['date'][0])]
                    
    #now lets get the time since first obs
    datedict = {}
    
    for y in ysedict.keys():
        print('')
        try: print(f"{y} {data['SNID'][data['ID'] == y][0]}")
        except: print(f"{y} None")
        print(f"time since first obs: {date_to_mjd(datetime.datetime.utcnow())-np.min(mjddict[y]):.0f} days")
        print(f"drop date: {mjd_to_date(np.min(mjddict[y])+40).split('T')[0]}")
        print(f"days until drop: {(np.min(mjddict[y])+40)-date_to_mjd(datetime.datetime.utcnow()):.0f}")
        print('---------')
        for e in np.array(ysedict[y])[np.argsort(mjddict[y])]:
            print(e)
                    
if __name__ == "__main__":
    main()
