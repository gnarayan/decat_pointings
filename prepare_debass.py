import numpy as np
import pandas as pd
import readobslogs as ro
import makeobservabilityplot as mop
from datetime import datetime
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import matplotlib.pyplot as plt
import os
import editjson as ej
import sys
import shutil

def split(word):
    return [char for char in word]


df = pd.read_csv('debass_sample.csv')


from datetime import datetime
currdate=datetime.today().strftime('%Y-%m-%d')
datestr = input(f'Please Enter Observing Date in the following format (YYYY-MM-DD)\ndefault {currdate}\n')
if datestr == '': datestr=currdate
basename = 'debass'+datestr[2:].replace('-','')
newjsonfolder = 'json_files/2025A/'+basename
os.makedirs(newjsonfolder, exist_ok=True)
obsdict = ro.run(verbose=True)

doonlysn=False
templates=False
try:
    startsn = sys.argv[1]
    if startsn == '--templates':
        templates=True
    if '--sn=' in startsn:
        onlysn=startsn.replace('--sn=','')
        doonlysn=True
        keepskipping=False
    else:
        if '20' in startsn:
            keepskipping = True
        else:
            keepskipping = False
except:
    keepskipping = False
#skiprows = 60
skiprows = 0
for i,row in df.iterrows():
    if doonlysn:
        if row['snid'] == onlysn:
            pass
        else:
            continue
    if keepskipping:
        if row['snid'] == startsn:
            keepskipping = False
        else:
            continue
    if i < skiprows: continue
    if str(row['snid']) == 'nan': continue
    if str(row['snid']) == 'NaN': continue
    if not templates:
        if 'FINISHED' in row['Following?']:
            os.system('rm jsons/2020B-0053_DEBASS_Brout/EVERYTHING/%s_P*.json'%(row['snid']))
            continue
    if 'ABANDON' in row['Following?']:
        os.system('rm jsons/2020B-0053_DEBASS_Brout/EVERYTHING/%s_P*.json'%(row['snid']))
        continue
    if 'LOST' in row['Following?']:
        os.system('rm jsons/2020B-0053_DEBASS_Brout/EVERYTHING/%s_P*.json'%(row['snid']))
        continue
    if 'NON IA YSE' in row['Following?']:  continue
    if 'YSE' in row['Following?']:  continue
    if '91' in row['TNS class']:  continue
    if ('ia' in row['TNS class'].lower()) | ('?' in row['TNS class']):
        #asdf
        os.system('clear')
        try:
            print(obsdict[row['snid']])
        except:
            print('NO OBSERVATIONS YET')

        print('SPEC CLASS: %s'%row['TNS class'])
        print('WiFeS Spectrum of transient: %s'%row['WiFeS Spectrum of transient'])
        print('Redshift: %s'%row['Redshift'])
        print()
        print('Comment: %s'%row['Comment'])
        print()
        print('Got Template? %s'%row['Got Template?'])
        print()
        print('WiFeS Time Series %s'%row['Following up with time-series spectroscopy'])
        print()

        if (row['YSE Field'] != '') & (row['YSE Field'] != '?') & (str(row['YSE Field'])!='nan'): print('THIS IS SHARED WITH YSE: %s'%row['YSE Field'])
        mop.doplot(datestr,ra=float(row['RA']),dec=float(row['DEC']),name=row['snid'],block=False)
        priority = input('\nPlease enter a priority for this object (1 2 3 TCTM)\n')
        if not priority in ['1','2','3','TCTM']:
            priority = input('Please enter a valid priority for this object (1 2 3 TCTM)\n')
            if not priority in ['1','2','3','TCTM']:
                priority = input('Please enter a valid priority for this object (1 2 3 TCTM)\n')
        filters = input('Please enter filters for this object (default: griz)\n')
        if filters == '':
            filters = ['g','r','i','z']
        else:
            filters = split(filters)
        #try:
        default_exptimes = ej.getfiltersexptimes('jsons/2020B-0053_DEBASS_Brout/TEMPLATE/%s.json'%row['snid'])
        #    maketemplate = False
        #except:
        #    print('could not find template, using generic.json. This is generally okay.')
        #    default_exptimes = ej.getfiltersexptimes('jsons/2020B-0053_DEBASS_Brout/TEMPLATE/generic.json')
        #    maketemplate = True
            
        exptimes = []
        for f in filters:
            if not f in default_exptimes.keys():
                if templates:
                    default_exptimes[f] = 30
                else:
                    default_exptimes[f] = 15
            if templates:
                exptime = input(f'Enter Exptime {f} (default 30.0)\n')
            else:
                exptime = input(f'Enter Exptime {f} (default {default_exptimes[f]})\n')
            if exptime == '':
                if templates:
                    exptime = 30.0
                else:
                    exptime = str(default_exptimes[f])
            exptimes.append(exptime)
        plt.clf()
        #if maketemplate:
        #    ej.edit('jsons/2020B-0053_DEBASS_Brout/TEMPLATE/generic.json',priority,filters,exptimes,
        #            'jsons/2020B-0053_DEBASS_Brout/TEMPLATE/%s.json'%(row['snid']))
        if os.path.exists('jsons/2020B-0053_DEBASS_Brout/TEMPLATE/%s.json'%row['snid']):
            os.system('rm jsons/2020B-0053_DEBASS_Brout/EVERYTHING/%s_P*.json'%(row['snid']))
            ej.edit('jsons/2020B-0053_DEBASS_Brout/TEMPLATE/%s.json'%row['snid'],priority,filters,exptimes,
                 'jsons/2020B-0053_DEBASS_Brout/EVERYTHING/%s_P%s.json'%(row['snid'],priority))
        else:
            print('WARNING: jsons/2020B-0053_DEBASS_Brout/TEMPLATE/%s.json'%row['snid']+'\nDoes Not Exist')
            input('press enter to continue')

shutil.copytree('jsons/2020B-0053_DEBASS_Brout/EVERYTHING', newjsonfolder, dirs_exist_ok=True)
