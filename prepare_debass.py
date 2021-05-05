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

def split(word):
    return [char for char in word]


df = pd.read_csv('debass_sample.csv')


from datetime import datetime
currdate=datetime.today().strftime('%Y-%m-%d')
datestr = input(f'Please Enter Observing Date in the following format (YYYY-MM-DD)\ndefault {currdate}\n')
if datestr == '': datestr=currdate
obsdict = ro.run(verbose=False)

try:
    startsn = sys.argv[1]
    keepskipping = True
except:
    keepskipping = False
#skiprows = 60
skiprows = 0
for i,row in df.iterrows():
    if keepskipping:
        if row['snid'] == startsn:
            keepskipping = False
        else:
            continue
    if i < skiprows: continue
    if 'FINISHED' in row['Following?']: continue
    if 'ABANDON' in row['Following?']:	continue
    if 'LOST' in row['Following?']:  continue
    if 'NON IA YSE' in row['Following?']:  continue
    if 'YSE' in row['Following?']:  continue
    if '91' in row['TNS class']:  continue
    if ('ia' in row['TNS class'].lower()) | ('?' in row['TNS class']):
        os.system('clear')
        try:
            print(obsdict[row['snid']])
        except:
            print('NO OBSERVATIONS YET')
        print('SPEC CLASS: %s'%row['TNS class'])
        print('Redshift: %s'%row['Redshift'])
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
        default_exptimes = ej.getfiltersexptimes('jsons/2020B-0053_DEBASS_Brout/BASE/%s.json'%row['snid'])
        exptimes = []
        for f in filters:
            if not f in default_exptimes.keys():
                default_exptimes[f] = 15
            exptime = input(f'Enter Exptime {f} (default {default_exptimes[f]})\n')
            if exptime == '': exptime = str(default_exptimes[f])
            exptimes.append(exptime)
        plt.clf()
        if os.path.exists('jsons/2020B-0053_DEBASS_Brout/BASE/%s.json'%row['snid']):
            os.system('rm jsons/2020B-0053_DEBASS_Brout/EVERYTHING/%s_P*.json'%(row['snid']))
            ej.edit('jsons/2020B-0053_DEBASS_Brout/BASE/%s.json'%row['snid'],priority,filters,exptimes,
                 'jsons/2020B-0053_DEBASS_Brout/EVERYTHING/%s_P%s.json'%(row['snid'],priority))
        else:
            print('WARNING: jsons/2020B-0053_DEBASS_Brout/BASE/%s.json'%row['snid']+'\nDoes Not Exist')
            input('press enter to continue')
