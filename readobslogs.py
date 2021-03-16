import pandas as pd
from glob import glob
import ccdmap as c
import numpy as np
import os

snids = list(c.ccdmap.keys())
m = pd.read_csv('fieldmaps.txt',delim_whitespace=True)
snids.extend(list(m['SNID'].to_numpy()))

os.system('rm obslogs/*~')

ysedict = {}
ysedict2 = {}
rysedict = {}
for row in m.iterrows():
    #print(row[1]['YSEID'])
    ysedict[str(row[1]['YSEID'])] = str(row[1]['SNID'])
    if str(row[1]['YSEID'])[-1] == 'a':
        ysedict2[str(row[1]['YSEID'])[:-1]+'b'] = str(row[1]['SNID'])
    if str(row[1]['YSEID'])[-1] == 'b':
        ysedict2[str(row[1]['YSEID'])[:-1]+'a'] = str(row[1]['SNID'])
    rysedict[str(row[1]['SNID'])] = str(row[1]['YSEID'])


    
dfs = []
for f in glob('obslogs/*'):
    date = f.split('/')[-1].split('.')[0]
    tdf = pd.read_csv(f,names=['expnum','ra','dec','ut','filt','exp','secz','type','object'],delim_whitespace=True,comment='#')
    tdf['date'] = [date for e in range(len(tdf))]
    dfs.append(tdf)

obsdict = {}
df = pd.concat(dfs)
for row in df.iterrows():
    r = row[1]
    if str(r['object']).split('_')[0] in snids:
        snid = r['object'].split('_')[0]
        if not snid in obsdict.keys():
            obsdict[snid] = {'dates':[],'expnums':[],'filts':[]}
        #print(snid,r['date'])
        obsdict[snid]['dates'].append(r['date'])
        obsdict[snid]['expnums'].append(r['expnum'])
        obsdict[snid]['filts'].append(r['filt'])
    elif str(r['object']) in ysedict.keys():
        snid = ysedict[r['object']]
        if not snid in obsdict.keys():
            obsdict[snid] = {'dates':[],'expnums':[],'filts':[]}
        obsdict[snid]['dates'].append(r['date'])
        obsdict[snid]['expnums'].append(r['expnum'])
        obsdict[snid]['filts'].append(r['filt'])
    elif str(r['object']) in ysedict2.keys():
        snid = ysedict2[r['object']]
        if not snid in obsdict.keys():
            obsdict[snid] = {'dates':[],'expnums':[],'filts':[]}
        obsdict[snid]['dates'].append(r['date'])
        obsdict[snid]['expnums'].append(r['expnum'])
        obsdict[snid]['filts'].append(r['filt'])

print('-'*25)
for k,v in obsdict.items():
    if k in rysedict.keys():
        print('SNID',k,'YSE',rysedict[k])
    else:
        print('SNID',k)

    for date in np.sort(np.unique(v['dates'])):
        ww = np.array(v['dates']) == date
        filts = np.array(v['filts'])[ww]
        print(date+':',filts,)
    print()
    #print('Dates',v['dates'])
    print('-'*25)
