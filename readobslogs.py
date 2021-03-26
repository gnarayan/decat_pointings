import pandas as pd
from glob import glob
import ccdmap as c
import numpy as np
import os

snids = list(c.ccdmap.keys())
m = pd.read_csv('fieldmaps.txt',delim_whitespace=True)
snids.extend(list(m['SNID'].to_numpy()))

ignoref = open('ignore.list','r').readlines()
ignore = [i.strip() for i in ignoref]


os.system('rm obslogs/*~')
os.system('rm 2021A/*/*~')

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
for f in glob('2021A/*/*nv'):
    date = f.split('/')[-1].split('.')[0]
    expnums = []
    objects = []
    filts = []
    ras = []
    decs = []
    print(f)
    for l in open(f,'r').readlines():
        if l[0] == '#': continue
        if l.split()[0] == 'MJD': continue
        #if l.split()[0] == 'ID': continue
        expnums.append(int(l.split()[0]))
        objects.append(str(l.split()[-1]))
        filts.append(str(l.split()[4]))
        ras.append(float(l.split()[1]))
        decs.append(float(l.split()[2]))
        #tdf = pd.read_csv(f,names=['expnum','ra','dec','ut','filt','exp','secz','type','object'],delim_whitespace=True,comment='#')
    dates = [date for e in range(len(expnums))]
    tdf = pd.DataFrame.from_dict({'expnum':expnums,'object':objects,'date':dates,'filt':filts,'ra':ras,'dec':decs})
    dfs.append(tdf)

obsdict = {}
df = pd.concat(dfs)
for row in df.iterrows():
    r = row[1]
    if str(r['object']).split('_')[0] in snids:
        snid = r['object'].split('_')[0]
        ra = r['ra']
        dec = r['dec']
        if not snid in obsdict.keys():
            obsdict[snid] = {'dates':[],'expnums':[],'filts':[],'ra':ra,'dec':dec}
        #print(snid,r['date'])
        obsdict[snid]['dates'].append(r['date'])
        obsdict[snid]['expnums'].append(r['expnum'])
        obsdict[snid]['filts'].append(r['filt'])
    elif str(r['object']) in ysedict.keys():
        snid = ysedict[r['object']]
        ra = r['ra']
        dec = r['dec']
        if not snid in obsdict.keys():
            obsdict[snid] = {'dates':[],'expnums':[],'filts':[],'ra':ra,'dec':dec}
        obsdict[snid]['dates'].append(r['date'])
        obsdict[snid]['expnums'].append(r['expnum'])
        obsdict[snid]['filts'].append(r['filt'])
    elif str(r['object']) in ysedict2.keys():
        snid = ysedict2[r['object']]
        ra = r['ra']
        dec = r['dec']
        if not snid in obsdict.keys():
            obsdict[snid] = {'dates':[],'expnums':[],'filts':[],'ra':ra,'dec':dec}
        obsdict[snid]['dates'].append(r['date'])
        obsdict[snid]['expnums'].append(r['expnum'])
        obsdict[snid]['filts'].append(r['filt'])

print('-'*25)
cnt = 0
for k,v in obsdict.items():
    if k in ignore: continue
    if k in rysedict.keys():
        print('SNID',k,'YSE',rysedict[k],'RA',v['ra'],'DEC',v['dec'])
    else:
        print('SNID',k,'RA',v['ra'],'DEC',v['dec'])

    cnt += 1
    for date in np.sort(np.unique(v['dates'])):
        ww = np.array(v['dates']) == date
        filts = np.array(v['filts'])[ww]
        print(date+':',filts,)
    print()
    #print('Dates',v['dates'])
    print('-'*25)
print('Total SNe %d'%cnt)
