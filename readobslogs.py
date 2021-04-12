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

outtxt = open('debass_sne.txt','w')

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
        #tdf = pd.read_csv(f,names=['expnum','ra','dec','ut','filt','exp','secz','type','object'],delim_whitespace=True,comment='#')
    dates = [date for e in range(len(expnums))]
    tdf = pd.DataFrame.from_dict({'expnum':expnums,'object':objects,'date':dates,'filt':filts,'ra':ras,'dec':decs,'teff':teffs})
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
            obsdict[snid] = {'dates':[],'expnums':[],'filts':[],'ra':ra,'dec':dec,'teffs':[]}
        #print(snid,r['date'])
        obsdict[snid]['dates'].append(r['date'])
        obsdict[snid]['expnums'].append(r['expnum'])
        obsdict[snid]['filts'].append(r['filt'])
        obsdict[snid]['teffs'].append(r['teff'])
    elif str(r['object']) in ysedict.keys():
        snid = ysedict[r['object']]
        ra = r['ra']
        dec = r['dec']
        if not snid in obsdict.keys():
            obsdict[snid] = {'dates':[],'expnums':[],'filts':[],'ra':ra,'dec':dec,'teffs':[]}
        obsdict[snid]['dates'].append(r['date'])
        obsdict[snid]['expnums'].append(r['expnum'])
        obsdict[snid]['filts'].append(r['filt'])
        obsdict[snid]['teffs'].append(r['teff'])
    elif str(r['object']) in ysedict2.keys():
        snid = ysedict2[r['object']]
        ra = r['ra']
        dec = r['dec']
        if not snid in obsdict.keys():
            obsdict[snid] = {'dates':[],'expnums':[],'filts':[],'ra':ra,'dec':dec,'teffs':[]}
        obsdict[snid]['dates'].append(r['date'])
        obsdict[snid]['expnums'].append(r['expnum'])
        obsdict[snid]['filts'].append(r['filt'])
        obsdict[snid]['teffs'].append(r['teff'])

print('-'*25)
cnt = 0
for k,v in obsdict.items():
    if k in ignore: continue
    if k in rysedict.keys():
        print('SNID',k,'YSE_Field',rysedict[k],'CCD',c.ccdmap[k],'RA',v['ra'],'DEC',v['dec'])
        outtxt.write(' '.join(['SNID',str(k),'YSE_Field',rysedict[k],'CCD',str(c.ccdmap[k]),'\n']))
    else:
        print('SNID',k,'CCD',c.ccdmap[k],'RA',v['ra'],'DEC',v['dec'])
        outtxt.write(' '.join(['SNID',str(k),'CCD',str(c.ccdmap[k]),'\n']))
    #cnt += 1
    if len(np.unique(v['dates'])) > 1:
        cnt += 1
    for date in np.sort(np.unique(v['dates'])):
        ww = np.array(v['dates']) == date
        filts = np.array(v['filts'])[ww]
        print(date+':',filts,'Avg Teff %.2f'%np.nanmean(np.array(v['teffs'])[ww]),)
        outstr = date+':'+str(filts)+'\n'
        outtxt.write(outstr)

    outtxt.write('\n')
    print()
    #print('Dates',v['dates'])
    print('-'*25)
    outtxt.write('-'*25)
    outtxt.write('\n')
outtxt.close()
print('Total SNe %d'%cnt)

print('-'*1000)
allexps = []
for k,v in obsdict.items():
    if k in ignore: continue
    #print(v['expnums'])
    allexps.extend(v['expnums'])
allexps = np.unique(allexps)

fout = open('debass_allexpnums.txt','w')
for exp in allexps:
    fout.write(str(exp)+'\n')
fout.close()

    
