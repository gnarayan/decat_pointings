import pandas as pd
from glob import glob
import ccdmap as c
import numpy as np
import os

def run(verbose=True):
    snids = list(c.ccdmap.keys())
    m = pd.read_csv('yse/fieldmaps.txt',delim_whitespace=True, comment='#')
    snids.extend(list(m['SNID'].to_numpy()))

    ignoref = open('debass/ignore.list','r').readlines()
    ignore = [i.strip() for i in ignoref]

    outtxt = open('debass/debass_sne.txt','w')
    outtxtexp = open('debass/debass_sne_wexpnums.txt','w')
    
    os.system('rm obslogs/*~ >& dump')
    os.system('rm 2021A/*/*~ >& dump')
    os.system('rm 2022A/*/*~ >& dump')
    os.system('rm 2022B/*/*~ >& dump')
    os.system('rm 2023A/*/*~ >& dump')
    os.system('rm 2023B/*/*~ >& dump')
    os.system('rm 2024A/*/*~ >& dump')
    os.system('rm 2024B/*/*~ >& dump')
    os.system('rm 2025A/*/*~ >& dump')

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
    l = glob('2021A/*/*qc*nv')
    l.extend(glob('2021B/*/*qc*nv'))
    l.extend(glob('2022A/*/*qc*nv'))
    l.extend(glob('2022B/*/*qc*nv'))
    l.extend(glob('2023A/*/*qc*nv'))
    l.extend(glob('2023B/*/*qc*nv'))
    l.extend(glob('2024A/*/*qc*nv'))
    l.extend(glob('2024B/*/*qc*nv'))
    l.extend(glob('2025A/*/*qc*nv'))
    for f in l:
        #print(f)
        datestr = f.split('/')[-1].split('.')[0]
        date = datestr[:4]+'-'+datestr[4:6]+'-'+datestr[6:8]
        datestrf = int(datestr)
        expnums = []
        objects = []
        filts = []
        ras = []
        decs = []
        datestrfs = []
        teffs = []
        if verbose: print(f)
        for l in open(f,'r').readlines():
            if l[0] == '#': continue
            if len(l.split()) < 2: continue
            if l.split()[0] == 'MJD': continue
            #if l.split()[0] == 'ID': continue
            expnums.append(int(l.split()[0]))
            objects.append(str(l.split()[-1]))
            print(str(l.split()[-1]))
            if len(l.split()) <= 10:
                teffs.append(np.nan)
            else:
                teffs.append(float(l.split()[10]))
            filts.append(str(l.split()[4]))
            ras.append(float(l.split()[1]))
            decs.append(float(l.split()[2]))
            #tdf = pd.read_csv(f,names=['expnum','ra','dec','ut','filt','exp','secz','type','object'],delim_whitespace=True,comment='#')
        dates = [date for e in range(len(expnums))]
        datestrfs = [datestrf for e in range(len(expnums))]
        tdf = pd.DataFrame.from_dict({'expnum':expnums,'object':objects,'date':dates,'filt':filts,'ra':ras,'dec':decs,'teff':teffs,'datestrf':datestrfs})
        dfs.append(tdf)

    obsdict = {}
    debassonlyexpnums = []
    df = pd.concat(dfs)
    for row in df.iterrows():
        r = row[1]
        if str(r['object']).split('_')[0] in snids:
            snid = r['object'].split('_')[0]
            ra = r['ra']
            dec = r['dec']
            if not snid in obsdict.keys():
                obsdict[snid] = {'dates':[],'expnums':[],'filts':[],'ra':ra,'dec':dec,'datestrfs':[],'teffs':[], 'objects':[]}
            #print(snid,r['date'])
            obsdict[snid]['dates'].append(r['date'])
            obsdict[snid]['expnums'].append(r['expnum'])
            obsdict[snid]['filts'].append(r['filt'])
            obsdict[snid]['teffs'].append(r['teff'])
            obsdict[snid]['datestrfs'].append(r['datestrf'])
            obsdict[snid]['objects'].append(r['object'])
            debassonlyexpnums.append(r['expnum'])
        elif str(r['object']) in ysedict.keys():
            snid = ysedict[r['object']]
            ra = r['ra']
            dec = r['dec']
            if not snid in obsdict.keys():
                obsdict[snid] = {'dates':[],'expnums':[],'filts':[],'ra':ra,'dec':dec,'datestrfs':[],'teffs':[], 'objects':[]}
            obsdict[snid]['dates'].append(r['date'])
            obsdict[snid]['expnums'].append(r['expnum'])
            obsdict[snid]['filts'].append(r['filt'])
            obsdict[snid]['teffs'].append(r['teff'])
            obsdict[snid]['datestrfs'].append(r['datestrf'])
            obsdict[snid]['objects'].append(r['object'])
        elif str(r['object']) in ysedict2.keys():
            snid = ysedict2[r['object']]
            ra = r['ra']
            dec = r['dec']
            if not snid in obsdict.keys():
                obsdict[snid] = {'dates':[],'expnums':[],'filts':[],'ra':ra,'dec':dec,'datestrfs':[],'teffs':[], 'objects':[]}
            obsdict[snid]['dates'].append(r['date'])
            obsdict[snid]['expnums'].append(r['expnum'])
            obsdict[snid]['filts'].append(r['filt'])
            obsdict[snid]['teffs'].append(r['teff'])
            obsdict[snid]['datestrfs'].append(r['datestrf'])
            obsdict[snid]['objects'].append(r['object'])

    if verbose: print('-'*25)
    cnt = 0
    keys = obsdict.keys()
    values = obsdict.values()
    datestrfs = [min(v['datestrfs']) for v in values]
    if verbose: print(np.array(list(keys)))
    if verbose: print(datestrfs)

    returndict = {}
    ss = np.argsort(datestrfs)
    for k in np.array(list(keys))[ss]:
        #for k,v in zip(keys[ss],values[ss]):
        v = obsdict[k]
        if k in ignore: continue

        if c.ccdmap.get(k, None) is None:
            if not m['SNID'].str.contains(k).any():
                continue
            else:
                ind = m['SNID'].str.contains(k)
                fixccd = m['candCCD'].values[ind][0]
                c.ccdmap[k] = fixccd

        if k in rysedict.keys():
            if verbose: print('SNID',k,'YSE_Field',rysedict[k],'CCD',c.ccdmap[k],'RA',v['ra'],'DEC',v['dec'])
            outtxt.write(' '.join(['SNID',str(k),'YSE_Field',rysedict[k],'CCD',str(c.ccdmap[k]),'\n']))
            outtxtexp.write(' '.join(['SNID',str(k),'YSE_Field',rysedict[k],'CCD',str(c.ccdmap[k]),'\n']))
            returndict[k] = ' '.join(['SNID',str(k),'YSE_Field',rysedict[k],'RA',str(v['ra']),'DEC',str(v['dec']),'\n\n'])
        else:
            if verbose: print('SNID',k,'CCD',c.ccdmap[k],'RA',v['ra'],'DEC',v['dec'])
            outtxt.write(' '.join(['SNID',str(k),'CCD',str(c.ccdmap[k]),'\n']))
            outtxtexp.write(' '.join(['SNID',str(k),'CCD',str(c.ccdmap[k]),'\n']))
            returndict[k] = ' '.join(['SNID',str(k),'RA',str(v['ra']),'DEC',str(v['dec']),'\n\n'])
        #cnt += 1
        if len(np.unique(v['dates'])) > 1:
            cnt += 1
        for date in np.sort(np.unique(v['dates'])):
            ww = np.array(v['dates']) == date
            filts = np.array(v['filts'])[ww]
            exps = np.array(v['expnums'])[ww]
            if verbose: print(date+':',filts,'Avg Teff %.2f'%np.nanmean(np.array(v['teffs'])[ww]),)
            outstr = date+':'+str(filts)+'\n'
            outstre = date+':'+str(filts)+str(exps)+'\n'
            outtxt.write(outstr)
            outtxtexp.write(outstre)
            returndict[k] += date+': ' + str(filts) + ' Avg Teff %.2f'%np.nanmean(np.array(v['teffs'])[ww])+'\n'

        outtxt.write('\n')
        outtxtexp.write('\n')
        if verbose: print()
        #print('Dates',v['dates'])
        if verbose: print('-'*25)
        outtxt.write('-'*25)
        outtxt.write('\n')
        outtxtexp.write('-'*25)
        outtxtexp.write('\n')
    outtxt.close()
    outtxtexp.close()
    if verbose: print('Total SNe %d'%cnt)

    if verbose: print('-'*1000)
    allexps = []
    for k,v in obsdict.items():
        if k in ignore: continue
        #print(v['expnums'])
        #for e,o in zip(v['expnums'],v['objects']):   
        allexps.extend(v['expnums'])
    allexps = np.unique(allexps)

    fout = open('debass/combined_expnums.txt','w')
    for exp in allexps:
        fout.write(str(exp)+'\n')
    fout.close()

    fout = open('debass/debass_expnums.txt','w')
    for exp in np.unique(debassonlyexpnums):
        fout.write(str(exp)+'\n')
    fout.close()

    return returndict
    
if __name__ == '__main__':
    run()
