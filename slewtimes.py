import numpy as np
import json as pyjson
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u

readout = 29.

#from here http://www.ctio.noao.edu/noao/node/5826
def slewtime_from_two_coords(ra1,dec1,ra2,dec2):
    delra = (ra1-ra2)*np.cos(np.mean([dec1,dec2]))
    deg_separation = np.sqrt((delra)**2+(dec1-dec2)**2)
    return max([0.,20.-readout+220./100.*deg_separation])#subtracted off readout and bounded by zero

def time_from_list_of_ras_decs_exptimes(ids,ras,decs,exptimes,sort=True):
    ids,ras,decs,exptimes = np.array(ids),np.array(ras),np.array(decs),np.array(exptimes)
    df = pd.DataFrame.from_dict({'ids':ids,'rra':np.round(ras,-1),'ra':ras,'dec':decs,'exptime':exptimes})
    if sort:
        df.sort_values(by=['rra', 'dec'],inplace=True)
        df = df.reset_index(drop=True)
    totaltime = 0.
    totaltime_p1 = 0
    lastid=''
    iobserve = open('iobserve.txt','w')
    pout = []
    for i,row in df.iterrows():
        #print(row['ids'])
        if i==0:
            slewtime = 0
        else:
            #print(row['ra'],row['dec'])
            if row['ids'] != lastid:
                iobserve.write('%s %d.0 %d.0\n'%(row['ids'],row['ra'],row['dec']))
                pout.append('%s [%d %d]'%(row['ids'],row['ra'],row['dec']))
                lastid=row['ids']
            slewtime = slewtime_from_two_coords(row['ra'],row['dec'],df['ra'][i-1],df['dec'][i-1])
        totaltime += row['exptime']+slewtime+readout
        if '_P1' in row['ids']:
            totaltime_p1 += row['exptime']+slewtime+readout
    iobserve.close()
    for p in pout:
        if '_P1' in p:
            print(p)
    print()
    for	p in pout:
        if '_P2' in p:
            print(p+' OPTIONAL')
    print()
    for p in pout:
        if '_P3' in p:
            print(p)
    print()
    for p in pout:
        if '_TCTM' in p:
            print(p)

    return totaltime,totaltime_p1
        
def get_ras_decs_exptimes_from_json(json):
    ids,fns, ras,decs,exptimes = [],[],[],[],[]
    with open(json) as json_file: 
        pointings = pyjson.load(json_file)
        lastra = np.nan
        lastdec = np.nan
        for pointing in pointings:
            try:
                raname = np.array(list(pointing.keys()))[np.array([k.lower() == 'ra' for k in pointing.keys()])][0]
                ra = pointing[raname]
                decname = np.array(list(pointing.keys()))[np.array([k.lower() == 'dec' for k in pointing.keys()])][0]
                dec = pointing[decname]
            except:
                print('ra/dec not found, using previous exposure ra/dec')
                ra = lastra
                dec = lastdec
            if ":" in str(ra):
                c = SkyCoord(str(ra)+' '+str(dec), unit=(u.hourangle, u.deg))
                ra = c.ra.degree
                dec = c.dec.degree
            
            ids.append(pointing['object'])
            ras.append(float(ra))
            decs.append(float(dec))
            fns.append(json.split('/')[-1].split('.json')[0])
            lastra = ra
            lastdec = dec
            expname = np.array(list(pointing.keys()))[np.array([k.lower() == 'exptime' for k in pointing.keys()])][0]
            exptimes.append(float(pointing[expname]))
    return fns,ras,decs,exptimes

def time_for_single_json(json):
    #json: filename
    ids, ras,decs,exptimes = get_ras_decs_exptimes_from_json(json)
    totaltime = time_from_list_of_ras_decs_exptimes(ids,ras,decs,exptimes)
    return totaltime

def total_time_from_jsons(jsons,sort=True):
    #jsons: list of filenames
    ids, ras,decs,exptimes = [],[],[],[]
    for json in jsons:
        #print(json)
        tids, tras,tdecs,texptimes = get_ras_decs_exptimes_from_json(json)
        ids.extend(tids)
        ras.extend(tras)
        decs.extend(tdecs)
        exptimes.extend(texptimes)
        #print(np.sum(texptimes)/60)
        #print(tids[0],tras[0],tdecs[0])
    totaltime =	time_from_list_of_ras_decs_exptimes(ids,ras,decs,exptimes,sort=sort)
    return totaltime
        

#USAGE:
#t = time_for_single_json('/Users/djbrout/Dropbox/decat_pointings/jsons/20210223/2021dep_211.053_27.626_PRIORITY2.json')
#print(t)
#t = total_time_from_jsons(['/Users/djbrout/Dropbox/decat_pointings/jsons/20210223/2021dep_211.053_27.626_PRIORITY2.json',
#                           '/Users/djbrout/Dropbox/decat_pointings/jsons/20210223/2021dnl_221.762_0.476_PRIORITY2.json',
#                           '/Users/djbrout/Dropbox/decat_pointings/jsons/20210223/2021dov_134.078_0.442_PRIORITY1.json'])
#print(t)
