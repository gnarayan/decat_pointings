import numpy as np
import json as pyjson

readout = 29.

#from here http://www.ctio.noao.edu/noao/node/5826
def slewtime_from_two_coords(ra1,dec1,ra2,dec2):
    delra = (ra1-ra2)*np.cos(np.mean([dec1,dec2]))
    deg_separation = np.sqrt((delra)**2+(dec1-dec2)**2)
    return max([0.,20.-readout+220./100.*deg_separation])#subtracted off readout and bounded by zero

def time_from_list_of_ras_decs_exptimes(ras,decs,exptimes):
    ras,decs,exptimes = np.array(ras),np.array(decs),np.array(exptimes) 
    totaltime = 0.
    for i,ra,dec,exptime in zip(range(len(ras)),ras,decs,exptimes):
        if i==0:
            slewtime = 0
        else:
            slewtime = slewtime_from_two_coords(ra,dec,ras[i-1],decs[i-1])
        totaltime += exptime+slewtime+readout
    return totaltime
        
def get_ras_decs_exptimes_from_json(json):
    ras,decs,exptimes = [],[],[]
    with open(json) as json_file: 
        pointings = pyjson.load(json_file) 
        for pointing in pointings:
            ras.append(float(pointing['RA']))
            decname = np.array(list(pointing.keys()))[np.array([k.lower() == 'dec' for k in pointing.keys()])][0]
            decs.append(float(pointing[decname]))
            expname = np.array(list(pointing.keys()))[np.array([k.lower() == 'exptime' for k in pointing.keys()])][0]
            exptimes.append(float(pointing[expname]))
    return ras,decs,exptimes

def time_for_single_json(json):
    #json: filename
    ras,decs,exptimes = get_ras_decs_exptimes_from_json(json)
    totaltime = time_from_list_of_ras_decs_exptimes(ras,decs,exptimes)
    return totaltime

def total_time_from_jsons(jsons):
    #jsons: list of filenames
    ras,decs,exptimes = [],[],[]
    for json in jsons:
        print(json)
        tras,tdecs,texptimes = get_ras_decs_exptimes_from_json(json)
        ras.extend(tras)
        decs.extend(tdecs)
        exptimes.extend(texptimes)
        print(np.sum(texptimes)/60)
    totaltime =	time_from_list_of_ras_decs_exptimes(ras,decs,exptimes)
    return totaltime
        

#USAGE:
#t = time_for_single_json('/Users/djbrout/Dropbox/decat_pointings/jsons/20210223/2021dep_211.053_27.626_PRIORITY2.json')
#print(t)
#t = total_time_from_jsons(['/Users/djbrout/Dropbox/decat_pointings/jsons/20210223/2021dep_211.053_27.626_PRIORITY2.json',
#                           '/Users/djbrout/Dropbox/decat_pointings/jsons/20210223/2021dnl_221.762_0.476_PRIORITY2.json',
#                           '/Users/djbrout/Dropbox/decat_pointings/jsons/20210223/2021dov_134.078_0.442_PRIORITY1.json'])
#print(t)
