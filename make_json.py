import numpy as np

#make_json.individual('json_out',
#                     ['525.F.a','525.F.b','525.E.a','525.A.a','525.A.b'], #'525.B.a','525.B.b'],
#                     [184.8,187.8,186.6,186.3,188.5], #,189.3,186.9],
#                     [6.1,7.2,12.3,14.5,14.5], #9.8,9.9],
#                     ['g50,r50','g50,r50','g50,r50','g50,r50','g50,r50'],
#                     ['2021A-0275','2021A-0275','2021A-0275','2021A-0275','2021A-0275'],
#                     ['','','','',''],
#                     ['object','object','object','object','object'],
#                     ['','','','',''])

sniddict = {
    '415.A.a':('2021cyn',2,36),
    '1511.B.a':('2021cym',4,16),
    '415.A.b':('2021ddh',1,36),
    '415.B.a':('2021dov',1,29),
    '418.F.a':('2021dev',1,11),
    '523.D.a':('2021dha',1,9),
    '523.C.a':('2021dib',2,15),
    '575.B.a':('2021ckc',2,6),
    '525.D.a':('2021chy',2,28),
    '474.E.a':('2021dgu',2,23),
    '577.B.a':('2021dhg',2,22),
    '577.A.a':('2021dch',1,22),
    '674.E.a':('2021dep',2,42),
    '376.2020esm.a':('2021cjb',2,48),
    '477.F.a':('2021dam',2,36),
    '428.B.a':('2021dlb',2,41),
    '632.C.a':('2021dcu',1,42),
    '477.E.a':('2021dur',1,16),
    '423.E.a.2021fxy':('2021fxy',2,34),
    '583.F.a':('2021gez',1,6),
    '531.E.a':('2021haz',1,36)}

_tnsapi = 'https://www.wis-tns.org/api/get'
_tnsapikey = 'ecd2dec8cee4ed72a39fe8467ddd405fec4eef14'

import requests
import json
from collections import OrderedDict

def get(url,json_list,api_key):
    try:
        # url for get obj
        get_url=url+'/object'
        # change json_list to json format
        json_file=OrderedDict(json_list)
        # construct the list of (key,value) pairs
        get_data=[('api_key',(None, api_key)),
                  ('data',(None,json.dumps(json_file)))]
        # get obj using request module
        response=requests.post(get_url, files=get_data)
        return response
    except Exception as e:
        return [None,'Error message : \n'+str(e)]

def format_to_json(source):
    # change data to json format and return
    parsed=json.loads(source,object_pairs_hook=OrderedDict)
    return parsed

def mklist():
    field,fra,fdec = np.loadtxt('fields.txt',unpack=True,dtype=str)
    
    for k in sniddict.keys():
        TNSGetSingle = [("objname",sniddict[k][0]),
                        ("photometry","0"),
                        ("spectra","0")]
        response=get(_tnsapi, TNSGetSingle, _tnsapikey)
        json_data = format_to_json(response.text)
        ra,dec = json_data['data']['reply']['radeg'],json_data['data']['reply']['decdeg']

        print(f"{k} {float(fra[k == field][0]):.6f} {float(fdec[k == field][0]):.6f} {sniddict[k][0]} {ra:.6f} {dec:.6f} {sniddict[k][2]} {sniddict[k][1]}")
        
def all():
    snid,ra,dec = np.loadtxt('fields.txt',unpack=True,dtype=str)
    ra,dec = ra.astype(float),dec.astype(float)
    #import pdb; pdb.set_trace()
    for s,r,d in zip(snid,ra,dec):
        individual('jsons/2021A-0275_YSE_Rest',[s],[r],[d],['i50,z50'],['2021A-0275'],[s],['object'],[''])
        
def individual(json_outpath,json_prefixs,pointRAs,pointDECs,obss,propids,objects,exptypes,programs):
    for json_prefix,pointRA,pointDEC,obs,propid,tobject,exptype,program in zip(json_prefixs,pointRAs,pointDECs,obss,propids,objects,exptypes,programs):
        json_out = []
        json_out.append('[')
        cntr = 0
        lenobs = len(obs.split(','))
        for ob in obs.split(','):
            cntr += 1
            filt = ob[0]
            exptime = ob[1:]
            json_out.append('\t{') 
            json_out.append('\t\t"count": 1,')
            json_out.append('\t\t"filter": "%s",'%filt)
            json_out.append('\t\t"exptime": %.1f,'%float(exptime.strip()))
            json_out.append('\t\t"RA": %f,' % round(float(pointRA), 5))
            json_out.append('\t\t"dec": %f,' % round(float(pointDEC), 5))
            json_out.append('\t\t"object": "%s",'%tobject)
            json_out.append('\t\t"program": "%s",'%program)
            json_out.append('\t\t"expType": "%s",'%exptype)
            json_out.append('\t\t"note": "None",')
            if tobject in sniddict.keys():
                json_out.append('\t\t"comment": "%s %s",'%(json_prefix,sniddict[tobject][0]))
            else:
                json_out.append('\t\t"comment": "%s",'%json_prefix)
            json_out.append('\t\t"wait": "False",')
            json_out.append('\t\t"propid": "%s"'%propid)
            if cntr == lenobs:
                json_out.append('\t}')
            else:            
                json_out.append('\t},')
        json_out.append(']')
        
        np.savetxt(json_outpath+'/'+json_prefix+'.json', json_out, fmt='%s')
        
if __name__ == "__main__":
    all()
    #mklist()
