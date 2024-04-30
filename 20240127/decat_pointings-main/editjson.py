import json
import os
from copy import copy


def edit(file,priority,filters,exptimes,outfile):
    with open(file) as f:
        data = json.load(f)
        pointing = copy(data[0])
        outlist = []
        for filt,e in zip(filters,exptimes):
            pointing['filter'] = filt
            pointing['exptime'] = float(e)
            pointing['comment'] = pointing['comment'].replace('_P1','_P%s'%priority).replace('_P2','_P%s'%priority).replace('_P3','_P%s'%priority)
            pointing['object'] = pointing['object'].replace('_P1','_P%s'%priority).replace('_P2','_P%s'%priority).replace('_P3','_P%s'%priority)
            outlist.append(copy(pointing))
    with open(outfile, 'w') as json_file:
        json.dump(outlist, json_file, indent=4)

def getfiltersexptimes(file):
    dict = {}
    with open(file) as f:
        data = json.load(f)
        for d in data:
            dict[d['filter']] = d['exptime']
    return dict
        
if __name__ == '__main__':
    edit('jsons/2020B-0053_DEBASS_Brout/TEMPLATE/2021inj.json',1,['g','r','i','z','Y'],['15','15','15','15','15'],'jsons/2020B-0053_DEBASS_Brout/test/2021inj_P1.json')
