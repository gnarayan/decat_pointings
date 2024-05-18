#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 19 23:34:52 2021

@author: arest
"""

import sys,os,argparse,re
import glob
import numpy as np
import pandas as pd

def unique(A):
    unique = []
    for a in A:
        if a not in unique:
            unique.append(a)
    return unique

def makepath(path,raiseError=True):
    if path == '':
        return(0)
    if not os.path.isdir(path):
        os.makedirs(path)
        if not os.path.isdir(path):
            if raiseError:
                raise RuntimeError('ERROR: Cannot create directory %s' % path)
            else:
                return(1)
    return(0)

def makepath4file(filename,raiseError=True):
    path = os.path.dirname(filename)
    if not os.path.isdir(path):
        return(makepath(path,raiseError=raiseError))
    else:
        return(0)




if __name__ == "__main__":
    files = glob.glob("./???.?.?.json")
    for i in range(len(files)): files[i]=os.path.basename(files[i])
    files.sort()
    t = pd.DataFrame(data=files,columns=['filename'])
    t['field']=t['filename'].str.extract(r'(^\w+)')
    fields = unique(t['field'])
    fields.sort()

    print(t)
    print(fields)
    
    for field in fields:
        ixs = t.index.values
        (keep,) = np.where(t.loc[ixs,'field'].eq(field))
        ixs_field = ixs[keep]
        #print(t.loc[ixs_field])
        
        newjson = ['[\n']
        for ix in ixs_field:
            lines = open(t.loc[ix,'filename'],'r').readlines()
            del(lines[-1])
            del(lines[0])
            if ix != ixs_field[-1]:
                lines[-1]=lines[-1].rstrip()+',\n'
            #print(lines)
            newjson.extend(lines)
        newjson[-1] = re.sub(',','',newjson[-1])
        newjson.append(']\n')
        outfilename = f'./combined/{field}.json'
        makepath4file(outfilename)
        print(f'writing {outfilename}')
        open(outfilename,'w').writelines(newjson)
      
        
        