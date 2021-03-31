#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 17:28:01 2021

@author: arest
"""
import sys, os, re, copy, shutil,io
import argparse
from pdastro import pdastroclass,AnotB,AandB
import pandas as pd
import numpy as np
from astropy import time

programlist = [
'2019A-0065_Shen1',
'2019B-0304_Martini',
'2020A-0906_eFEDS',
'2021A-0037_Shen2',
'2021A-0275_YSE',
'2020B-0053_DEBASS',
'2021A-0113_DDF',
'2021A-0148_DESI',
'2020A-0415_EtaCar'
]

program2fieldpattern = {
'2019A-0065_Shen1':      ['^SN\-C3','^S\-CVZ','SN\-X\d'],
'2019B-0304_Martini':   ['E1','E3'],
'2020A-0906_eFEDS':     ['^eFEDS'],
'2021A-0037_Shen2':      ['^CO\d$'],
'2021A-0275_YSE':       ['^\d\d\d\.\w+\.[abcde]'],
'2020B-0053_DEBASS':    ['^2021\w+'],
'2021A-0113_DDF':       ['^COSMOS','^DECaPS.*'],
'2021A-0148_DESI':      ['^TILEID\:\s+\d+'], 
'2020A-0415_EtaCar':    ['^ec\d\d\d\d'],
'2019A-0305_Drlica_TRADE':['^DELVE'],
'2021A-0244_Miller_TRADE':['^n2997'],
'2021A-0010_Rector_TRADE':['^Cha'],
'2021A-0149_Zenteno_TRADE':['^BLA'],
'STANDARDS':             ['^E','^SDSS','^LTT','C26202'],
'TECHSETUP':               ['^pointing','^MaxVis']
    }


class calcTimeclass(pdastroclass):
    def __init__(self):
        pdastroclass.__init__(self)
        
        self.qcinv = pdastroclass()
        self.minimal_overhead = 28 # in seconds
        
        self.verbose=0
        self.debug=0
        
        self.warnings = []
        
        self.t = pd.DataFrame(columns=['blockID','assigned_program','program','UTfirst','UTlast','dt_block_h','dt_prevgap_sec','dt_nextgap_sec','dt_gaps_sec','dt_block_full_h'])
        self.summary = pdastroclass()
        self.summary.t = pd.DataFrame(columns=['assigned_program','t_total'])
        
        self.programcol_formatter='{:<24}'.format

    def addwarning(self,warningstring):
        print(warningstring)
        self.warnings.append(warningstring)
        return(0)


    def add_arguments(self, parser=None, usage=None, conflict_handler='resolve'):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)
            
        parser.add_argument('qcinvfile')

        parser.add_argument('-s','--save', nargs='*', help="Save the tables. if no argument is specified, then the input name is used as basename")
        parser.add_argument('-r','--reassign', nargs=2, action="append", help="reassign one program block to another program")
            
        parser.add_argument('-v', '--verbose', action='count', default=0)
        parser.add_argument('-d', '--debug', action='count', default=0)

        return(parser)
    
    def create_fieldpattern2program(self):
        self.fieldpattern2program ={}
        self.fieldpatterns=[]
        for program in program2fieldpattern:
            for fieldpattern in program2fieldpattern[program]:
                self.fieldpatterns.append(fieldpattern)
                self.fieldpattern2program[fieldpattern]={}
                self.fieldpattern2program[fieldpattern]['program']=program
                self.fieldpattern2program[fieldpattern]['compiled']=re.compile(fieldpattern)

        return(0)            
    
    def readqcinv(self,filename):
        # load the qcinv file and fix it: remove \n and remove extra headers
        if self.verbose: print('loading ',filename)
        self.filename = filename
        lines = open(filename,'r').readlines()
        for i in range(len(lines)-1,-1,-1):
            lines[i] = re.sub('\\n$','',lines[i])
            if re.search('^\s*\#',lines[i]):
                if i==0:
                    lines[i] = re.sub('^\#',' ',lines[i])
                else:
                    del(lines[i])

        #insert dummy line: otherwise the last column width is set to the width of the last column in the first row
        s = re.sub('\w+$','dummydummydummydummydummy',lines[1])
        lines.insert(1,s)
        
        # 'time' is too long for the column width and butts into the secz column,
        # which confuses pandas reading, and it combines the time and secz column.
        # hack: rename time to tim
        lines[0]=re.sub('time','tim ',lines[0])
        
        # Remove the extra "MJD = ..." line at the end
        if re.search('^MJD',lines[-1]): 
            lines.pop(-1)

        # parse the qcinv file
        self.qcinv.t = pd.read_fwf(io.StringIO('\n'.join(lines)))
        #self.qcinv.write(indices=range(1,10))
        
        # rename tim to time again...
        self.qcinv.t.rename(columns={'tim': 'time'},inplace=True)
        
        #remove dummy line
        self.qcinv.t = self.qcinv.t[1:]
        #self.qcinv.t.drop(index=[0],inplace=True)

        if self.verbose:
            print('file loaded, first 10 lines:')
            self.qcinv.write(indices=range(1,11))
        return(0)
    
    def assignPrograms(self):
        self.create_fieldpattern2program()
        self.qcinv.t['program']=None
        self.qcinv.t['blockID']=0
        
        ixs = self.qcinv.getindices()
        
        blockID = 1
        
        special_programs = {}
        for p in ['UNKNOWN','TECHSETUP','STANDARDS']:
            special_programs[p]={}
            special_programs[p]['counter']=0
            special_programs[p]['pattern']=re.compile('^%s' % p)
        
        for i in range(len(ixs)):
            ix = ixs[i]

            foundflag=False
            for fieldpattern in self.fieldpattern2program:
                m = self.fieldpattern2program[fieldpattern]['compiled']
                if m.search(self.qcinv.t.loc[ix,'Object']):
                    program = self.fieldpattern2program[fieldpattern]['program']
                    if self.verbose>2: print('FOUND! pattern %s matches %s, program %s' % (fieldpattern,self.qcinv.t.loc[ix,'Object'],program))
                    if program=='2020B-0053_DEBASS':
                        if self.qcinv.t.loc[ix,'time']>15:
                            program='2021A-0275_YSE'

                    if program in special_programs:
                        m = special_programs[program]['pattern']
                        if i==0 or (not m.search(self.qcinv.t.loc[ixs[i-1],'program'])):
                            special_programs[program]['counter']+=1
                        program += '%d' %  special_programs[program]['counter']
 
                    self.qcinv.t.loc[ix,'program']=program
                    foundflag=True
                    break
            if not foundflag:
                self.addwarning('WARNING: Could not find the program for %s in line %d' % (self.qcinv.t.loc[ix,'Object'],ix))
                program='UNKNOWN' 
                m = special_programs[program]['pattern']
                if i==0 or (not m.search(self.qcinv.t.loc[ixs[i-1],'program'])):
                    special_programs[program]['counter']+=1
                program += '%d' %  special_programs[program]['counter']
                self.qcinv.t.loc[ix,'program']=program
            
            # if not the first entry AND if different than previous row's program: inc blockID
            if ix!=ixs[0] and self.qcinv.t.loc[ixs[i-1],'program']!=self.qcinv.t.loc[ix,'program']:
                blockID+=1
            self.qcinv.t.loc[ix,'blockID']=blockID

            
    def calcTimes(self):
        blockIDs = self.qcinv.t['blockID'].unique()
        m = re.compile('^2')
        
        # set format of dt_h*.
        self.qcinv.t['ut_decimal']=np.nan
        self.qcinv.default_formatters['ut_decimal']='{:.4f}'.format
        self.t['dt_block_h']=self.t['dt_prevgap_sec']=self.t['dt_nextgap_sec']=self.t['dt_gaps_sec']=self.t['dt_block_full_h']=np.nan
        self.default_formatters['dt_block_h']='{:.4f}'.format
        self.default_formatters['dt_prevgap_sec']='{:.0f}'.format
        self.default_formatters['dt_nextgap_sec']='{:.0f}'.format
        self.default_formatters['dt_gaps_sec']='{:.0f}'.format
        self.default_formatters['dt_block_full_h']='{:.4f}'.format
        self.default_formatters['assigned_program']=self.programcol_formatter
        
        ix_all = self.qcinv.getindices()
        
        re_techsetup_standards = re.compile('^TECHSETUP|^STANDARDS')
        
        # I just choose a random date. The date itself is not important, it's 
        # just that all UT times 2?:?? have the date before this random date,
        # so that the time difference is correct
        t0 = time.Time('2020-01-02T00:00:00',scale='utc',format='isot')
        
        # get time difference of ut to t0, save it in dt_h in hours
        for ix in ix_all:
            # choose date depending on whether the time is before or after UT midnight
            if m.search(self.qcinv.t.loc[ix,'ut']):
                s = '2020-01-01T'+self.qcinv.t.loc[ix,'ut']+':00'
            else:
                s = '2020-01-02T'+self.qcinv.t.loc[ix,'ut']+':00'
                
            t1 = time.Time(s,scale='utc',format='isot')
            dt = (t1-t0)
            self.qcinv.t.loc[ix,'ut_decimal'] = dt.to_value('hr')
        
        # get info for each block
        for i in range(len(blockIDs)):
            ixs = self.qcinv.ix_inrange('blockID', blockIDs[i],blockIDs[i])
            if len(ixs)==0:
                self.newrow({'blockID':blockIDs[i]})
                self.addwarning('WARNING: could not find any entries for blockID %d' % blockIDs[i])
                continue
            
            # first get the difference between last and first UT
            dt_block_h = self.qcinv.t.loc[ixs[-1],'ut_decimal']-self.qcinv.t.loc[ixs[0],'ut_decimal'] 
            
            # exposure time and nominal overhead for last im
            lastim_dt = (self.qcinv.t.loc[ixs[-1],'time']+self.minimal_overhead)/3600.0
            
            # now add in exposure time of last ecposure and nominal overhead
            dt_block_h += lastim_dt

            # get info for gap
            if i < len(blockIDs)-1:
                ixs_next = self.qcinv.ix_inrange('blockID', blockIDs[i]+1,blockIDs[i]+1)
                dt_nextgap_sec = self.qcinv.t.loc[ixs_next[0],'ut_decimal']-(self.qcinv.t.loc[ixs[-1],'ut_decimal']+lastim_dt)
            else:
                dt_nextgap_sec=0.0
                
            self.newrow({'blockID':blockIDs[i],
                         'assigned_program':self.qcinv.t.loc[ixs[0],'program'],
                         'UTfirst':self.qcinv.t.loc[ixs[0],'ut'],
                         'UTlast':self.qcinv.t.loc[ixs[-1],'ut'],
                         'dt_block_h':dt_block_h,
                         'dt_prevgap_sec':1.0,                    
                         'dt_nextgap_sec':dt_nextgap_sec*3600.0                       
                         })
        
        # assign previous gap times
        # Also, take into account the edge cases of first and last block
        ixs_blocks = self.getindices()
        for i in range(len(ixs_blocks)):
            
            if i==0: 
                self.t.loc[ixs_blocks[i],'dt_prevgap_sec']=0.0
                continue
            
            if i==len(ixs_blocks)-1:
                self.t.loc[ixs_blocks[i],'dt_nextgap_sec']=0.0
            
            self.t.loc[ixs_blocks[i],'dt_prevgap_sec']=self.t.loc[ixs_blocks[i-1],'dt_nextgap_sec']
        
        # If TECHSETUP and STANDARDS is block: eat all the time, so set nextgap from previous block and prevgap from next block to 0.0 !
        for i in range(len(ixs_blocks)-1):
            if re_techsetup_standards.search(self.t.loc[ixs_blocks[i],'assigned_program']):
                if i>0: self.t.loc[ixs_blocks[i-1],'dt_nextgap_sec']=0.0
                if i<len(ixs_blocks)-1: self.t.loc[ixs_blocks[i+1],'dt_prevgap_sec']=0.0
        
        # Add up block time and gap time:
        # TECHSETUP, STANDARDS: eat all gap time
        # all other programs: eat half of each gap time (previous and next)
        for i in range(len(ixs_blocks)):
            if re_techsetup_standards.search(self.t.loc[ixs_blocks[i],'assigned_program']):
                self.t.loc[ixs_blocks[i],'dt_gaps_sec'] = 1.0*(self.t.loc[ixs_blocks[i],'dt_nextgap_sec']+self.t.loc[ixs_blocks[i],'dt_prevgap_sec'])
            else:
                self.t.loc[ixs_blocks[i],'dt_gaps_sec'] = 0.5*(self.t.loc[ixs_blocks[i],'dt_nextgap_sec']+self.t.loc[ixs_blocks[i],'dt_prevgap_sec'])
            
            self.t.loc[ixs_blocks[i],'dt_block_full_h'] = self.t.loc[ixs_blocks[i],'dt_block_h']+self.t.loc[ixs_blocks[i],'dt_gaps_sec']/3600.0
                
        self.t['program']='-'
                
            
    def reassign_programs(self,reassign):
        if reassign is not None:
            
            for (inprogram,outprogram) in reassign:
                ixs = self.ix_equal('assigned_program',inprogram)
                if len(ixs)==0:
                    self.addwarning('WARNING! could not find any entries for program %s, so that it can be reassigned to %s' % (inprogram,outprogram))
                    continue
                self.t.loc[ixs,'program']=inprogram
                self.t.loc[ixs,'assigned_program']=outprogram
        
    def mkSummary(self):
        current_programs = self.t['assigned_program'].unique()
        current_mainprograms = AandB(current_programs,programlist)
        extra_programs = AnotB(current_programs,current_mainprograms)
        extra_programs.sort()
        current_mainprograms.sort()
        programs = list(current_mainprograms)
        programs.extend(extra_programs)
        for program in programs:
            ixs = self.ix_equal('assigned_program',program)
            t_total = self.t.loc[ixs,'dt_block_full_h'].sum()
            self.summary.newrow({'assigned_program':program,
                                 't_total':t_total})
        self.summary.default_formatters['assigned_program']=self.programcol_formatter
        self.summary.write()
        
    def savetables(self,basename=None):
        if basename is not None:
            if len(basename)==0:
                basename=self.filename
            elif len(basename)==1:
                basename=basename[0]
            else:
                raise RuntimeError('more than one argument for --save is not allowed!')
            print('Saving tables with basename',basename)
            self.qcinv.write(basename+'.times.txt',verbose=2)
            self.write(basename+'.blocks.txt',verbose=2)
            self.summary.write(basename+'.summary.txt',verbose=2)
            
            
        return(0)

if __name__ == "__main__":
    calcTime = calcTimeclass()
    usagestring='USAGE: calcTime.py qcinv_filename'
    parser=calcTime.add_arguments(usage=usagestring)
    args = parser.parse_args()
    
    calcTime.verbose=args.verbose
    calcTime.debug=args.debug

    calcTime.readqcinv(args.qcinvfile)    
    calcTime.assignPrograms()
    calcTime.calcTimes()
    calcTime.reassign_programs(args.reassign)
    if args.verbose:
        calcTime.qcinv.write()
        calcTime.write()

    calcTime.mkSummary()
    
    calcTime.savetables(basename=args.save)
    
    if len(calcTime.warnings)>0:
        print('THERE WERE WARNINGS!!!')
        for s in calcTime.warnings: print(s)
