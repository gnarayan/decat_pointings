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
from astroplan import Observer
from astropy import units as u
ctio = Observer.at_site("CTIO")

from mk_semester_summary import semester_summary_class
from semesterinfo import semesterinfoclass,default_semester


class calcTimeclass(pdastroclass):
    def __init__(self):
        pdastroclass.__init__(self)
        
        self.qcinv = pdastroclass()
        self.minimal_overhead = 28 # in seconds
        
        self.verbose=0
        self.debug=0
        
        self.warnings = []
        
        self.t = pd.DataFrame(columns=['blockID','assigned_program','program','UTfirst','UTlast','dt_block_h','dt_prevgap_sec','dt_nextgap_sec','dt_gaps_sec','dt_block_full_h','twi','dt_charged_h'])
        self.nightsummary = pdastroclass(columns=['assigned_program','t_dark','t_twi1','t_twi2','t_twi3','t_down'])
        #self.nightsummary.t = pd.DataFrame(columns=['assigned_program','t_dark','t_twi','t_downtime'])
        self.nightsummary.default_formatters = {'t_dark':'{:.4f}'.format,
                                                't_twi1':'{:.4f}'.format,
                                               't_twi2':'{:.4f}'.format,
                                               't_twi3':'{:.4f}'.format,
                                               't_down':'{:.4f}'.format
                                               }
        
        self.semestersummary = semester_summary_class()
        
        

        self.programcol_formatter='{:<24}'.format
        
        self.horizons = [18,15,12]
        #self.twi_charge_fraction = (1.0,2/3,1/3,0.0)
        self.twi_charge_fraction = (1.0,1.0,1.0,1.0)
         
        self.downtime_reassign = []
        
        self.semesterinfo = semesterinfoclass()

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
        parser.add_argument('--downtime', nargs='+', action="append", help="specify downtime: Last_exp_ID program1 time1 program2 time2 ...")
        parser.add_argument('--semester', default=default_semester, help='Specify the semester (default=%(default)s)')
        parser.add_argument('--semester_summaryfile', default=None, help='Specify the filename for the semester summary. If not specified, it will be <semester>/hours_summary.txt (default=%(default)s)')
        parser.add_argument('-a','--add2semestersummary', default=False,action='store_true',help='Save the hours into the semester summary')
        parser.add_argument('-v', '--verbose', action='count', default=0)
        parser.add_argument('-d', '--debug', action='count', default=0)

        return(parser)
    
    def create_fieldpattern2program(self):
        self.fieldpattern2program ={}
        self.fieldpatterns=[]
        for program in self.semesterinfo.program2fieldpattern:
            for fieldpattern in self.semesterinfo.program2fieldpattern[program]:
                self.fieldpatterns.append(fieldpattern)
                self.fieldpattern2program[fieldpattern]={}
                self.fieldpattern2program[fieldpattern]['program']=program
                self.fieldpattern2program[fieldpattern]['compiled']=re.compile(fieldpattern)

        return(0)            
    
    def readqcinv(self,filename):
        print('filename argument:',filename)
        if not os.path.isfile(filename) and re.search('^\d+$',filename):
            newfilename = '%s/%s/%s.qcinv' % (self.semesterinfo.semester,filename,filename)
            if not os.path.isfile(newfilename):
                raise RuntimeError('filename argument %s extended to filename %s, but doesn\'t exist!' % (filename,newfilename))
            filename=newfilename
        
        # load the qcinv file and fix it: remove \n and remove extra headers
        if self.verbose: print('loading ',filename)
        self.qcinv.filename = filename
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
        lines[0]=re.sub('time','xx  ',lines[0])
        
        # remove any empty lines at the end of the lines
        while lines[-1]=='':
            lines.pop(-1)
            
        # Remove the extra "MJD = ..." line at the end
        if re.search('^MJD',lines[-1]): 
            m = re.search('^MJD\s*=\s*(\w+)',lines[-1])
            if m is None:
                raise RuntimeError('Could not get the MJD from the qcinv file! the last line is %s, but should be something like "MJD = 59146 (Oct 23/Oct 24)"' % (lines[-1]))
            MJD = int(m.groups()[0])
            print('MJD: %d' % MJD)
            
            # define tonight for CTIO.
            # Subtract 0.4 from MJD: if the MJD is not *before* the night starts, it tonight
            # returns as start value the given MJD. Subtracting 0.4 makes sure MJD is before the night starts.
            self.tonight = ctio.tonight(time.Time(MJD-0.4, scale='utc',format='mjd'))
            print('Night Start:',self.tonight[0].to_value('isot'))
            print('Night End  :',self.tonight[1].to_value('isot'))

            self.twi={}
            for horizon in self.horizons:
                self.twi[horizon] =  ctio.tonight(time.Time(MJD-0.4, scale='utc',format='mjd'),horizon=-horizon*u.deg)
                if self.verbose: print('%d deg twilight: %s %s' % (-horizon,self.twi[horizon][0].to_value('isot'),self.twi[horizon][1].to_value('isot')))
            # Remove line with MJD
            lines.pop(-1)
        else:
            raise RuntimeError('Could not get the MJD from the qcinv file! the last line should be something like "MJD = 59146 (Oct 23/Oct 24)"')
        
        # parse the qcinv file
        self.qcinv.t = pd.read_fwf(io.StringIO('\n'.join(lines)))
        #self.qcinv.write(indices=range(1,10))
        
        # rename tim to time again...
        self.qcinv.t.rename(columns={'xx': 'time'},inplace=True)
       
        #remove dummy line
        self.qcinv.t = self.qcinv.t[1:]
        #self.qcinv.t.drop(index=[0],inplace=True)

        if self.verbose:
            print('file loaded, first 10 lines:')
            self.qcinv.write(indices=range(1,11))
        return(0)
    
    def downtime2qcinv(self,downtimelist):
        """
        Go through the list of Downtimes, and add rows to qcinv list accordingly.
        The downtimes are automatically enumerated.

        Parameters
        ----------
        downtimelist : TYPE
            1st entry is expID after which downtime will be added.
            The next entries are are pairs of program+downtime. 

        Returns
        -------
        None.

        """
        if args.downtime is None:
            return(0)
        counter_downtime = 1
        downtimerows = pdastroclass(columns=['expid','time','Object','utdate','twi'])

        for downtime in downtimelist:
            print('\n############\n###DOWNTIME: ',downtime)
            # The first entry of the downtime list is the expID after which the downtime should be added
            expID = int(downtime[0])
            
            # if expID = 0, then downtime is starting at -18 deg twilight
            # in this block, we determine utdown and utdownmax: the UT range of the downtime.
            if expID == 0 or expID == 'None':
                utdown = self.twi[18][0]
                print('Downtime starting at evening -18 deg twilight: %s' % utdown.to_value('isot'))
                ix = None
                
                if len(self.qcinv.t.index.values)>0:
                    ix_next = self.qcinv.t.index.values[0]
                else: 
                    ix_next = None
            else:
                print('Downtime starting after expid %d' % expID)
                ixs = self.qcinv.ix_equal('expid',expID)
                if len(ixs)==0:
                    raise RuntimeError('could not find expid %d' % expID)
                elif len(ixs)>1:
                    raise RuntimeError('more than one expid %d' % expID)
                    
                # ix and ix_next are the two indices that sandwich the downtime.
                ix=ixs[0]
                # get the indix of the next row, None if last row
                if ix<len(self.qcinv.t.index.values)-1:
                    ix_next = ix+1
                else:
                    ix_next = None
                
                utdown = time.Time(self.qcinv.t.loc[ix,'utdate']) + (self.qcinv.t.loc[ix,'time']+self.minimal_overhead) * u.s
                #print('VVVV',(self.qcinv.t.loc[ix,'time']+self.minimal_overhead),self.qcinv.t.loc[ix,'utdate'],utdown.to_value('isot'))

            if ix_next is None:
                # morning -18 deg twilight
                utdownmax = self.twi[18][1]
            else:
                utdownmax = time.Time(self.qcinv.t.loc[ix_next,'utdate'])
            
            print('UT range: %s  -  %s' % (utdown.to_value('isot'),utdownmax.to_value('isot')))
            
            # Now loop through all entries of downtime for this downtime block 
            for i in range(1,len(downtime),2):
                print('\n# Downtime entry %d: ' % i,downtime[i],downtime[i+1])
                if utdown.to_value('isot')>utdownmax.to_value('isot'):
                    raise RuntimeError('downtime UT %s is later than the maximum %s!' % (utdown.to_value('isot'),utdownmax.to_value('isot')))
                
                downtimename = 'DOWNTIME%02d' % counter_downtime
                program=downtime[i]
                self.downtime_reassign.append((downtimename,program))
                
                t = downtime[i+1]
                if t.lower() == 'rest':
                    tsec = (utdownmax - utdown).to_value('sec')
                else:
                    if re.search('^[0-9\.]+$',t) is None:
                        (t,unit) = re.search('(^[0-9\.]+)([a-zA-Z]+$)',t).groups()
                    else:
                        unit='sec'
                    print('downtime: assigning %s the time of %s %s' % (program,t,unit))
    
                    if unit.lower() in ['sec','s']:
                        tsec=float(t)
                    elif unit.lower() in ['m','min','minutes','minute']:
                        tsec=float(t)*60.0
                    elif unit.lower() in ['h','hour','hours']:
                        tsec=float(t)*3600.0
                    else:
                        raise RuntimeError('Could not understand time unit of %s, only sec,min,hour allowed!')
                
                downtimerows.newrow({'expid':0,'twi':0,'utdate':utdown,'time':tsec,'program':downtimename,'Object':'DOWNTIME'})
                
                utdown = utdown + tsec*u.s
                if utdown.to_value('isot')>utdownmax.to_value('isot'):
                    raise RuntimeError('downtime UT %s (ut start + downtime) is later than the maximum %s!' % (utdown.to_value('isot'),utdownmax.to_value('isot')))
                
                #reset_index
                    
                #unit = downtime[i+2]
                counter_downtime += 1

        # Any downtime?
        if len(downtimerows.t)>0:
            
            if self.verbose:
                print('Downtime table:')
                downtimerows.write()
            # append the downtime table to qcinv table!
            self.qcinv.t = pd.concat([self.qcinv.t,downtimerows.t],ignore_index=True)        
             
            # only fill the new rows!
            ix_down = self.qcinv.ix_equal('Object','DOWNTIME')
            self.fill_qcinv_table(qcinv_indices = ix_down)
     
            #reset the index of the qcinv table, sorted by utdate!
            ixs = self.qcinv.ix_sort_by_cols(['utdate'])
            #self.qcinv.t.loc[ixs].reset_index(inplace=True)
            self.qcinv.t = self.qcinv.t.loc[ixs].reset_index(drop=True)
            
            if self.verbose:
                print('new qcinv table with Downtimes added!')
                self.qcinv.write()
        else:
            raise RuntimeError('Bug? there should be some downtime?')
        
        return(0)
        

    def fill_qcinv_table(self, qcinv_indices=None):
        """
        Fill out utdate, ut_decimal, and twi columns

        Parameters
        ----------
        qcinv_indices : TYPE, optional
            if specified, only passed indices are filled. The default is None.

        Returns
        -------
        None.

        """        

        if qcinv_indices is None:
            self.qcinv.t['utdate']=None
            self.qcinv.t['ut_decimal']=np.nan
            self.qcinv.t['twi']=0
            self.qcinv.default_formatters['ut_decimal']='{:.4f}'.format
            self.qcinv.t['program']=None
            self.qcinv.t['blockID']=0
           
        ix_all = self.qcinv.getindices(indices=qcinv_indices)
        
        # get just the dates (not hours) for the beginning and end of the night
        datestart = self.tonight[0].to_value('isot')[:10]
        dateend = self.tonight[1].to_value('isot')[:10]
        
        # reference t 
        t0 = time.Time(dateend+'T00:00:00.00',scale='utc',format='isot')
        
        
        for ix in ix_all:
            if self.qcinv.t.loc[ix,'utdate'] is None:
                # first try if the datestart is the correct date to use.
                tobs = time.Time(datestart+'T'+self.qcinv.t.loc[ix,'ut']+':00', scale='utc')
                dt = tobs - self.tonight[0]
                
                # if it is not after the start of the night, try dateend
                if dt.to_value('hr')<0.0:
                    tobs = time.Time(dateend+'T'+self.qcinv.t.loc[ix,'ut']+':00', scale='utc')
                    dt = tobs - self.tonight[0]
                    # comething is wrong!!
                    if dt.to_value('hr')<0.0:
                        raise RuntimeError('dt = %f, could not figure out the UT date for %s that is past the startdate %s' % (dt.to_value('hr'),self.qcinv.t.loc[ix,'ut'],self.tonight[0].to_value('isot')))
            
                # Make sure tobs is before the end of the night
                dt = self.tonight[1] - tobs
                if dt.to_value('hr')<0.0:
                    raise RuntimeError('dt = %f, could not figure out the UT date for %s that is before the enddate %s' % (dt.to_value('hr'),self.qcinv.t.loc[ix,'ut'],self.tonight[1].to_value('isot')))
                
                self.qcinv.t.loc[ix,'utdate']=tobs.to_value('isot')
            else:
                tobs = time.Time(self.qcinv.t.loc[ix,'utdate'], scale='utc')
            
            # now get the relative time in decimal hours with respect to t0
            dt = (tobs-t0)
            self.qcinv.t.loc[ix,'ut_decimal'] = dt.to_value('hr')
            
            # check for twilight!
            twi_zone = 0
            for i in range(len(self.horizons)):
                horizon = self.horizons[i]
                if (tobs-self.twi[horizon][0]).to_value('hr')<0.0:
                    twi_zone = i+1
                #print('vvv',(tobs-self.twi[horizon][1]).to_value('hr'),self.qcinv.t.loc[ix,'time']/3600.0)
                if (tobs-self.twi[horizon][1]).to_value('hr')+self.qcinv.t.loc[ix,'time']/3600.0>0.0:
                    twi_zone = i+1
            self.qcinv.t.loc[ix,'twi'] = twi_zone
                    
                
            #print(dt.to_value('hr'))
        if self.verbose>2: self.qcinv.write()
        return(0)

    
    def assignPrograms(self):
        self.create_fieldpattern2program()
        
        ixs = self.qcinv.getindices()
        
        blockID = 1
        
        special_programs = {}
        for p in ['UNKNOWN','TECHSETUP','STANDARDS']:
            special_programs[p]={}
            special_programs[p]['counter']=0
            special_programs[p]['pattern']=re.compile('^%s' % p)
        
        for i in range(len(ixs)):
            ix = ixs[i]

            # check for DOWNTIME, already prefilled!
            if (self.qcinv.t.loc[ix,'program'] is not None): 
                if re.search('^DOWNTIME',self.qcinv.t.loc[ix,'program']):
                    # all good!! Downtime!
                    pass
                else:
                    raise RuntimeError('Program preset to %s, but not DOWNTIME???' % self.qcinv.t.loc[ix,'program'])
            else:
                # found the program: check the search patterns for each Object
                foundflag=False
                for fieldpattern in self.fieldpattern2program:
                    m = self.fieldpattern2program[fieldpattern]['compiled']
                    if m.search(self.qcinv.t.loc[ix,'Object']):
                        program = self.fieldpattern2program[fieldpattern]['program']
                        if self.verbose>2: print('FOUND! pattern %s matches %s, program %s' % (fieldpattern,self.qcinv.t.loc[ix,'Object'],program))
                        if program=='2020B-0053_DEBASS':
                            if self.qcinv.t.loc[ix,'time']>25:
                                program='2021A-0275_YSE'
    
                        if program in special_programs:
                            m = special_programs[program]['pattern']
                            if i==0 or (not m.search(self.qcinv.t.loc[ixs[i-1],'program'])):
                                special_programs[program]['counter']+=1
                            program += '%d' %  special_programs[program]['counter']
     
                        self.qcinv.t.loc[ix,'program']=program
                        foundflag=True
                        break
                # Program not found? Assign to UNKNOWNX
                if (not foundflag): 
                    if (self.qcinv.t.loc[ix,'program'] is  None): 
                        program='UNKNOWN' 
                        m = special_programs[program]['pattern']
                        if i==0 or (not m.search(self.qcinv.t.loc[ixs[i-1],'program'])):
                            special_programs[program]['counter']+=1
                        program += '%d' %  special_programs[program]['counter']
                        self.addwarning('WARNING: Could not find the program for %s in line %d (block %s)' % (self.qcinv.t.loc[ix,'Object'],ix,program))
                        self.qcinv.t.loc[ix,'program']=program
             
            # New block ID?
            if i>0: 
                if self.qcinv.t.loc[ixs[i-1],'twi']!=self.qcinv.t.loc[ix,'twi'] or self.qcinv.t.loc[ixs[i-1],'program']!=self.qcinv.t.loc[ix,'program']:
                # if not the first entry AND if different than previous row's program: inc blockID
                #if ix!=ixs[0] and self.qcinv.t.loc[ixs[i-1],'program']!=self.qcinv.t.loc[ix,'program']:
                    blockID+=1
            self.qcinv.t.loc[ix,'blockID']=blockID     
        #sys.exit(0)
        return(0)
        
    def calcTimes(self):
        self.t['dt_block_h']=self.t['dt_prevgap_sec']=self.t['dt_nextgap_sec']=self.t['dt_gaps_sec']=self.t['dt_block_full_h']=self.t['dt_charged_h']=np.nan
        self.default_formatters['dt_block_h']='{:.4f}'.format
        self.default_formatters['dt_prevgap_sec']='{:.0f}'.format
        self.default_formatters['dt_nextgap_sec']='{:.0f}'.format
        self.default_formatters['dt_gaps_sec']='{:.0f}'.format
        self.default_formatters['dt_block_full_h']='{:.4f}'.format
        self.default_formatters['assigned_program']=self.programcol_formatter

        re_techsetup_standards_downtime = re.compile('^TECHSETUP|^STANDARDS|^DOWNTIME')
        re_downtime = re.compile('^DOWNTIME')
        blockIDs = self.qcinv.t['blockID'].unique()
        # get info for each block
        for i in range(len(blockIDs)):
            ixs = self.qcinv.ix_inrange('blockID', blockIDs[i],blockIDs[i])
            if len(ixs)==0:
                self.newrow({'blockID':blockIDs[i]})
                self.addwarning('WARNING: could not find any entries for blockID %d' % blockIDs[i])
                continue
            
            twi_zone = self.qcinv.t.loc[ixs,'twi'].unique()
            if len(twi_zone)!=1:
                raise RuntimeError('Could not determine twi_zone')
            
            # first get the difference between last and first UT
            dt_block_h = self.qcinv.t.loc[ixs[-1],'ut_decimal']-self.qcinv.t.loc[ixs[0],'ut_decimal'] 
            
            # exposure time and nominal overhead for last im
            #self.qcinv.write(indices=ixs)
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
                         'dt_nextgap_sec':dt_nextgap_sec*3600.0,                   
                         'twi':twi_zone[0]                      
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
            if re_downtime.search(self.t.loc[ixs_blocks[i],'assigned_program']):
                self.t.loc[ixs_blocks[i],'dt_nextgap_sec']=0.0
                self.t.loc[ixs_blocks[i],'dt_prevgap_sec']=0.0
            if re_techsetup_standards_downtime.search(self.t.loc[ixs_blocks[i],'assigned_program']):
                if i>0: self.t.loc[ixs_blocks[i-1],'dt_nextgap_sec']=0.0
                if i<len(ixs_blocks)-1: self.t.loc[ixs_blocks[i+1],'dt_prevgap_sec']=0.0
 
        
        # Add up block time and gap time:
        # TECHSETUP, STANDARDS: eat all gap time
        # all other programs: eat half of each gap time (previous and next)
        for i in range(len(ixs_blocks)):
            if re_techsetup_standards_downtime.search(self.t.loc[ixs_blocks[i],'assigned_program']):
                self.t.loc[ixs_blocks[i],'dt_gaps_sec'] = 1.0*(self.t.loc[ixs_blocks[i],'dt_nextgap_sec']+self.t.loc[ixs_blocks[i],'dt_prevgap_sec'])
            else:
                self.t.loc[ixs_blocks[i],'dt_gaps_sec'] = 0.5*(self.t.loc[ixs_blocks[i],'dt_nextgap_sec']+self.t.loc[ixs_blocks[i],'dt_prevgap_sec'])
            
            self.t.loc[ixs_blocks[i],'dt_block_full_h'] = self.t.loc[ixs_blocks[i],'dt_block_h']+self.t.loc[ixs_blocks[i],'dt_gaps_sec']/3600.0
            self.t.loc[ixs_blocks[i],'dt_charged_h'] = self.t.loc[ixs_blocks[i],'dt_block_full_h'] *  self.twi_charge_fraction[self.t.loc[ixs_blocks[i],'twi']]
                
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
        print('################\n### SUMMARY:')
        current_programs = self.t['assigned_program'].unique()
        allmainprograms = list(self.semesterinfo.programlist.keys())
        current_mainprograms = AandB(current_programs,allmainprograms)
        extra_programs = AnotB(current_programs,current_mainprograms)
        extra_programs.sort()
        current_mainprograms.sort()
        programs = list(current_mainprograms)
        programs.extend(extra_programs)
        for program in programs:
            ixs = self.ix_equal('assigned_program',program)
            ixs_dark = self.ix_equal('twi',0,indices=ixs)
            ixs_twi1  = self.ix_inrange('twi',1,1,indices=ixs)
            ixs_twi2  = self.ix_inrange('twi',2,2,indices=ixs)
            ixs_twi3  = self.ix_inrange('twi',3,3,indices=ixs)
            ixs_downtime  = self.ix_inrange('program','DOWNTIME0','DOWNTIME9999',indices=ixs)
            
            t_dark = self.t.loc[ixs_dark,'dt_charged_h'].sum()
            t_twi1 = self.t.loc[ixs_twi1,'dt_charged_h'].sum()
            t_twi2 = self.t.loc[ixs_twi2,'dt_charged_h'].sum()
            t_twi3 = self.t.loc[ixs_twi3,'dt_charged_h'].sum()
            t_downtime  = self.t.loc[ixs_downtime,'dt_charged_h'].sum()
            self.nightsummary.newrow({'assigned_program':program,
                                 't_dark':t_dark,
                                 't_twi1':t_twi1,
                                 't_twi2':t_twi2,
                                 't_twi3':t_twi3,
                                 't_down':t_downtime,
                                 })
        self.nightsummary.default_formatters['assigned_program']=self.programcol_formatter
        ixs_all = self.nightsummary.getindices()
        
    
        ix_used = self.nightsummary.newrow({'assigned_program':'total time used (hr)'})
        ix_available = self.nightsummary.newrow({'assigned_program':'total time avail. (hr)'})
        for timecol in ['t_dark','t_twi1','t_twi2','t_twi3','t_down']:
            ix_notnull = self.nightsummary.ix_remove_null(timecol,indices=ixs_all)
            self.nightsummary.t.loc[ix_used,timecol]=self.nightsummary.t.loc[ix_notnull,timecol].sum()
        
        self.nightsummary.t.loc[ix_available,'t_dark']=(self.twi[18][1]-self.twi[18][0]).to_value('sec')/3600.0
        self.nightsummary.t.loc[ix_available,'t_twi1']=((self.twi[18][0]-self.twi[15][0]).to_value('sec') + (self.twi[15][1]-self.twi[18][1]).to_value('sec'))/3600.0
        self.nightsummary.t.loc[ix_available,'t_twi2']=((self.twi[15][0]-self.twi[12][0]).to_value('sec') + (self.twi[12][1]-self.twi[15][1]).to_value('sec'))/3600.0

        self.nightsummary.write()
        print('Total dark time assigned : %.4f hours,\nTotal dark time available: %.4f hours\nDifference(available-assigned)=%.1f seconds' % (self.nightsummary.t.loc[ix_used,'t_dark'],self.nightsummary.t.loc[ix_available,'t_dark'],(self.nightsummary.t.loc[ix_available,'t_dark']-self.nightsummary.t.loc[ix_used,'t_dark'])*3600.0))
        
    def savetables(self,basename=None):
        if basename is not None:
            if len(basename)==0:
                basename=self.qcinv.filename
            elif len(basename)==1:
                basename=basename[0]
            else:
                raise RuntimeError('more than one argument for --save is not allowed!')
            print('Saving tables with basename',basename)
            self.qcinv.write(basename+'.times.txt',verbose=2)
            self.write(basename+'.blocks.txt',verbose=2)
            self.nightsummary.write(basename+'.nightsummary.txt',verbose=2)
            
            
        return(0)
    
    def setsemester(self,semester,semester_summaryfile=None):
        self.semesterinfo.setsemester(semester)
        self.semestersummary.semesterinfo = self.semesterinfo
        self.semester_summaryfile = self.semestersummary.get_semester_summary_filename(semester_summaryfile)
        self.semestersummary.initcols()
        if self.verbose>2: 
            print('Semester summary file:',self.semester_summaryfile)
        
    def load_semestersummary(self,filename=None):
        if filename is None:
            if self.semester_summaryfile is None:
                raise RuntimeError('Semester summary filename is not specified for loading!')
            filename = self.semester_summaryfile
        if os.path.isfile(filename):
            print('Loading',filename)
            self.semestersummary.load_spacesep(filename)
        else:
            print('WARNING! semester summary file %s does not exist yet!' % filename)
            return(1)
        return(0)

    def write_semestersummary(self,filename=None):
        if filename is None:
            if self.semester_summaryfile is None:
                raise RuntimeError('Semester summary filename is not specified for writing!')
            filename = self.semester_summaryfile
        print('Saving',filename)
        self.semestersummary.t['date']=self.semestersummary.t['date'].astype('int')
        ixs = self.semestersummary.ix_sort_by_cols(['date'])
        self.semestersummary.write(filename, indices=ixs, overwrite=True)
        return(0)
            
            
    def add2semestersummary(self,semester_summaryfile=None):
        self.semester_summaryfile = self.semestersummary.get_semester_summary_filename(semester_summaryfile)
        self.load_semestersummary(self.semester_summaryfile)
        
        # get the date for the semester summary entry
        date = re.sub('\.qcinv','',os.path.basename(self.qcinv.filename))
        print('Date:',date)
        if not re.search('^\d+$',date):
            raise RuntimeError('Could not determine the date from %s!' % self.qcinv.filename)
        
        # Add the results of this date to the semester summary
        self.semestersummary.add2summary(self.nightsummary,date)
        if self.verbose:
            print('### SEMESTER SUMMARY:')
            self.semestersummary.write()
            
        # save it!
        self.write_semestersummary(semester_summaryfile)
            

if __name__ == "__main__":
    calcTime = calcTimeclass()
    usagestring='USAGE: calcTime.py qcinv_filename'
    parser=calcTime.add_arguments(usage=usagestring)
    args = parser.parse_args()
    
    calcTime.verbose=args.verbose
    calcTime.debug=args.debug

    calcTime.setsemester(args.semester,args.semester_summaryfile)

    calcTime.readqcinv(args.qcinvfile)  
    calcTime.fill_qcinv_table()
    calcTime.downtime2qcinv(args.downtime)
    calcTime.assignPrograms()
    calcTime.calcTimes()
    calcTime.reassign_programs(args.reassign)
    calcTime.reassign_programs(calcTime.downtime_reassign)
    if args.verbose>1:
        calcTime.qcinv.write()
        calcTime.write()

    calcTime.mkSummary()
    
    calcTime.savetables(basename=args.save)
    if args.add2semestersummary:
        calcTime.add2semestersummary(args.semester_summaryfile)
        calcTime.semestersummary.summarystatistics()
        calcTime.semestersummary.showtables()
    
    if len(calcTime.warnings)>0:
        print('THERE WERE WARNINGS!!!')
        for s in calcTime.warnings: print(s)
