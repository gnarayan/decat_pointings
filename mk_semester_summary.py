#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 11:36:43 2021

@author: arest
"""
import sys, os, re, copy, shutil,io
import argparse
from pdastro import pdastroclass,AnotB,AandB
from semesterinfo import semesterinfoclass,default_semester
import numpy as np
from astropy import time
from astroplan import Observer
from astropy import units as u 
ctio = Observer.at_site("CTIO")

class semester_summary_class(pdastroclass):
    def __init__(self):
        pdastroclass.__init__(self)

        self.semesterinfo = semesterinfoclass()

        self.verbose=0
        self.debug=0
        
        self.timetypes = ['t_dark','t_twi1','t_twi2','t_down']
        self.cols4time = {}
        
    def initcols(self):
        self.t['night']=self.t['t_dark14']=self.t['t_dark_used']=self.t['t_unacc']=self.t['t_dark_avail']=None
        self.default_formatters['t_dark14']='{:.4f}'.format
        self.default_formatters['t_dark_used']='{:.4f}'.format
        self.default_formatters['t_unacc']='{:.4f}'.format
        self.default_formatters['t_dark_avail']='{:.4f}'.format
        programs = list(self.semesterinfo.programlist.keys())
        for timetype in self.timetypes:
            cols = list(map(lambda orig_string: orig_string + '_'+timetype, programs))
            self.cols4time[timetype]=cols
            self.t[cols]=np.nan
            for col in cols:
               self.default_formatters[col]='{:.4f}'.format
        

    def get_semester_summary_filename(self,semester_summaryfile=None):
        if semester_summaryfile is None:
            if self.semesterinfo.semester is None:
                raise RuntimeError('Semester is not defined yet!!!')
            semester_summaryfile = './%s/hours_summary.txt' % (self.semesterinfo.semester)
        return(semester_summaryfile)
            

    def add_arguments(self, parser=None, usage=None, conflict_handler='resolve'):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)
            
        parser.add_argument('-v', '--verbose', action='count', default=0)
        parser.add_argument('-d', '--debug', action='count', default=0)
#        parser.add_argument('--semester', default=default_semester, help='Specify the semester (default=%(default)s)')
#        parser.add_argument('--hours_summaryfile', default=default_hours_summaryfile, help='Specify the filename for the hours summary (default=%(default)s)')

        return(parser)
    
    def initnights(self,nights):

        for night in nights:
            self.initnight(night)
        self.t['night']=self.t['night'].astype('string')
        return(0)

    def initnight(self,night):
        night=str(night)

        # get the correct index from semester summary for this night.
        # If it doesn't exist yet, add a new row.
        ixs_semestersummary = self.ix_equal('night',night)
        if len(ixs_semestersummary)==0:
            # new row
            print('night %s is a new entry in semester summary' % night)
            ix_sem = self.newrow({'night':night})
        elif len(ixs_semestersummary)==1:
            ix_sem = ixs_semestersummary[0]
        else:
            raise RuntimeError('BUG! more than one entry for night %s' % night)

        night_iso = '%s-%s-%s' % (night[:4],night[4:6],night[6:8])        
        tonight = ctio.tonight(time.Time(night_iso, scale='utc',format='iso')+0.6,horizon=-14*u.deg)
        t_dark14 =(tonight[1]-tonight[0]).to_value('sec')/3600.0*self.semesterinfo.nights[night]
        self.t.loc[ix_sem,'t_dark14']=t_dark14
        
        return(ix_sem)

    def add2summary(self,nightsummary,night):
        # get the indices of the night summary of the main programs 
        ixs_nightsummary = []
        for program in self.semesterinfo.programlist.keys():
            ixs_nightsummary.extend(nightsummary.ix_equal('assigned_program',program))
            #ixs_nightsummary = nightsummary.ix_inrange('assigned_program','20','21')

        ix_sem = self.initnight(night)


        # get the correct index from semester summary for this night.
        # If it doesn't exist yet, add a new row.
        #ixs_semestersummary = self.ix_equal('night',int(night))
        #if len(ixs_semestersummary)==0:
            # new row
        #    print('night %s is a new entry in semester summary' % night)
       #     ix_sem = self.newrow({'night':night})
       # elif len(ixs_semestersummary)==1:
       #     ix_sem = ixs_semestersummary[0]
       # else:
       #     raise RuntimeError('BUG! more than one entry for night %s' % night)
            
        ix_tdark = nightsummary.ix_equal('assigned_program','total time used (hr)')
        ix_tdark_avail = nightsummary.ix_equal('assigned_program','total time avail. (hr)')
        self.t.loc[ix_sem,'t_dark_used'] = nightsummary.t.loc[ix_tdark[0],'t_dark']
        self.t.loc[ix_sem,'t_dark_avail'] = nightsummary.t.loc[ix_tdark_avail[0],'t_dark']
        
        cols =  list(nightsummary.t.loc[ixs_nightsummary,'assigned_program'])
        for timetype in self.timetypes:
            cols2 = list(map(lambda orig_string: orig_string + '_' + timetype, cols))
            self.t.loc[ix_sem,cols2]=list(nightsummary.t.loc[ixs_nightsummary,timetype])
        
        return(0)
    
    def summarystatistics(self):
        #n_nights = 0.0
        #n_nights_obs = 0.0
        #n_nights_left = 0.0
        #for night in self.t['night'].astype('string'):
        #    n_nights_obs += self.semesterinfo.nights[night]
            
        #for night in self.semesterinfo.nights:
        #    n_nights += self.semesterinfo.nights[night]
        #    if not(night in list(self.t['night'].astype('string'))):
        #        n_nights_left  += self.semesterinfo.nights[night]
 
        #print('\nTotal %.1f nights, %.1f already observed, %.1f left' % (n_nights,n_nights_obs,n_nights_left))
        
        self.t['t_unacc'] = self.t['t_dark14']-self.t['t_dark_used']
        
        ix_nights_observed = self.ix_remove_null(['t_dark_used'])
        ix_nights_not_observed = AnotB(self.getindices(),ix_nights_observed)
        n_nights = len(self.getindices())
        n_nights_obs = len(ix_nights_observed)
        n_nights_left = len(ix_nights_not_observed)
        
        t_sumprograms_h = 0.0
        for program in self.semesterinfo.programlist: t_sumprograms_h+=self.semesterinfo.programlist[program]
        print('\nTotal %d nights with observing time, %d already observed, %d left.' % (n_nights,n_nights_obs,n_nights_left))
        print('Total %.3f hours (%.3f hours assigned to programs), %.3f hours already observed, %.3f hours left.' % (self.t['t_dark14'].sum(),t_sumprograms_h,self.t.loc[ix_nights_observed,'t_dark14'].sum(),self.t.loc[ix_nights_not_observed,'t_dark14'].sum()))
        print('%.3f hours already observed: %.3f hours accounted for, %.3f hours not accounted for.' % (self.t.loc[ix_nights_observed,'t_dark14'].sum(),self.t.loc[ix_nights_observed,'t_dark_used'].sum(),
                                                                                                        self.t.loc[ix_nights_observed,'t_unacc'].sum()))        

        # Calculate the stats for each program
        ix_tot_used = self.newrow({'night':'total time used'})
        ix_tot = self.newrow({'night':'total time'})
        ix_tot_left = self.newrow({'night':'total time left'})
        ix_excess = self.newrow({'night':'total time over-usage'})
        ix_tot_used_pernight = self.newrow({'night':'time used per night'})
        ix_tot_left_pernight = self.newrow({'night':'time left per night'})
        #ix_timeleft = self.newrow({'night':'time left'})
        for timetype in self.timetypes:
            for col in self.cols4time[timetype]:
                self.t.loc[ix_tot_used,col] = self.t[col].sum()

        timetype = 't_dark'
        for col in self.cols4time[timetype]:
            self.t.loc[ix_tot,col] = self.semesterinfo.programlist[re.sub('_'+timetype,'',col)]
            self.t.loc[ix_tot_left,col] = self.t.loc[ix_tot,col] - self.t.loc[ix_tot_used,col]
            self.t.loc[ix_tot_used_pernight,col] = self.t.loc[ix_tot_used,col]/n_nights_obs
            self.t.loc[ix_tot_left_pernight,col] = self.t.loc[ix_tot_left,col]/n_nights_left
            self.t.loc[ix_excess,col] = self.t.loc[ix_tot_used,col]-self.t.loc[ix_tot,col]*n_nights_obs/n_nights
 

    def showtables(self,verbose=0):
        print('\n### t_dark',verbose)
        self.write(columns = ['night','t_dark14','t_dark_used','t_unacc']+self.cols4time['t_dark'])
        if verbose>2:
            for timetype in ['t_twi1','t_twi2','t_down']:
                print('\n### %s' % timetype)
                self.write(columns = ['night','t_dark14','t_dark_used','t_unacc']+self.cols4time[timetype])
        
        
if __name__ == "__main__":
    summary = semester_summary_class()
    usagestring='USAGE: mk_semester_summary.py'
    parser=summary.add_arguments(usage=usagestring)
    args = parser.parse_args()
    
    summary.verbose=args.verbose
    summary.debug=args.debug

