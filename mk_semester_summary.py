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

class semester_summary_class(pdastroclass):
    def __init__(self):
        pdastroclass.__init__(self)

        self.semesterinfo = semesterinfoclass()
        self.timecols = ['t_dark','t_twi1','t_twi2','t_down']

        self.verbose=0
        self.debug=0
        
    def initcols(self):
        self.t['date']=None
        programs = list(self.semesterinfo.programlist.keys())
        for timecol in self.timecols:
            cols = list(map(lambda orig_string: orig_string + '_'+timecol, programs))
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

    def add2summary(self,nightsummary,date):
        print(self.t)
        print(self.t.info())
        
        #
        ixs_nightsummary = []
        for program in self.semesterinfo.programlist.keys():
            ixs_nightsummary.extend(nightsummary.ix_equal('assigned_program',program))
            #ixs_nightsummary = nightsummary.ix_inrange('assigned_program','20','21')


        # get the correct index from semester summary for this date.
        # If it doesn't exist yet, add a new row.
        ixs_semestersummary = self.ix_equal('date',int(date))
        if len(ixs_semestersummary)==0:
            # new row
            print('date %s is a new entry in semester summary' % date)
            ix_sem = self.newrow({'date':date})
        elif len(ixs_semestersummary)==1:
            ix_sem = ixs_semestersummary[0]
        else:
            raise RuntimeError('BUG! more than one entry for date %s' % date)
            
        cols =  list(nightsummary.t.loc[ixs_nightsummary,'assigned_program'])
        for timecol in self.timecols:
            #bla = ['_'+timecol]*len(cols)
            #print('cccc',bla)
            #cols2 = cols + bla
            cols2 = list(map(lambda orig_string: orig_string + '_' + timecol, cols))
            self.t.loc[ix_sem,cols2]=list(nightsummary.t.loc[ixs_nightsummary,timecol])
        
        return(0)

    
if __name__ == "__main__":
    summary = semester_summary_class()
    usagestring='USAGE: mk_semester_summary.py'
    parser=summary.add_arguments(usage=usagestring)
    args = parser.parse_args()
    
    summary.verbose=args.verbose
    summary.debug=args.debug

