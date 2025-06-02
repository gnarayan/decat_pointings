#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 19:44:17 2024

@author: arest
"""

import argparse,os,sys,re,copy,shutil,glob

from astropy.time import Time, TimeDelta
from datetime import datetime
import pandas as pd
import numpy as np
from subprocess import Popen, PIPE, STDOUT

import pytz

import yaml

from pdastro import pdastroclass,AnotB,unique

from astroplan import Observer,moon
from astropy import units as u
from pytz import timezone
import matplotlib.pyplot as plt

from airmass import AirmassCalculator,ctio_location
from moon import MoonCalculator

from astropy.coordinates import (EarthLocation, SkyCoord, AltAz, get_sun,
                                 get_moon, Angle, Longitude)


def calc_slewtimes_2arrays(coord1,coord2,decflip_time=180.0):
    """
    Calculate the slew times for a 2 sets of SkyCoor (The slew time does NOT include the overhead!, so for short separations, slew time=0.0!):
        slew_time[i] is the time needed to slew from coord1[i] to coord2[i]

    Parameters
    ----------
    coord1 : SkyCoor list
        The first skycoor list.
    coord2 : SkyCoor list
        The second skycoor list.
    decflip_time: float, optional
        Extra time needed to do the declination flip. The default is 180.0.
    Returns
    -------
    slew times, Delination flip, and separation

    """
    # get the separation in degree
    sep = coord1.separation(coord2).to(u.deg).value
    
    slew_rate = 5 / 11.0
    # slew_time does NOT include the overhead!, so for short separations, slew time=0.0!
    slew_time = sep / slew_rate

    # readout is about 29 s, and slew time graph has y-intercept of ~29 s
    # So 29 s is added either way
    # slew_time = sep / slew_rate + 29.0

    # Determine if there is a Declination flip at -30.17 from coord1[i] to coord2[i] 
    #south1 = np.where(coord1.dec < -30.17 * u.deg,1,0)
    #south2 = np.where(coord2.dec < -30.17 * u.deg,1,0)
    south1 = np.where(coord1.dec < -30.8 * u.deg,1,0)
    south2 = np.where(coord2.dec < -30.8 * u.deg,1,0)
    decflip = np.bitwise_xor(south1,south2)
    
    # if declination flip, add decflip_time to slew_time
    slew_time += decflip * decflip_time

    return(slew_time,decflip,sep)
    
def calc_slewtimes_1array(coord,decflip_time=180.0):
    """
    Calculate the slew times for a 1 sets of SkyCoor (The slew time does NOT include the overhead!, so for short separations, slew time=0.0!):
        slew_time[i] is the time needed to slew from coord[i-1] to coord[i]
    

    Parameters
    ----------
    coord : SkyCoor list
    decflip_time : float, optional
        Extra time needed to do the declination flip. The default is 180.0.

    Returns
    -------
    slew times, Delination flip, and separation

    """
    # create the two sets of coordinates: shift by 1!
    coord1 = coord[:-1]
    coord2 = coord[1:]
    
    # calc the slew times, Delination flip, and separation
    (slew_time,decflip,sep)=calc_slewtimes_2arrays(coord1,coord2,decflip_time=decflip_time)
    
    # insert 0s at the beginning of the lists
    slew_time = np.insert(slew_time,0,0.0)
    decflip = np.insert(decflip,0,0)
    sep = np.insert(sep,0,0.0)
    
    return(slew_time,decflip,sep)


class obsplan_baseclass:

    def get_arguments(self, args, configfile = None):
        '''

        Parameters
        ----------
        args : list
            pass the command line arguments to self.params.
        configfile : string, optional
            Config filename. The default is None. If None, then
            $SYNDIFF_CFGFILE is used if exists.

        Returns
        -------
        None.

        '''

        def subenvvarplaceholder(paramsdict):
            """ Loop through all string parameters and substitute environment variables. environment variables have the form $XYZ """
            envvarpattern=re.compile('\$(\w+)')

            for param in paramsdict:
                print(param)
                if isinstance(paramsdict[param], str):
                    envvarnames=envvarpattern.findall(paramsdict[param])
                    if envvarnames:
                        for name in envvarnames:
                            if not (name in os.environ):
                                raise RuntimeError("environment variable %s used in config file, but not set!" % name)
                            else:
                                envval=os.environ[name]
                                subpattern='\$%s' % (name)
                                paramsdict[param] = re.sub(subpattern,envval,paramsdict[param])
                elif isinstance(paramsdict[param], dict):
                #elif type(dict[param]) is types.DictType:
                    # recursive: sub environment variables down the dictiionary
                    subenvvarplaceholder(paramsdict[param])
            return(0)


        # get the parameters from the config file
        if args.configfile is not None:
            #cfgparams = yaml.load_file(args.configfile)
            if not os.path.isfile(args.configfile):
                raise RuntimeError(f'config file {args.configfile} does not exist!')
            print(f'Loading config file {args.configfile}')
            cfgparams = yaml.load(open(args.configfile,'r'), Loader=yaml.FullLoader)
            self.params.update(cfgparams)

            subenvvarplaceholder(self.params)

            if args.verbose>2:
                print('\n### CONFIG FILE PARAMETERS:')
                for p in cfgparams:
                    print('config file args: setting %s to' % (p),cfgparams[p])

        # Go through optional parameters.
        # 'None' does not overwrite previously set parameters (either default or config file)
        if args.verbose>2:
            print('\n### OPTIONAL COMMAND LINE PARAMETERS (shown because verbose>2):')

        argsdict = vars(args)
        for arg in argsdict:

            # skip config file
            if arg=='configfile': continue

            if argsdict[arg] is not None:
                if args.verbose>2:
                    print('optional args: setting %s to %s' % (arg,argsdict[arg]))
                self.params[arg]=argsdict[arg]
            else:
                if arg not in  self.params:
                    self.params[arg]=None

        if 'verbose' in self.params:
            self.verbose = self.params['verbose']

        if self.verbose>2:
            print('\n### FINAL PARAMETERS (shown because verbose>2):')
            for p in self.params:
                print(p,self.params[p])

        return(0)
    
    def get_param_list(self,paramval):
        if paramval is None:
            return([])
        elif isinstance(paramval,tuple) or isinstance(paramval,list):
            return(paramval)
        elif isinstance(paramval,str):
            if re.search('\s+|\,',paramval):
                paramlist = re.split('\s+|\,',paramval)
                return(paramlist)
            else:
                return([paramval])
        else:
            return([paramval])

    
    def __init__(self):
        self.verbose=0
        self.ordercounter_stepsize = 10

        self.reset()

    def reset(self):

        self.programtable = pdastroclass()
        self.json_ordertable = pdastroclass()

        # self.params will be populated with the arguments
        self.params = {}
                
        self.YYMMDD=None
        self.semester=None

        # main directories
        self.json_dir = None
        self.obsplan_filename = None
        self.obsplan_dir = None
        
        self.jsontable = pdastroclass(columns=['jsonID','program','priority','json_short','ra','dec','UT','CLT','order','t_exp[m]','t_slew[m]','t_tot[m]','sep_deg','decflip','airmass','sep_moon','json_filename'])
                                               #'before_UT','after_UT','UT1_airmass2','UT2_airmass2','json_filename'])
        
        self.exposures_columns2copy=['expID','object','ra','dec','filter','exptype','count','exptime','proposer','propid']
        self.exposurestable = pdastroclass(columns=self.exposures_columns2copy)

        self.programtable = None

        self.ExtraTime_min = None

    def define_optional_arguments(self,parser=None,usage=None,conflict_handler='resolve'):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)

        if not('OBSPLAN_ROOTDIR' in os.environ):
            raise RuntimeError('The environment variable OBSPLAN_ROOTDIR needs to be defined!')
        obsplan_rootdir = os.environ["OBSPLAN_ROOTDIR"]

        defaultcfgfilename = f'{os.environ["OBSPLAN_ROOTDIR"]}/obsplan.cfg'
        #defaultjsonrootdir = f'{os.environ["OBSPLAN_ROOTDIR"]}/json_files'
        #defaultobsplanrootdir  = f'{os.environ["OBSPLAN_ROOTDIR"]}/obsplans'

        #parser.add_argument('-s','--sector', type=int, default=None, help='TESS Sector.')
        #parser.add_argument('--ccd', type=int, default=None, choices=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16],help='TESS CCD.')
        #parser.add_argument('-o','--outsubdir', type=str, default=None,help='subdir that gets added to the output directory')

        parser.add_argument('-v','--verbose', default=0, action='count')
        parser.add_argument('-c','--configfile', type=str, default=defaultcfgfilename, help=f'primary config file. default is set to {defaultcfgfilename}. Use -vvv to see the full list of all set parameters.')
        #parser.add_argument('--jsonrootdir', type=str, default=defaultjsonrootdir, help=f'root directory for the json files: <jsonrootdir>/<semester>/<programname>_<progID>/YYMMD; default is {defaultjsonrootdir}')
        parser.add_argument('--obsplanrootdir', type=str, default=obsplan_rootdir, help='root directory for the obsplan files: <obsplanrootdir>/<semester>/YYMMD/YYMMDD_obsplan.txt (default=%(default)s)')
        parser.add_argument('-d','--date', type=str, default='now', help='date in YYMMDD format. If "now", then the current date is used (default=%(default)s)')
        parser.add_argument('--add2starttime_minutes', type=float, default=None, help='Add this amount of minutes to the start time (default=%(default)s)')
        parser.add_argument('-1','--firsthalf', default=False, action='store_true', help='first half observing')        
        parser.add_argument('-2','--secondhalf', default=False, action='store_true', help='second half observing')        
        parser.add_argument('-p','--makeplots', default=False, action='store_true', help='make the airmass and moon plots and save them as *airmass.png file and *moon.png, respectively (same basename as obsplan)')        
        parser.add_argument('-s','--showplots', default=False, action='store_true', help='make the airmass and moon plots, save it as *airmass.png file and *moon.png, respectively (same basename as obsplan), and also pop it up.')        
        return(parser)
    
    def initialize4date(self,YYMMDDdate):
        # determine the date in YYMMDD format
        if YYMMDDdate.lower()=='now':
            t = Time.now()
            m = re.search('^\d\d(\d\d)\-(\d\d)\-(\d\d)T',t.to_value('isot'))
            if m is not None:
                self.YYMMDD = ''.join(m.groups())
            else:
                raise RuntimeError('Could not parse {} for YYMMDD'.format(t.to_value('isot')))
        else:
            self.YYMMDD = YYMMDDdate
            if len(self.YYMMDD)!=6:
                raise RuntimeError(f'"{self.YYMMDD}" should have 6 numbers, but has a length of {len(self.YYMMDD)}!')
            m = re.search('^\d+$',self.YYMMDD)
            if m is None:
                raise RuntimeError(f'"{self.YYMMDD}" should only have numbers!')
        
        # get the semester
        YY = int(self.YYMMDD[:2])
        MM = int(self.YYMMDD[2:4])
        if MM==1:
            YY-=1
            MM=13
        MM-=1
        self.semester = f'20{YY}'
        if MM<7:
            self.semester+='A'
        else:
            self.semester+='B'
        
        print(f'YYMMDD: {self.YYMMDD}')
        print(f'semester: {self.semester}')
        
        date = f'20{self.YYMMDD[:2]}-{self.YYMMDD[2:4]}-{self.YYMMDD[4:6]}T18:00:00'

        t = Time(date)        

        self.observer = Observer.at_site(self.params['observatory'],timezone=self.params['timezone'])
        self.horizon = 14

        print('\n##### night info:')
        #twi =  self.observer.tonight(t+0.5,horizon=-horizon*u.deg)
        twi =  self.observer.tonight(t,horizon=-self.horizon*u.deg)
        #self.UT_eveningtwi =  twi[0].to_value('isot')
        #self.UT_morningtwi =  twi[1].to_value('isot')
        self.UT_eveningtwi =  twi[0]
        self.UT_morningtwi =  twi[1]
        self.UT_night_midpoint = twi[0]+0.5*(twi[1]-twi[0])
        self.UT_eveningtwi.format = "iso"
        self.UT_morningtwi.format = "iso"

        self.get_night_info_description()
        
        # define directories and filenames
        self.json_dir = f'{self.params["obsplanrootdir"]}/{self.params["jsonfiles_subdir"]}/{self.semester}'
        programtablename = f'{self.json_dir}/{self.semester}_programs.txt'
        generictargetorder_filename = f'{self.json_dir}/{self.semester}_generictargetorder.txt'
        
        #self.obsplan_dir = f'{self.params["obsplanrootdir"]}/{self.params["obsplan_subdir"]}/{self.semester}/{self.YYMMDD}'
        self.obsplan_dir = f'{self.params["obsplanrootdir"]}'
        if self.params["obsplan_subdir"] is not None: self.obsplan_dir += f'/{self.params["obsplan_subdir"]}'
        self.obsplan_dir += f'/{self.semester}/{self.YYMMDD}'
        
        self.obsplan_filename = f'{self.obsplan_dir}/{self.YYMMDD}_obsplan.txt'
        self.obsplan4human_filename = f'{self.obsplan_dir}/{self.YYMMDD}_obsplan4human.txt'
        self.airmass_filename = f'{self.obsplan_dir}/{self.YYMMDD}_airmass.png'
        self.moon_filename = f'{self.obsplan_dir}/{self.YYMMDD}_moon.png'
        self.alltargets = f'{self.obsplan_dir}/{self.YYMMDD}_alltargets.txt'
        self.airmass_cache = f'{self.params["obsplanrootdir"]}/{self.params["obsplan_subdir"]}/airmass_grid'
        print(f'json_dir: {self.json_dir}')
        print(f'obsplan_dir: {self.obsplan_dir}')
        print(f'obsplan: {self.obsplan_filename}')
        print(f'obsplan4human: {self.obsplan4human_filename}')
        print(f'alltargets: {self.alltargets}')

        print(f'\n### programtablename: {programtablename}')
        self.programtable = pdastroclass()
        self.programtable.load(programtablename)
        if self.verbose>1: self.programtable.write()
    
        print(f'\n### generic target order filename: {generictargetorder_filename}')
        self.generictargetorder = pdastroclass()
        self.generictargetorder.load(generictargetorder_filename)
        if self.verbose>1: self.generictargetorder.write()
        
    def get_night_info_description(self):
        output = ''
        output += f'UT {-self.horizon} deg twilight: {self.UT_eveningtwi.to_value("isot")} {self.UT_morningtwi.to_value("isot")}\n' #% (twi[0].to_value('isot'),twi[1].to_value('isot'))
        output += f'UT middle of the night: {self.UT_night_midpoint.to_value("isot")}\n\n'# % (self.UT_night_midpoint.to_value('isot')))

        tznames = self.params['timezones'] 
        tz2use = self.params['timezones2use'] 
        output_midpoint = ''
        for tz_abbr in tz2use:
            output += f'timezone {tz_abbr}: {tznames[tz_abbr]}'
            tz = timezone(tznames[tz_abbr])
            eveningtwi = self.UT_eveningtwi.to_datetime(timezone=tz)
            morningtwi = self.UT_morningtwi.to_datetime(timezone=tz)
            night_midpoint = eveningtwi+0.5*(morningtwi-eveningtwi)
            output += f'{tz_abbr:3s} {-self.horizon} deg twilight: {eveningtwi.strftime("%m/%d/%Y,   %H:%M:%S")} to {morningtwi.strftime("%H:%M:%S")}\n' #% (twi[0].to_value('isot'),twi[1].to_value('isot'))
            output_midpoint +=  f'{tz_abbr:3s} middle of the night: {night_midpoint.strftime("%H:%M:%S")}\n'
        output+=output_midpoint

        return(output)
            

    
    def jsondir4program(self,programname,YYMMDD=None,priority=None):
        if YYMMDD is None:
            YYMMDD = self.YYMMDD
        jsondir4program = f'{self.json_dir}/{programname}/{YYMMDD}'
        if (priority is not None) and (priority!=0):
            jsondir4program += f'/{priority}'
        return(jsondir4program)
            
    def find_new_jsonfiles(self, programnames=None, YYMMDD=None, priorities=None, 
                           # delete json files entries in the jsontable that do not exist anymore!
                           delete_if_not_exists=True):
        if programnames is None:
            programnames = list(self.programtable.t['name'])
        # if programnames is a single string or number, make it a list!
        if isinstance(programnames,str):
            programnames=[programnames]
        if priorities is None:
            priorities=[0,1,2,3,4,5]
        else:
            # if priorities is a single string or number, make it a list!
            if isinstance(priorities,str) or isinstance(priorities,int):
                priorities=[int(priorities)]
            
        ixs_new = []
            
        print('bbb',programnames)
        for programname in programnames:
            
            # these are the indices for  all the json files already in the table for the given program
            ixs_program = self.jsontable.ix_equal('program',programname)
            
            for priority in priorities:
                
                # these are the indices for  all the json files already in the table for the given program and priority
                ixs_prio = self.jsontable.ix_equal('priority',priority,indices=ixs_program)
                jsonfiles_already_in_table = self.jsontable.t.loc[ixs_prio,'json_short']
                
                jsondir4program = self.jsondir4program(programname,YYMMDD=YYMMDD,priority=priority)
                if not os.path.isdir(jsondir4program):
                    if self.verbose>2: print(f'Skipping {jsondir4program} since doesn\'t exist')
                    continue
                    
                
                if self.verbose>2: 
                    print(f'Looking for json files in {jsondir4program}, priority {priority}   {jsondir4program}/*.json')
                filenames = glob.glob(f'{jsondir4program}/*.json')
                if len(filenames)==0:
                    if self.verbose>2: 
                        print(f'No json files for program {jsondir4program}, priority {priority}')
                if len(filenames)>0:
                    tmp = pdastroclass(columns=['json_filename','json_short'])
                    tmp.t['json_filename'] = filenames
                        
                    tmp.replace_regex('json_filename','json_short','.*\/','')
                    #tmp.write()
    
                    print(f'jsonfiles_already_in_table: {jsonfiles_already_in_table}')
                    newfiles = AnotB(list(tmp.t['json_short']),list(jsonfiles_already_in_table))
                    print(f'newfiles: {newfiles}')
                    files_do_not_exist = AnotB(list(jsonfiles_already_in_table),list(tmp.t['json_short']))
                    print(f'files_do_not_exist: {files_do_not_exist}')
                    
                    if self.verbose: 
                        print(f'{programname}: {len(newfiles)} new json files, {len(files_do_not_exist)} json files that don\'t exist anymore')
                        
                    for filename in newfiles:
                        ix_tmp = int(tmp.t[tmp.t['json_short']==filename].index.values)
                        ix = self.jsontable.newrow({'program':programname,
                                                    'priority':priority,
                                                    'json_short':filename,
                                                    'json_filename':tmp.t.loc[ix_tmp,'json_filename']})
                        ixs_new.append(ix)
                else:
                    pass
                
        if len(ixs_new)>0:
            self.jsontable.t.loc[ixs_new,['order','decflip']]=(-1,0)
            
            #####   ADD OFFSET HERE!!!
            ixs_notnull = self.jsontable.ix_not_null('jsonID')
            if len(ixs_notnull)==0:
                jsonID_start = 0
            else:
                jsonID_start = int(np.max(self.jsontable.t.loc[ixs_notnull,'jsonID']))+1
            
            #self.jsontable.t.loc[ixs_new,'jsonID']=range(len(ixs_new))
            self.jsontable.t.loc[ixs_new,'jsonID']=range(jsonID_start,jsonID_start+len(ixs_new))

        self.jsontable.t['jsonID'] = self.jsontable.t['jsonID'].astype(int)


        if self.verbose>2:
            self.jsontable.write()
            
        return(ixs_new)
                #for filename in filenames:
                #    if self.verbose>1: 
                #        print(f'{programname}: {os.path.basename(filename)}')
            

    def update_tables(self, ixs_json = None, programs = []):
        ixs_json = self.jsontable.getindices(indices=ixs_json)
        
        # if only for limited list of programs: get the ixs_jsons!
        if (programs is not None) and len(programs)>0:
            ixs_json = []
            for program in programs:
                ixs_json.extend(self.jsontable.ix_equal('program',program))
        
        ixs_delete = []
        
        for ix_json in ixs_json:
            filename = self.jsontable.t.loc[ix_json,'json_filename']

            # Load the json file and change the relevant columns into lower case
            exptable = pdastroclass()
            if not os.path.isfile(filename):
                print(f'WARNING: json file {filename} does not exist, markoing it for deletion!')
                ixs_delete.append(ix_json)
                continue
            print(f'Updating json file {filename}')
            exptable.t = pd.read_json(filename)
            columns2copy = []
            for col in exptable.t.columns:
                if col.lower() in self.exposures_columns2copy:
                    if col.lower()!=col:
                        exptable.t=exptable.t.rename(columns={col:col.lower()})
                    columns2copy.append(col.lower())
            
            # make sure the essential columns are there
            for col in ['ra','dec','count','exptime']:
                if col not in columns2copy: raise RuntimeError(f'"{col}" not in {columns2copy}')
            
            ixs_notnull_exp = exptable.ix_not_null('ra')
            ixs_notnull_exp = exptable.ix_not_null('dec',indices=ixs_notnull_exp)
            
            exptable.assert_radec_cols_decimal_degrees(racol='ra',deccol='dec',indices=ixs_notnull_exp)
            
            exptable.t['t_exp[m]'] = (exptable.t['exptime']+np.full(len(exptable.t['exptime']),float(self.params['overhead_sec']))) * exptable.t['count']/60.0
            columns2copy.append('t_exp[m]')
            
            #print(exptable.t)
            self.jsontable.t.loc[ix_json,'t_exp[m]'] = float(np.sum(exptable.t['t_exp[m]']))
            self.jsontable.t.loc[ix_json,'ra'] = np.mean(exptable.t.loc[ixs_notnull_exp,'ra'])
            self.jsontable.t.loc[ix_json,'dec'] = np.mean(exptable.t.loc[ixs_notnull_exp,'dec'])
            
            self.exposurestable.t = pd.concat([self.exposurestable.t,exptable.t[columns2copy]],axis=0, ignore_index=True)
        
        if len(ixs_delete)>0:
            print('WARNING: the following json files do not seem to exist anymore!!!')
            self.jsontable.write(indices=ixs_delete)
            do_it = input('Do you want to continue and delete the above entries? [y/n]?  ')
            if do_it.lower() in ['y','yes']:
                self.jsontable.t = self.jsontable.t.drop(ixs_delete,axis=0)
            elif do_it.lower() in ['n','no']:
                print('OK, stopping....')
                sys.exit(0)
            else:
                print(f'Hmm, \'{do_it}\' is neither yes or no. Don\'t know what to do, so stopping ....')
                sys.exit(0)

        self.jsontable.t['order'] = self.jsontable.t['order'].astype(int)
        self.jsontable.t['ra'] = self.jsontable.t['ra'].astype(float)
        self.jsontable.t['dec'] = self.jsontable.t['dec'].astype(float)
        self.jsontable.t['t_exp[m]'] = self.jsontable.t['t_exp[m]'].astype(float)
        self.jsontable.t['t_slew[m]'] = self.jsontable.t['t_slew[m]'].astype(float)
        self.jsontable.t['t_tot[m]'] = self.jsontable.t['t_tot[m]'].astype(float)
        self.jsontable.t['sep_deg'] = self.jsontable.t['dec'].astype(float)
        #self.jsontable.radeccols_to_coord('ra','dec','coord')
        self.jsontable.default_formatters = {}
        self.jsontable.default_formatters['t_exp[m]']='{:,.2f}'.format
        self.jsontable.default_formatters['ra']='{:,.2f}'.format
        self.jsontable.default_formatters['dec']='{:,.2f}'.format
        self.jsontable.default_formatters['t_slew[m]']='{:,.2f}'.format
        self.jsontable.default_formatters['t_tot[m]']='{:,.2f}'.format
        self.jsontable.default_formatters['sep_deg']='{:,.2f}'.format
        #self.jsontable.t.style.format({'t_exp[m]':'{:.2f}'})
        if self.verbose>3:
            self.jsontable.write(columns=['program','json_short','ra','dec','priority','sep_deg','t_slew[m]','t_exp[m]','t_tot[m]'])
        
        
    def wrap_up_and_summary(self):

        (ixs, ixs_ordered, ixs_notordered) = self.final_ixs()
        cols = self.get_param_list(self.params['obsplan_columns_short'])
        
        self.jsontable.t.loc[ixs_notordered,'UT']=np.nan
        for tz in self.params['timezones2use']:
            self.jsontable.t.loc[ixs_notordered,tz]=np.nan
        self.jsontable.t.loc[ixs_notordered,'t_slew[m]']=np.nan
        self.jsontable.t.loc[ixs_notordered,'t_tot[m]']=np.nan
        self.jsontable.t.loc[ixs_notordered,'sep_deg']=np.nan
        self.jsontable.t.loc[ixs_notordered,'airmass']=np.nan
        self.jsontable.t.loc[ixs_notordered,'sep_moon']=np.nan
        
        
        self.jsontable.write(indices=ixs,columns=cols)
        
        
        self.programtable.t['t_prio0[m]']=np.NaN 
        self.programtable.t['t_prioX[m]']=np.NaN 
        
        programs=unique(self.jsontable.t.loc[ixs_ordered,'program'])
        
        sum_tot=0.0
        ixs_prio0 = self.jsontable.ix_equal('priority',0,indices=ixs_ordered)
        for program in programs:
            ixs = self.jsontable.ix_equal('program',program,indices=ixs_prio0)
            sum_program = np.sum(self.jsontable.t.loc[ixs,'t_tot[m]'])
            ix_program=self.programtable.ix_equal('name',program)
            if len(ix_program)<1:
                raise RuntimeError(f'program {program} not found in program table!')
            self.programtable.t.loc[ix_program[0],'t_prio0[m]']=sum_program
            if program!='twi':
                sum_tot += sum_program
            #print(f'{program}: {sum_program:.2f} min')
        print(f'total priority 0: {sum_tot/60.0:.2f} h')
        
        sum_tot=0.0
        ixs_prioX = self.jsontable.ix_inrange('priority',1,None,indices=ixs_ordered)
        for program in programs:
            ixs = self.jsontable.ix_equal('program',program,indices=ixs_prioX)
            sum_program = np.sum(self.jsontable.t.loc[ixs,'t_tot[m]'])
            ix_program=self.programtable.ix_equal('name',program)
            if len(ix_program)<1:
                raise RuntimeError(f'program {program} not found in program table!')
            self.programtable.t.loc[ix_program[0],'t_prioX[m]']=sum_program
            if program!='twi':
                sum_tot += sum_program
            #print(f'{program}: {sum_program:.2f} min')
        print(f'total priority 1+: {sum_tot/60.0:.2f} h')
        
        #self.programtable.t['']
        self.programtable.default_formatters['t_prio0[m]']='{:,.1f}'.format
        self.programtable.default_formatters['t_prioX[m]']='{:,.1f}'.format
        self.programtable.write()
        
        print(f'Extra Time left at the end of the night: {self.ExtraTime_min:.2f} minutes')
        if self.ExtraTime_min<0.0:
            print(f'\n####################################\n### WARNING!! OBSPLAN IS {-self.ExtraTime_min:.2f} minutes TOO LONG!!!\n####################################')


    def update_json_tables_prog4(self, programname, YYMMDD=None, priority=None,
                           Nmax=None, skip_if_exists=True, reload_targets=False):
        jsondir4program = self.jsondir4program(programname,YYMMDD=None,priority=None)
        if not os.path.isdir(jsondir4program):
            print(f'Skipping {jsondir4program} since doesn\'t exist')
            return(0)
            
        if self.verbose>2: 
            print(f'Looking for json files in {jsondir4program}')
        filenames = glob.glob(f'{jsondir4program}/*.json')
        for filename in filenames:
            if self.verbose>1: 
                print(f'{programname}: {os.path.basename(filename)}')
                
            # Load the json file and change the relevant columns into lower case
            jsontable = pdastroclass()
            jsontable.t = pd.read_json(filename)
            columns2copy = []
            for col in jsontable.t.columns:
                if col.lower() in self.exposures_columns2copy:
                    if col.lower()!=col:
                        jsontable.t=jsontable.t.rename(columns={col:col.lower()})
                    columns2copy.append(col.lower())
            
            # make sure the essential columns are there
            for col in ['ra','dec','count','exptime']:
                if col not in columns2copy: raise RuntimeError(f'"{col}" not in {columns2copy}')
                
            jsontable.assert_radec_cols_decimal_degrees(racol='ra',deccol='dec')
                    
            print(jsontable.t)
            self.exposurestable.t = pd.concat([self.exposurestable.t,jsontable.t[columns2copy]],axis=0, ignore_index=True)
            #sys.exit(0)

        self.exposurestable.default_formatters['ra']='{:.3f}'.format
        self.exposurestable.default_formatters['dec']='{:.3f}'.format
        if self.verbose>3:
            self.exposurestable.write()

        return(0)             

    def update_json_tables4(self, programnames=None, YYMMDD=None, priority=None,
                           Nmax=None, skip_if_exists=True, reload_targets=False):
        if programnames is None:
            programnames = list(self.programtable.t['name'])
        print(f'updating programs {" ".join(programnames)}')
        for programname in programnames:
            self.update_json_tables_prog(programname,YYMMDD=YYMMDD,priority=priority,
                                         Nmax=Nmax,skip_if_exists=skip_if_exists,reload_targets=reload_targets)

    def update_json_tables(self, programnames=None, YYMMDD=None, priorities=None,
                           Nmax=None, skip_if_exists=True, reload_targets=False):
        
        ixs_new = self.find_new_jsonfiles(programnames=programnames, YYMMDD=YYMMDD, priorities=priorities)
        return(ixs_new)
        
    def _slew_time(self, coord1, coord2):
        """Calculate slew time between two coordinates."""
        sep = coord1.separation(coord2).to(u.deg).value
        slew_rate = 5 / 11.0
        # readout is about 29 s, and slew time graph has y-intercept of ~29 s
        # So 29 s is added either way
        slew_time = sep / slew_rate + 29.0

        # account for equatorial mount flip
        if (coord1.dec < -30.17 * u.deg) and (coord2.dec > -30.17 * u.deg):
            slew_time += 180.0
        elif (coord1.dec > -30.17 * u.deg) and (coord2.dec < -30.17 * u.deg):
            slew_time += 180.0

        return pd.Timedelta(slew_time, unit="s")
        
    def loadobsplan(self,filename=None):
        if filename is None:
            filename = self.obsplan_filename
            
        print(f'Loading obsplan {filename}')
        self.jsontable.load(filename)
        
        # the following allows to re-order the schedule outside the script in an editor: 
        # the new order is the order as read in!
        ixs_current = self.jsontable.ix_inrange('order',0,None)
        if len(ixs_current)>0:
            self.jsontable.t.loc[ixs_current,'order']=range(self.ordercounter_stepsize,(len(ixs_current)+1)*self.ordercounter_stepsize,self.ordercounter_stepsize)
        
    def final_ixs(self):
        ixs_ordered = self.jsontable.ix_sort_by_cols('order',indices=self.jsontable.ix_inrange('order',0,None))
        #ixs_notordered = AnotB(self.jsontable.getindices(),ixs_ordered)
        ixs_notordered =  self.jsontable.ix_sort_by_cols('ra',indices=AnotB(self.jsontable.getindices(),ixs_ordered))
        ixs = np.concatenate((ixs_ordered,ixs_notordered))
        
        return(ixs,ixs_ordered,ixs_notordered)
    
       
    def saveobsplan(self,filename=None, obsplan4human_filename=None):

        if filename is None:
            filename = self.obsplan_filename
        
        if obsplan4human_filename is None:
            obsplan4human_filename = self.obsplan4human_filename

        (ixs, ixs_ordered, ixs_notordered) = self.final_ixs()
        cols = self.get_param_list(self.params['obsplan_columns'])

        #if self.verbose: 
        #    self.jsontable.write(indices=ixs)
            
        print(f'Saving obsplan to {filename}')
        self.jsontable.write(filename,columns=cols, indices=ixs)
        
        print(f'Saving obsplan for humans to {obsplan4human_filename}')
        with open(obsplan4human_filename, 'w') as f:
            
            output=self.get_night_info_description()
            f.write(output+'\n\n')

            
            cols = self.get_param_list(self.params['obsplan_columns_short'])

            f.write('\n#####################\n### OBSPLAN:\n#####################\n')
            (errorflag,output) = self.jsontable.write(columns=cols, indices=ixs_ordered,return_lines=True)
            print(f'{type(output)} X1: {output}')
            f.write(output+'\n')

            f.write('\n#####################\n### EXTRA TARGETS:\n#####################\n')
            (errorflag,output) = self.jsontable.write(columns=cols, indices=ixs_notordered,return_lines=True)
            print(f'{type(output)} X2: {output}')
            f.write(output+'\n')

            f.close()
        
        return(0)
        #self.jsontable.write(obsplan4human_filename,columns=cols, indices=ixs)
        
            

    def jsonlist_radec(self):
        ixs = self.jsontable.ix_equal('priority',0)
        ixs_sorted = self.jsontable.ix_sort_by_cols('ra',indices=ixs)
        
        self.jsontable.write(columns=['program','json_short','ra','dec'],indices=ixs_sorted)

    def generic_target_ordering(self, ixs_json=None, priorities=[0]):
        ixs_generic_order_table = self.generictargetorder.getindices()
        
        if ixs_json is None:
            if len(priorities)==0:
                ixs_json = self.jsontable.getindices()
            else:
                ixs_json = []
                for prio in priorities:
                    ixs_json.extend(self.jsontable.ix_equal('priority',prio))
        
        # Don't start with 0 so that we can still add entries before the current first entry
        ordercounter = self.ordercounter_stepsize
        
        for ix in ixs_generic_order_table:
            print(f'generic entry: {self.generictargetorder.t.loc[ix,"json_short"]}')
            ixs_match = self.jsontable.ix_equal('json_short',self.generictargetorder.t.loc[ix,'json_short'],indices=ixs_json)
            ixs_match = self.jsontable.ix_matchregex('json_short',self.generictargetorder.t.loc[ix,'json_short'],indices=ixs_json)
            if len(ixs_match)>0:
                for ix_match in ixs_match:
                    self.jsontable.t.loc[ix_match,'order']=ordercounter
                    ordercounter+=self.ordercounter_stepsize
            else:
                if self.verbose>3: print(f'could not find an entry in the generic order table for {self.generictargetorder.t.loc[ix,"json_short"]}')
                continue

        ixs = self.jsontable.ix_inrange('order',0,None)
        ixs_ordered = self.jsontable.ix_sort_by_cols('order',indices=ixs)
        ixs_notordered = AnotB(ixs_json,ixs)
        cols = self.get_param_list(self.params['obsplan_columns'])
        
        if self.verbose: 
            print(f'{len(ixs_ordered)} json files found in the generic order file, {len(ixs_notordered)} not ordered')
            if self.verbose>2:
                print('### json files not ordered yet:')
                self.jsontable.write(indices=ixs_notordered)
        

            #print('### ordered json files:')
            #self.jsontable.write(columns=cols, indices=ixs_ordered)

        return(ixs_ordered,ixs_notordered)

    def calc_slew(self,ixs_ordered=None):
        if ixs_ordered is None:
            ixs_ordered = self.jsontable.ix_sort_by_cols('order',indices=self.jsontable.ix_inrange('order',0,None))
            
        ra = np.array(self.jsontable.t.loc[ixs_ordered,'ra'])
        dec = np.array(self.jsontable.t.loc[ixs_ordered,'dec'])
        
        coord = SkyCoord(ra=ra,dec=dec,unit='deg')
        
        #coord1 = SkyCoord(self.jsontable.t.loc[ixs_ordered[:-1],'ra']*u.degree,dec=self.jsontable.t.loc[ixs_ordered[:-1],'dec']*u.degree)
        #coord2 = SkyCoord(ra=self.jsontable.t.loc[ixs_ordered[1:],'ra']*u.degree,dec=self.jsontable.t.loc[ixs_ordered[1:],'dec']*u.degree)
        
        (slew_time,decflip,sep) = calc_slewtimes_1array(coord, decflip_time=self.params['decflip_seconds'])
        self.jsontable.t.loc[ixs_ordered,'t_slew[m]'] = slew_time/60.0
        self.jsontable.t.loc[ixs_ordered,'t_tot[m]'] = self.jsontable.t.loc[ixs_ordered,'t_exp[m]']  + self.jsontable.t.loc[ixs_ordered,'t_slew[m]']  
        self.jsontable.t.loc[ixs_ordered,'sep_deg'] = sep
        self.jsontable.t.loc[ixs_ordered,'decflip'] = decflip
        #self.jsontable.write(indices=ixs_ordered,columns=['json_short','ra','dec','order','priority','decflip','sep_deg','t_slew[m]','t_exp[m]','t_tot[m]'])

        

    def order_targets_DELME(self,ixs_json,ixs_ordered=None):
        if ixs_ordered is None:
            ixs_ordered = self.jsontable.ix_sort_by_cols('order',indices=self.jsontable.ix_inrange('order',0,None))
        else:
            ixs_ordered = self.jsontable.getindices(ixs_ordered)
        
            
        counter=0
        for ix_json in ixs_json:
            #if self.jsontable.t.loc[ix_json,'program']=='twi':
            #    continue
            
            ra = np.array(self.jsontable.t.loc[ixs_ordered,'ra'])
            dec = np.array(self.jsontable.t.loc[ixs_ordered,'dec'])
            coord = SkyCoord(ra=ra,dec=dec,unit='deg')

            coord_obj =SkyCoord(ra=np.full(len(ixs_ordered),self.jsontable.t.loc[ix_json,'ra']),
                                dec=np.full(len(ixs_ordered),self.jsontable.t.loc[ix_json,'dec']),unit='deg')
            (slew_time,decflip,sep) = calc_slewtimes_2arrays(coord,coord_obj, decflip_time=self.params['decflip_seconds'])
            self.jsontable.t.loc[ixs_ordered,'_slew2obj'] = slew_time/60.0
            self.jsontable.t.loc[ixs_ordered,'_sep'] = sep
            self.jsontable.t.loc[ixs_ordered,'_decflip'] = decflip
            self.jsontable.t.loc[ixs_ordered[1:],'_extra'] = np.array(self.jsontable.t.loc[ixs_ordered[:-1],'_slew2obj']) + np.array(self.jsontable.t.loc[ixs_ordered[1:],'_slew2obj']) - np.array(self.jsontable.t.loc[ixs_ordered[1:],'t_slew[m]'])
           
            
            #self.jsontable.write(indices=ixs_ordered,columns=['json_short','ra','dec','order','priority','decflip','sep_deg','t_slew[m]','_slew2obj','_extra','_sep','_decflip','t_exp[m]','t_tot[m]'])
            
            # find the index for which the extra slew time is minimal
            # the first entry of ixs_ordered is NaN, so do [1:]
            ix_minval = self.jsontable.t.loc[ixs_ordered[1:],'_extra'].idxmin()
            # now get the position of that minval in ixs_ordered
            z = np.argwhere(ixs_ordered==ix_minval)[0][0]

            if self.verbose>3:
                print(f'\n####\n#### json file to be inserted:')
                self.jsontable.write(indices=[ix_json],columns=['json_short','ra','dec','order','priority','decflip','sep_deg','t_slew[m]','_slew2obj','_extra','_sep','_decflip','t_exp[m]','t_tot[m]'])
    
                print(f'#### json file for the new entry to be inserted before: ix: {ix_minval}: {ixs_ordered}')
                self.jsontable.write(indices=[ix_minval],columns=['json_short','ra','dec','order','priority','decflip','sep_deg','t_slew[m]','_slew2obj','_extra','_sep','_decflip','t_exp[m]','t_tot[m]'])
                
                print(f'#### full table before insertion')
                self.jsontable.write(indices=ixs_ordered,columns=['json_short','ra','dec','order','priority','decflip','sep_deg','t_slew[m]','_slew2obj','_extra','_sep','_decflip','t_exp[m]','t_tot[m]'])                    
            
            
            
            extra = self.jsontable.t.loc[ix_minval,'_extra']
            # placement can be first, middle, and last
            placement='middle'
            if self.jsontable.t.loc[ix_json,'program']=='twi':
                if self.jsontable.t.loc[ixs_ordered[0],'_slew2obj'] < self.jsontable.t.loc[ixs_ordered[-1],'_slew2obj']:
                    placement = 'first'
                else:
                    placement = 'last'
            else:    
                if (self.jsontable.t.loc[ixs_ordered[0],'_slew2obj']<extra) or (self.jsontable.t.loc[ixs_ordered[-1],'_slew2obj']<extra):
                    if self.jsontable.t.loc[ixs_ordered[0],'_slew2obj'] < self.jsontable.t.loc[ixs_ordered[-1],'_slew2obj']:
                        placement = 'first'
                    else:
                        placement = 'last'
                    if (placement == 'first') and (self.jsontable.t.loc[ixs_ordered[0],'program']=='twi'):
                        placement='middle'
                    if (placement == 'last') and (self.jsontable.t.loc[ixs_ordered[-1],'program']=='twi'):
                        placement='middle'
                    
            if placement=='first':
                if self.verbose>3: print(f'EVENING {self.jsontable.t.loc[ixs_ordered[0],"_slew2obj"]} {extra} {self.jsontable.t.loc[ixs_ordered[0],"program"]}')
                self.jsontable.t.loc[ix_json,'order'] = self.jsontable.t.loc[ixs_ordered[0],'order']-1
                ixs_ordered = np.concatenate(([ix_json],ixs_ordered))
                #sys.exit(0)
                #tmp = [ix_json]
                #tmp.extend(ixs_ordered)
                #ixs_ordered = tmp
            elif placement=='last':
                if self.verbose>3: print(f'MORNING {self.jsontable.t.loc[ixs_ordered[-1],"_slew2obj"]} {extra} {self.jsontable.t.loc[ixs_ordered[-1],"program"]}')
                self.jsontable.t.loc[ix_json,'order'] = self.jsontable.t.loc[ixs_ordered[-1],'order']+self.ordercounter_stepsize
                ixs_ordered = np.concatenate((ixs_ordered,[ix_json]))
                #ixs_ordered.append(ix_json)
                #ixs_slice = ixs_ordered[-2:]
            else:
                if self.verbose>3: print(f'NORMAL {self.jsontable.t.loc[ixs_ordered[0],"_slew2obj"]} {extra}')
                self.jsontable.t.loc[ix_json,'order'] = self.jsontable.t.loc[ix_minval,'order']-1
                # ... and insert it into the correct spot into the ordered list
                ixs_ordered = np.insert(ixs_ordered,z,ix_json)
                
            
                    
                                        
            """
            if ((self.jsontable.t.loc[ixs_ordered[0],'_slew2obj']<extra) and (self.jsontable.t.loc[ixs_ordered[0],'program']!='twi')):
                print(f'EVENING {self.jsontable.t.loc[ixs_ordered[0],"_slew2obj"]} {extra} {self.jsontable.t.loc[ixs_ordered[0],"program"]}')
                self.jsontable.t.loc[ix_json,'order'] = self.jsontable.t.loc[ixs_ordered[0],'order']-1
                tmp = [ix_json]
                tmp.extend(ixs_ordered)
                ixs_ordered = tmp
            elif (self.jsontable.t.loc[ixs_ordered[-1],'_slew2obj']<extra) and (self.jsontable.t.loc[ixs_ordered[-1],'program']!='twi'):
                print(f'MORNING {self.jsontable.t.loc[ixs_ordered[-1],"_slew2obj"]} {extra} {self.jsontable.t.loc[ixs_ordered[-1],"program"]}')
                self.jsontable.t.loc[ix_json,'order'] = self.jsontable.t.loc[ixs_ordered[-1],'order']+self.ordercounter_stepsize
                ixs_ordered.append(ix_json)
                #ixs_slice = ixs_ordered[-2:]
            else:
                print(f'NORMAL {self.jsontable.t.loc[ixs_ordered[0],"_slew2obj"]} {extra}')
                self.jsontable.t.loc[ix_json,'order'] = self.jsontable.t.loc[ix_minval,'order']-1
                # ... and insert it into the correct spot into the ordered list
                ixs_ordered = np.insert(ixs_ordered,z,ix_json)
            """
                    
            
            """
            
            #print(f'{z}')
            #z = ixs_ordered.index(ix_minval)
            #print(f'{ixs_ordered[z-1:z+2]}')
            
            # check if it would be better to insert the entry at the beginning, but only if the current first entry is *not* from the twilight program, or if the new entry is from twilight
            # It is better to insert the new entry before all other entries if the slew to the first entry would be less than the _extra if it 
            # is inserted between the first and second entry
            if ((z==1) and self.jsontable.t.loc[ix_json,'program']=='twi') or ((z==1) and self.jsontable.t.loc[ixs_ordered[0],'program']!='twi' and (self.jsontable.t.loc[ixs_ordered[z-1],'_slew2obj'] < self.jsontable.t.loc[ixs_ordered[z],'_extra'])):
                self.jsontable.t.loc[ix_json,'order'] = self.jsontable.t.loc[ixs_ordered[z-1],'order']-1
                tmp = [ix_json]
                tmp.extend(ixs_ordered)
                ixs_ordered = tmp
                #ixs_slice = ixs_ordered[:2]
            # check if it would be better to insert the entry at the beginning, but only if the first entry is *not* from the twilight program
            # It is better to insert the new entry before all other entries if the slew to the first entry would be less than the _extra if it 
            # is inserted between the first and second entry
            elif ((z==len(ixs_ordered)-1) and self.jsontable.t.loc[ix_json,'program']=='twi') or (z==len(ixs_ordered)-1 and self.jsontable.t.loc[ixs_ordered[-1],'program']!='twi' and ((self.jsontable.t.loc[ixs_ordered[z],'_slew2obj'] < self.jsontable.t.loc[ixs_ordered[z-1],'_extra']))):
                self.jsontable.t.loc[ix_json,'order'] = self.jsontable.t.loc[ix_minval,'order']+self.ordercounter_stepsize
                ixs_ordered.append(ix_json)
                #ixs_slice = ixs_ordered[-2:]
            else:
                # set the 'order' of the entry to before the one with the minimum _extra
                self.jsontable.t.loc[ix_json,'order'] = self.jsontable.t.loc[ix_minval,'order']-1
                # ... and insert it into the correct spot into the ordered list
                ixs_ordered = np.insert(ixs_ordered,z,ix_json)
                
                # we need to recalculate the slew times only for the new entry, and the entry with the minimum _extra
                #print(f'ix: {ix_minval}: {ixs_ordered}')
                #ixs_slice = ixs_ordered[z-2:z+1]
                
                #print(f'{ixs_slice}')
                #self.jsontable.write(indices=ixs_slice,columns=['json_short','ra','dec','order','priority','decflip','sep_deg','t_slew[m]','_slew2obj','_extra','_sep','_decflip','t_exp[m]','t_tot[m]'])
                
                #self.calc_slew(ixs_slice)
                #self.jsontable.write(indices=ixs_ordered,columns=['json_short','ra','dec','order','priority','decflip','sep_deg','t_slew[m]','_slew2obj','_extra','_sep','_decflip','t_exp[m]','t_tot[m]'])
            """
            
            # recalculate the slew times
            #self.calc_slew(ixs_slice)
            self.calc_slew(ixs_ordered)
            if self.verbose>3:
                print(f'#### full table after insertion')
                self.jsontable.write(indices=ixs_ordered,columns=['json_short','ra','dec','order','priority','decflip','sep_deg','t_slew[m]','_slew2obj','_extra','_sep','_decflip','t_exp[m]','t_tot[m]'])                    
            # reset the order values of all the ordered entries
            #neworder=np.array(range(self.ordercounter_stepsize,(len(ixs_ordered)*self.ordercounter_stepsize)+self.ordercounter_stepsize,self.ordercounter_stepsize))
            neworder=np.array(range(len(ixs_ordered)))*self.ordercounter_stepsize+self.ordercounter_stepsize
            self.jsontable.t.loc[ixs_ordered,'order'] = neworder
            
            self.jsontable.t = self.jsontable.t.drop(columns=['_slew2obj','_extra','_sep','_decflip'])
            #if self.jsontable.t.loc[ix_json,'json_short']=='SheppardTwilightMar29morning.json':
            #    sys.exit(0)
            #if counter==200:
            #    sys.exit(0)
            counter+=1
            
        return(ixs_ordered)

    def add_unordered_targets(self,ixs_json,ixs_ordered=None):
        if ixs_ordered is None:
            ixs_ordered = self.jsontable.ix_sort_by_cols('order',indices=self.jsontable.ix_inrange('order',0,None))
        else:
            ixs_ordered = self.jsontable.getindices(ixs_ordered)
        
            
        counter=0
        for ix_json in ixs_json:
            #if self.jsontable.t.loc[ix_json,'program']=='twi':
            #    continue
            
            ra = np.array(self.jsontable.t.loc[ixs_ordered,'ra'])
            dec = np.array(self.jsontable.t.loc[ixs_ordered,'dec'])
            coord = SkyCoord(ra=ra,dec=dec,unit='deg')

            coord_obj =SkyCoord(ra=np.full(len(ixs_ordered),self.jsontable.t.loc[ix_json,'ra']),
                                dec=np.full(len(ixs_ordered),self.jsontable.t.loc[ix_json,'dec']),unit='deg')
            (slew_time,decflip,sep) = calc_slewtimes_2arrays(coord,coord_obj, decflip_time=self.params['decflip_seconds'])
            self.jsontable.t.loc[ixs_ordered,'_slew2obj'] = slew_time/60.0
            self.jsontable.t.loc[ixs_ordered,'_sep'] = sep
            self.jsontable.t.loc[ixs_ordered,'_decflip'] = decflip
            self.jsontable.t.loc[ixs_ordered[1:],'_extra'] = np.array(self.jsontable.t.loc[ixs_ordered[:-1],'_slew2obj']) + np.array(self.jsontable.t.loc[ixs_ordered[1:],'_slew2obj']) - np.array(self.jsontable.t.loc[ixs_ordered[1:],'t_slew[m]'])
           
            
            #self.jsontable.write(indices=ixs_ordered,columns=['json_short','ra','dec','order','priority','decflip','sep_deg','t_slew[m]','_slew2obj','_extra','_sep','_decflip','t_exp[m]','t_tot[m]'])
            
            # find the index for which the extra slew time is minimal
            # the first entry of ixs_ordered is NaN, so do [1:]
            ix_minval = self.jsontable.t.loc[ixs_ordered[1:],'_extra'].idxmin()
            # now get the position of that minval in ixs_ordered
            z = np.argwhere(ixs_ordered==ix_minval)[0][0]

            if self.verbose>3:
                print('\n####\n#### json file to be inserted:')
                self.jsontable.write(indices=[ix_json],columns=['json_short','ra','dec','order','priority','decflip','sep_deg','t_slew[m]','_slew2obj','_extra','_sep','_decflip','t_exp[m]','t_tot[m]'])
    
                print('#### json file for the new entry to be inserted before: ix: {ix_minval}: {ixs_ordered}')
                self.jsontable.write(indices=[ix_minval],columns=['json_short','ra','dec','order','priority','decflip','sep_deg','t_slew[m]','_slew2obj','_extra','_sep','_decflip','t_exp[m]','t_tot[m]'])
                
                print('#### full table before insertion')
                self.jsontable.write(indices=ixs_ordered,columns=['json_short','ra','dec','order','priority','decflip','sep_deg','t_slew[m]','_slew2obj','_extra','_sep','_decflip','t_exp[m]','t_tot[m]'])                    
            
            
            
            extra = self.jsontable.t.loc[ix_minval,'_extra']
            # placement can be first, middle, and last
            placement='middle'
            if self.jsontable.t.loc[ix_json,'program']=='twi':
                if self.jsontable.t.loc[ixs_ordered[0],'_slew2obj'] < self.jsontable.t.loc[ixs_ordered[-1],'_slew2obj']:
                    placement = 'first'
                else:
                    placement = 'last'
            else:    
                if (self.jsontable.t.loc[ixs_ordered[0],'_slew2obj']<extra) or (self.jsontable.t.loc[ixs_ordered[-1],'_slew2obj']<extra):
                    if self.jsontable.t.loc[ixs_ordered[0],'_slew2obj'] < self.jsontable.t.loc[ixs_ordered[-1],'_slew2obj']:
                        placement = 'first'
                    else:
                        placement = 'last'
                    if (placement == 'first') and (self.jsontable.t.loc[ixs_ordered[0],'program']=='twi'):
                        placement='middle'
                    if (placement == 'last') and (self.jsontable.t.loc[ixs_ordered[-1],'program']=='twi'):
                        placement='middle'
                    
            if placement=='first':
                print(f'EVENING {self.jsontable.t.loc[ixs_ordered[0],"_slew2obj"]} {extra} {self.jsontable.t.loc[ixs_ordered[0],"program"]}')
                self.jsontable.t.loc[ix_json,'order'] = self.jsontable.t.loc[ixs_ordered[0],'order']-1
                ixs_ordered = np.concatenate(([ix_json],ixs_ordered))
                #sys.exit(0)
                #tmp = [ix_json]
                #tmp.extend(ixs_ordered)
                #ixs_ordered = tmp
            elif placement=='last':
                print(f'MORNING {self.jsontable.t.loc[ixs_ordered[-1],"_slew2obj"]} {extra} {self.jsontable.t.loc[ixs_ordered[-1],"program"]}')
                self.jsontable.t.loc[ix_json,'order'] = self.jsontable.t.loc[ixs_ordered[-1],'order']+self.ordercounter_stepsize
                ixs_ordered = np.concatenate((ixs_ordered,[ix_json]))
                #ixs_ordered.append(ix_json)
                #ixs_slice = ixs_ordered[-2:]
            else:
                print(f'NORMAL {self.jsontable.t.loc[ixs_ordered[0],"_slew2obj"]} {extra}')
                self.jsontable.t.loc[ix_json,'order'] = self.jsontable.t.loc[ix_minval,'order']-1
                # ... and insert it into the correct spot into the ordered list
                ixs_ordered = np.insert(ixs_ordered,z,ix_json)
                
                    
            # recalculate the slew times
            self.calc_slew(ixs_ordered)
            if self.verbose>3:
                print('#### full table after insertion')
                self.jsontable.write(indices=ixs_ordered,columns=['json_short','ra','dec','order','priority','decflip','sep_deg','t_slew[m]','_slew2obj','_extra','_sep','_decflip','t_exp[m]','t_tot[m]'])                    
            # reset the order values of all the ordered entries
            neworder=np.array(range(len(ixs_ordered)))*self.ordercounter_stepsize+self.ordercounter_stepsize
            self.jsontable.t.loc[ixs_ordered,'order'] = neworder
            
            self.jsontable.t = self.jsontable.t.drop(columns=['_slew2obj','_extra','_sep','_decflip'])
            #if self.jsontable.t.loc[ix_json,'json_short']=='SheppardTwilightMar29morning.json':
            #    sys.exit(0)
            #if counter==200:
            #    sys.exit(0)
            counter+=1
            
        return(ixs_ordered)


    def calc_obsUTtime(self, ixs_ordered=None, UT_starttime=None, UT_endtime=None, firsthalf=False, secondhalf=False,
                       add2starttime_minutes=None):
        if ixs_ordered is None:
            ixs_ordered = self.jsontable.ix_sort_by_cols('order',indices=self.jsontable.ix_inrange('order',0,None))
            
        if UT_starttime is None:
            if secondhalf:
                UT_starttime = self.UT_night_midpoint
            else:
                UT_starttime = self.UT_eveningtwi.to_value('isot')
                       
        if UT_endtime is None:
            if firsthalf:
                UT_endtime = self.UT_night_midpoint
            else:
                UT_endtime = self.UT_morningtwi.to_value('isot')
                
        # Make sure the UT times are not in string format, but Time format
        if isinstance(UT_starttime,str):
            UT_starttime = Time(UT_starttime)
        if isinstance(UT_endtime,str):
            UT_endtime = Time(UT_endtime)

        if add2starttime_minutes is not None:
            UT_starttime += pd.Timedelta(add2starttime_minutes, unit="m")
             

        tznames = self.params['timezones'] 
        tz2use = self.params['timezones2use'] 
        tz = {}
        for tz_abbr in tz2use:
            tz[tz_abbr] = timezone(tznames[tz_abbr])

        for i in range(len(ixs_ordered)):
            if i==0:
                ut = UT_starttime
                if self.jsontable.t.loc[ixs_ordered[i],'program']=='twi':
                    # the twilight script *ends* at twilight
                    ut_twi = ut - pd.Timedelta(self.jsontable.t.loc[ixs_ordered[i],'t_tot[m]'], unit="m")
                    #ut_twi = self.UT_eveningtwi - pd.Timedelta(self.jsontable.t.loc[ixs_ordered[i],'t_tot[m]'], unit="m")
                    self.jsontable.t.loc[ixs_ordered[i],'UT'] = ut_twi.to_value('isot')
                    for tz_abbr in tz2use:
                        self.jsontable.t.loc[ixs_ordered[i],tz_abbr] = ut_twi.to_datetime(timezone=tz[tz_abbr]).strftime('%H:%M')
                    i+=1
            else:
                ut = Time(self.jsontable.t.loc[ixs_ordered[i-1],'UT']) + pd.Timedelta(self.jsontable.t.loc[ixs_ordered[i-1],'t_tot[m]'], unit="m")
                if i==len(ixs_ordered)-1:
                    if self.jsontable.t.loc[ixs_ordered[i],'program']=='twi':
                        self.jsontable.t.loc[ixs_ordered[i],'UT'] = UT_endtime
                        ut_last_obs_finished = ut
                        ut = UT_endtime
                    else:
                        ut_last_obs_finished = ut + pd.Timedelta(self.jsontable.t.loc[ixs_ordered[i],'t_tot[m]'], unit="m")
                    self.ExtraTime_min = TimeDelta(UT_endtime - ut_last_obs_finished).to_value('min')
                    #print('ggg',dt2.to_value('hr'),dt2.to_value('min'))
                        
            self.jsontable.t.loc[ixs_ordered[i],'UT'] = ut.to_value('isot')
            #self.jsontable.t.loc[ixs_ordered[i],'CLT'] = ut.to_datetime(timezone=ctiotz).strftime('%H:%M')
            for tz_abbr in tz2use:
                self.jsontable.t.loc[ixs_ordered[i],tz_abbr] = ut.to_datetime(timezone=tz[tz_abbr]).strftime('%H:%M')

        # make sure that the evening twilights end at -14 deg twilight
        if self.jsontable.t.loc[ixs_ordered[0],'program']=='twi':
            ut_twi = self.UT_eveningtwi - pd.Timedelta(self.jsontable.t.loc[ixs_ordered[0],'t_tot[m]'], unit="m")
            self.jsontable.t.loc[ixs_ordered[0],'UT'] = ut_twi.to_value('isot')
            for tz_abbr in tz2use:
                self.jsontable.t.loc[ixs_ordered[0],tz_abbr] = ut_twi.to_datetime(timezone=tz[tz_abbr]).strftime('%H:%M')

        # make sure that the morning twilights start at -14 deg twilight
        if self.jsontable.t.loc[ixs_ordered[-1],'program']=='twi':
            ut_twi = self.UT_morningtwi
            self.jsontable.t.loc[ixs_ordered[-1],'UT'] = ut_twi.to_value('isot')
            for tz_abbr in tz2use:
                self.jsontable.t.loc[ixs_ordered[-1],tz_abbr] = ut_twi.to_datetime(timezone=tz[tz_abbr]).strftime('%H:%M')

        
        #self.jsontable.write(indices=ixs_ordered,columns=['jsonID','json_short','ra','dec','order','priority','decflip','sep_deg','t_slew[m]','t_exp[m]','t_tot[m]','UT','CLT'])
        print(f'Extra Time left at the end of the night: {self.ExtraTime_min:.2f} minutes')
        if self.ExtraTime_min<0.0:
            print(f'\n####################################\n### WARNING!! OBSPLAN IS {-self.ExtraTime_min:.2f} minutes TOO LONG!!!\n####################################')


    def init_airmass(self):
        self.airmass_calc = AirmassCalculator(self.airmass_cache)
        self.airmass_calc.generate_airmass_grid()


    def load_airmass(self):
        self.airmass_calc = AirmassCalculator(self.airmass_cache)
        self.airmass_calc.load_airmass_grid()
        
        # timezone load
        datetime_date = datetime.strptime(self.YYMMDD, "%y%m%d")
        self.observer_shift = self.observer.timezone.utcoffset(datetime_date).total_seconds() * u.s

        ctio_observer = ctio_location()
        self.datetime = Time(datetime_date, format="datetime")
        dummy_sidereal_time = ctio_observer.local_sidereal_time(self.datetime).to(u.hourangle)
        self.sidereal_shift = dummy_sidereal_time.hour * u.hour
        self.midnight = self.observer.midnight(self.datetime, which="next")


    def calc_airmass(self):
        (ixs, ixs_ordered, ixs_notordered) = self.final_ixs()
        
        for ix in ixs_ordered:
            t = Time(self.jsontable.t.loc[ix,'UT'])
            ra = self.jsontable.t.loc[ix,'ra']*u.deg
            dec = self.jsontable.t.loc[ix,'dec']*u.deg
            airmass = self.airmass_calc.query(t,ra,dec)
            self.jsontable.t.loc[ix,'airmass']=airmass
        
        self.jsontable.t['airmass'] = self.jsontable.t['airmass'].astype(float)
        self.jsontable.default_formatters['airmass']='{:,.2f}'.format
        
    def airmass_plot(self):
        (ixs, ixs_ordered, ixs_notordered) = self.final_ixs()
        self.jsontable.t.loc[ixs_ordered,'t_start_hidden']=Time(list(self.jsontable.t.loc[ixs_ordered,'UT']))+TimeDelta(list(self.jsontable.t.loc[ixs_ordered,'t_slew[m]'].astype('float'))*u.min)
        self.jsontable.t.loc[ixs_ordered,'t_end_hidden']=Time(list(self.jsontable.t.loc[ixs_ordered,'UT']))+TimeDelta(list(self.jsontable.t.loc[ixs_ordered,'t_tot[m]'].astype('float'))*u.min)

        fig = self.airmass_calc.plot_airmass(
            self.jsontable.t.loc[ixs_ordered], 
            self.UT_eveningtwi, #self.start_time, 
            self.UT_morningtwi, #self.end_time, 
            [self.sidereal_shift, self.observer_shift]
        )
        print(f'Saving airmass plot to {self.airmass_filename}')
        fig.savefig(self.airmass_filename, bbox_inches="tight", dpi=300)
        #plt.show()
        plt.close()

    def show_plots(self):
        cmd = f'open {self.airmass_filename} {self.moon_filename}'
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        cmd_in,cmd_out = p.stdin,p.stdout
        cmd_out.readlines()
                
        return(0)

    def calc_moon_distances(self):
        moon_calc = MoonCalculator()
        (ixs, ixs_ordered, ixs_notordered) = self.final_ixs()
        
        for ix in ixs_ordered:
            t = Time(self.jsontable.t.loc[ix,'UT'])
            ra = self.jsontable.t.loc[ix,'ra']*u.deg
            dec = self.jsontable.t.loc[ix,'dec']*u.deg
            sep_moon = moon_calc.get_moon_distance(ra,dec,t)
            self.jsontable.t.loc[ix,'sep_moon']=sep_moon
        
        self.jsontable.t['sep_moon'] = self.jsontable.t['sep_moon'].astype(float)
        self.jsontable.default_formatters['sep_moon']='{:,.1f}'.format

    def moon_plot(self):
        print(f'MOON PLOT START...')
        moon_calc = MoonCalculator()
        (ixs, ixs_ordered, ixs_notordered) = self.final_ixs()
        self.jsontable.t.loc[ixs_ordered,'t_start_hidden']=Time(list(self.jsontable.t.loc[ixs_ordered,'UT']))+TimeDelta(list(self.jsontable.t.loc[ixs_ordered,'t_slew[m]'].astype('float'))*u.min)
        self.jsontable.t.loc[ixs_ordered,'t_end_hidden']=Time(list(self.jsontable.t.loc[ixs_ordered,'UT']))+TimeDelta(list(self.jsontable.t.loc[ixs_ordered,'t_tot[m]'].astype('float'))*u.min)
        
        #self.jsontable.write(columns=['UT','t_slew[m]','t_tot[m]','t_start_hidden','t_end_hidden'],indices=ixs_ordered)

        fig = moon_calc.plot_moon_distances(
            self.jsontable.t.loc[ixs_ordered], 
            self.UT_eveningtwi, #self.start_time, 
            self.UT_morningtwi, #self.end_time, 
            [self.sidereal_shift, self.observer_shift]
        )
        print(f'Saving moon distance plot to {self.moon_filename}')
        fig.savefig(self.moon_filename, bbox_inches="tight", dpi=300)
        #plt.show()
        plt.close()
        print(f'MOON PLOT END...')


if __name__ == '__main__':

    test = obsplan_baseclass()
    parser = test.define_optional_arguments()
    args = parser.parse_args()

    # the arguments are saved in query.params
    test.get_arguments(args)
    test.initialize4date(test.params['date'])
    
    a = TimeDelta([100.0,200.0]*u.min)
    print('aaa',a.to_value('min'))
    b = Time(['2025-04-01T02:28:31.510','2025-04-01T02:30:31.510'])
    print('b',b)
    c = b+a
    print('c',c)
    
    sys.exit(0)
    
    test.load_airmass()
    ra = np.array([100.0,200.0])
    dec = np.array([-30.0,30.0])
    times = Time(['2025-04-01T01:42:11.154','2025-04-01T02:16:04.813'])
    
    ra=100.0*u.deg
    dec=-10.0*u.deg
    #times=Time('2025-04-01T01:42:11.154')
    
    coord = SkyCoord(ra,dec,unit='deg')
    print('fff',coord)
    airmass = test.airmass_calc.query(times,ra,dec)
    print('fff111',airmass)
    
    
    ras = np.array([100.0,200.0]*u.deg)
    decs = np.array([-30.0,30.0]*u.deg)
    times = Time(['2025-04-01T01:42:11.154','2025-04-01T02:16:04.813'])
    airmass = test.airmass_calc._query_direct(times,ras,decs)
    print('fff222',airmass)
    sys.exit(0)
    
    test.update_json_tables()
    test.update_tables()
    
    test.generic_target_ordering()
    
    test.calc_obsUTtime()
    
    test.wrap_up_and_summary()
    test.saveobsplan()
    #test.jsonlist_radec()
    sys.exit(0)
    test.set_maindirs()
    print(f'\n# Outrootdir: {test.outrootdir()}')
