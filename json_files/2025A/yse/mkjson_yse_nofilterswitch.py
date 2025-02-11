#!/usr/bin/env python
import re,sys,string,math,os,types,time,fcntl,shutil,random
# put the tools directory into the path
if 'PIPE_PYTHONSCRIPTS' in os.environ:
    sys.path.append(os.environ['PIPE_PYTHONSCRIPTS']+'/tools')
import optparse
from texttable import txttableclass,sex2deg,unique,makepath4file
from astropy.time import Time
from astroplan import Observer,moon
from astropy import units as u
from pytz import timezone
from astropy.coordinates import (EarthLocation, SkyCoord, AltAz, get_sun,
                                 get_moon, Angle, Longitude)
ctio = Observer.at_site("CTIO",timezone='America/Santiago')
#ctio = Observer.at_site("CTIO",timezone='UTC')



def rmfile(filename,raiseError=1,gzip=False):
    " if file exists, remove it "
    if os.path.lexists(filename):
        os.remove(filename)
        if os.path.isfile(filename):
            if raiseError == 1:
                raise RuntimeError('ERROR: Cannot remove %s' % filename)
            else:
                return(1)
    if gzip and os.path.lexists(filename+'.gz'):
        os.remove(filename+'.gz')
        if os.path.isfile(filename+'.gz'):
            if raiseError == 1:
                raise RuntimeError('ERROR: Cannot remove %s' % filename+'.gz')
            else:
                return(2)
    return(0)

 
class mkjsonclass(txttableclass):
    def __init__(self):
        txttableclass.__init__(self)
        self.jsonscripttemplate=['"count": 1', 
                                 '"expType": "object"',  
                                 '"object": "XXX_OBJECT_XXX"', 
                                 '"filter": "XXX_FILTER_XXX"', 
                                 '"expTime": XXX_EXPTIME_XXX']

        self.setdither1()
        
        # readout time in seconds, used to estimate the time for finishing the script
        self.readouttime_sec = 30

        self.pointingcollist = ['Name','ID','pointing']
        self.racollist = ['ra','Ra','RA']
        self.deccollist = ['dec','Dec','DEC']
        self.racol = 'ra'
        self.deccol = 'dec'
        self.pointingcol = 'pointing'
        
        self.horizons = [16,15,14,12,10]


    def setrectangledither(self,delta=80):
        self.dither = {}
        self.dither[0]=['"RA": "XXX_RA_XXX"','"dec": "XXX_DEC_XXX"']
        self.dither[1]=['"deltaDEC": "S%d.000000"' % (delta),'"deltaRA": "E%d.000000"' % (delta)]
        self.dither[2]=['"deltaDEC": "N%d.000000"' % (2*delta),'"deltaRA": "E0.000000"']
        self.dither[3]=['"deltaDEC": "N0.000000"','"deltaRA": "W%d.000000"' % (2*delta)]
        self.dither[4]=['"deltaDEC": "S%d.000000"' % (2*delta),'"deltaRA": "W0.000000"']

        self.options.Ndithers=5

        self.rapattern=re.compile('XXX_RA_XXX')
        self.decpattern=re.compile('XXX_DEC_XXX')
        self.pointingpattern=re.compile('XXX_OBJECT_XXX')
        self.filterpattern=re.compile('XXX_FILTER_XXX')
        self.exptimepattern=re.compile('XXX_EXPTIME_XXX')


    def setdither1(self):
        self.dither = {}
        self.dither[0]=['"RA": "XXX_RA_XXX"','"dec": "XXX_DEC_XXX"']
        self.dither[1]=['"deltaDEC": "S60.000000"','"deltaRA": "E60.000000"']
        self.dither[2]=['"deltaDEC": "S20.000000"','"deltaRA": "W40.000000"']
        self.dither[3]=['"deltaDEC": "N100.000000"','"deltaRA": "W80.000000"']
        self.dither[4]=['"deltaDEC": "N60.000000"','"deltaRA": "W20.000000"']
        #
        self.dither[5]=['"deltaDEC": "S20.000000"','"deltaRA": "E20.000000"']
        self.dither[6]=['"deltaDEC": "S20.000000"','"deltaRA": "W40.000000"']
        self.dither[7]=['"deltaDEC": "S20.000000"','"deltaRA": "W20.000000"']
        self.dither[8]=['"deltaDEC": "N30.000000"','"deltaRA": "E50.000000"']
        self.dither[9]=['"deltaDEC": "N20.000000"','"deltaRA": "E20.000000"']

        self.rapattern=re.compile('XXX_RA_XXX')
        self.decpattern=re.compile('XXX_DEC_XXX')
        self.pointingpattern=re.compile('XXX_OBJECT_XXX')
        self.filterpattern=re.compile('XXX_FILTER_XXX')
        self.exptimepattern=re.compile('XXX_EXPTIME_XXX')

    def updatedither(self,idither,dither):
        if dither != None:
            if not re.search('^[E|W]\d+\.*\d*$',dither[0]):
                print('%s is not in the form EX.XXX or WX.XXX, where X.XXX is the shift in arcsec' % dither[0])
                sys.exit(0)
            if not re.search('^[S|N]\d+\.*\d*$',dither[1]):
                print('%s is not in the form SX.XXX or NX.XXX, where X.XXX is the shift in arcsec' % dither[1])
                sys.exit(0)

            self.dither[idither]=['"deltaDEC": "%s"' % (dither[1]),'"deltaRA": "%s"' % (dither[0])]
            

    def loadfile(self,filename):
        if not os.path.isfile(filename): raise RuntimeError('file %s does nto exist!' % filename)
        try:
            txttableclass.loadfile(self,filename)
        except:
            data = open(filename,'r').readlines()
            print('No header in file %s? usind first line: %s' % (filename,data[0]))
            self.parsetable(data[0],data[1:])
           
        if not(self.racol in self.cols):
            foundflag=False
            for ra in self.racollist:
                if (ra in self.cols):
                    self.racol=ra
                    foundflag=True
            if not foundflag: raise RuntimeError('Could not find RA column from %s in header columns %s' % (','.join(self.racollist),','.join(self.cols)))
 
        if not(self.deccol in self.cols):
            foundflag=False
            for dec in self.deccollist:
                if (dec in self.cols):
                    self.deccol=dec
                    foundflag=True
            if not foundflag: raise RuntimeError('Could not find Dec column from %s in header columns %s' % (','.join(self.deccollist),','.join(self.cols)))
         
        if not(self.pointingcol in self.cols):
            foundflag=False
            for pointing in self.pointingcollist:
                if (pointing in self.cols):
                    self.pointingcol=pointing
                    foundflag=True
            if not foundflag: raise RuntimeError('Could not find pointing column from %s in header columns %s' % (','.join(self.pointingcollist),','.join(self.cols)))
    
        self.configcols([self.pointingcol],'s','%s',visible=1)
        for key in self.allrowkeys:
            self.setentry(key,self.racol,sex2deg(self.getentry(key,self.racol),ra=True))
            self.setentry(key,self.deccol,sex2deg(self.getentry(key,self.deccol),ra=False))
        self.configcols([self.racol,self.deccol],'f','%.4f',visible=1)
        self.configcols(['Nstars'],'d','%d',visible=1)

    def add_options(self, parser=None, usage=None):
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")
        parser.add_option('--pointingfile'  , default='./YSE.pointings.txt' , type="string",
                          help='RA column name (default=%default)')
        parser.add_option('--racol'  , default='ra' , type="string",
                          help='RA column name (default=%default)')
        parser.add_option('--deccol'  , default='dec' , type="string",
                          help='Dec column name (default=%default)')
        parser.add_option('--pointingcol'  , default='pointing' , type="string",
                          help='pointing column name (default=%default)')
        parser.add_option('-p','--priority'  , default=0 , action="count",
                          help='priority (default=%default)')
        parser.add_option('-o','--outrootdir'  , default="." , type="string",
                          help='base directory (default=%default)')
        parser.add_option('-s','--outsuffix'  , default=None , type="string",
                          help='add outsuffix before the .json. If "filter", then the filters are added (default=%default)')
        parser.add_option('-d','--outdate'  , default='now' , type="string",
                          help='outdate is used as the subdir for the output, needs to be in YYDDMM format. If "now", then the current date in YYMMDD is used (default=%default)')
        parser.add_option('--outsubdir'  , default=None , type="string",
                          help='outsubdir can be used to override outdate and use not a date string as subdir (default=%default)')
        parser.add_option('-e','--exptime'  , default='50' , type="string",
                          help='base exposure times of single dither. If comma-separated list, then corresponds to filters (default=%default)')
        parser.add_option('-f','--filter'  , default='r' , type="string",
                          help='filters, comma-separated list (default=%default)')
        parser.add_option('--propid'  , default='2024B-763968' , type="string",
                          help='proposal ID (default=%default)')
        parser.add_option('-n','--Ndithers'  , default=1 , type="int",
                          help='base number of dithers (default=%default)')
        parser.add_option('--N1'  , default=10000 , type="int",
                          help='Nstars value for which exptime gets halfed (default=%default)')
        parser.add_option('--N2'  , default=30000 , type="int",
                          help='Nstars value for which exptime gets 1/3 (default=%default)')
#        parser.add_option('--N1'  , default=15000 , type="int",
#                          help='Nstars value for which exptime gets halfed (default=%default)')
#        parser.add_option('--N2'  , default=40000 , type="int",
#                          help='Nstars value for which exptime gets 1/3 (default=%default)')
        parser.add_option('--pointings'  , default=None , type="string",
                          help='comma-separated list of pointings (default=%default)')
        parser.add_option('--sortbyra'  , default=False , action="store_true",
                          help='sort the output table by RA')
        
        parser.add_option('--d1'  , default=None , nargs=2, type="string",
                          help='1st dither in RA and DEC direction (default=%default)')
        parser.add_option('-r','--rectangledither',default=None , type="int",
                          help='center+rectangle dither pathern. 80" is good to close chip gaps (default=%default)')
        parser.add_option('-b','--bigdither',default=False, action="store_true",
                          help='center+rectangle big dither pathern with delta=1200" (default=%default)')
        parser.add_option('--d2'  , default=None , nargs=2, type="string",
                          help='2nd dither in RA and DEC direction (default=%default)')
        parser.add_option('--d3'  , default=None , nargs=2, type="string",
                          help='3rd dither in RA and DEC direction (default=%default)')
        parser.add_option('--d4'  , default=None , nargs=2, type="string",
                          help='4th dither in RA and DEC direction (default=%default)')
        parser.add_option('--d5'  , default=None , nargs=2, type="string",
                          help='5th dither in RA and DEC direction (default=%default)')
        parser.add_option('--d6'  , default=None , nargs=2, type="string",
                          help='6th dither in RA and DEC direction (default=%default)')
        parser.add_option('--d7'  , default=None , nargs=2, type="string",
                          help='7th dither in RA and DEC direction (default=%default)')
        parser.add_option('--d8'  , default=None , nargs=2, type="string",
                          help='8th dither in RA and DEC direction (default=%default)')
        parser.add_option('--d9'  , default=None , nargs=2, type="string",
                          help='9th dither in RA and DEC direction (default=%default)')
        return(parser)

    def replace_placeholders(self,cmds,key,exptime=10.0,camerafilter='r'):
        rastring = '%.4f' % self.getentry(key,self.racol)
        decstring = '%.4f' % self.getentry(key,self.deccol)
        pointing = self.getentry(key,self.pointingcol)
        exptimestring = '%d' % exptime
        for i in range(len(cmds)):
            cmd = self.rapattern.sub(rastring,cmds[i])
            cmd = self.decpattern.sub(decstring,cmd)
            cmd = self.pointingpattern.sub(pointing,cmd)
            cmd = self.filterpattern.sub(camerafilter,cmd)
            cmds[i] = self.exptimepattern.sub(exptimestring,cmd)
        return(cmds)

    def mkjsonscript4exposure(self,key,dithercounter=0,exptime=10.0,camerafilter='r',changepositionflag=True):
        cmds=[]
        cmds.extend(self.jsonscripttemplate)
        if changepositionflag:
            cmds.extend(self.dither[dithercounter])
        else:
            if dithercounter==0:
                cmds.extend(self.dither[dithercounter])
            else:
                cmds.extend(['"deltaDEC": "S0.000000"','"deltaRA": "E0.000000"'])            
        cmds.append('"propid": "%s"' % self.options.propid)
        cmds = self.replace_placeholders(cmds,key,exptime=exptime,camerafilter=camerafilter)
        jsonscript = '{\n' + ',\n'.join(cmds)+ '\n}'
        return(jsonscript)

    def mkjsonscript4object(self,key,filters=['r'],exptimes=[100],repeatfilters=[]):
        exposures4object = []

        #exptime = self.options.exptime
        #Ndithers = self.options.Ndithers

        ditherfactor = 1
        if self.options.N2!=None and self.getentry(key,'Nstars') and self.getentry(key,'Nstars')>self.options.N2:
            #exptime= int(exptime/3)
            ditherfactor = 3 
        elif self.options.N1!=None and self.getentry(key,'Nstars') and self.getentry(key,'Nstars')>self.options.N1:
            #exptime= int(exptime/2)
            ditherfactor = 2
        
        Ndithers = self.options.Ndithers*ditherfactor
            

        t = 0.0
        for n in range(Ndithers):
            for i in range(len(filters)):
                if len(exptimes)>1:
                    exptime = exptimes[i]
                else:
                    exptime = exptimes[0]
                exptime=int(exptime/ditherfactor)
                camerafilter=filters[i] 
                print('pointing %s, #%2d, exptime=%4d filter:%s' % (self.getentry(key,self.pointingcol),n+1,exptime,camerafilter))
                exposures4object.append(self.mkjsonscript4exposure(key,dithercounter=n,exptime=exptime,camerafilter=camerafilter,changepositionflag=(i==0)))
                t += (exptime+self.readouttime_sec)
                if camerafilter in repeatfilters:
                    print('pointing %s, #%2d, exptime=%4d filter:%s' % (self.getentry(key,self.pointingcol),n+1,exptime,camerafilter))
                    exposures4object.append(self.mkjsonscript4exposure(key,dithercounter=n,exptime=exptime,camerafilter=camerafilter,changepositionflag=False))
                    t += (exptime+self.readouttime_sec)
                    
        return(exposures4object,t)

    def mkjsonscript4all(self,targets,combineFlag=False,outrootdir=None,outsuffix=None,outsubdir=None,priority=0,repeatfilters=[],sortbyra=False):

        def outfile(basename,outrootdir=None,outsuffix=None,outsubdir=None,priority=0):
            outfilename = outrootdir
            if outfilename is None:
               outfilename ='.'
            if outsubdir is not None:
               outfilename += f'/{outsubdir}'
            if priority>0:
               outfilename += f'/{priority}'
               
            outfilename +=f'/{basename}' 
            if outsuffix is not None:
               outfilename += f'.{outsuffix}'
            outfilename += '.json'
            outfilename = os.path.abspath(outfilename)
            makepath4file(outfilename)
            return(outfilename)
        
        
        
        keysused = []

        targets2keyhash = {}
        keys=self.rowkeys()
        for key in keys:
            if int(self.getentry(key,'skip'))==1:
                continue
            group = self.getentry(key,'group')
            if group in targets:
                if group not in targets2keyhash:
                    targets2keyhash[group]=[]
                targets2keyhash[group].append(key)

            pointing = self.getentry(key,self.pointingcol)
            if pointing in targets:
                if pointing not in targets2keyhash:
                    targets2keyhash[pointing]=[]
                targets2keyhash[pointing].append(key)
                
        filters = self.options.filter.split(',')
        exptimes = [int(e) for e in self.options.exptime.split(',')]

        if len(exptimes)!=1 and len(exptimes)!=len(filters):
            raise RuntimeError("ERROR: # of filters (%s) need to much # of exposure times (%s)!" % (self.options.filter,self.options.exptime))

        self.configcols(['Moon_sep'],'f','%.1f',visible=True)
        self.configcols(['Moon_illum'],'f','%.2f',visible=True)

        for target in targets2keyhash:
            print(f'\n#############################\n### target {target}')
            exposures = []
            keys = targets2keyhash[target]
            keysused.extend(keys)
            tall = 0.0

            for filt in filters:
                print(f'### filter {filt}')
                for key in keys:
                    #target =  SkyCoord(34*u.deg,-4*u.deg)
                    sep = self.moon.separation(SkyCoord(self.getentry(key,self.racol)*u.deg,self.getentry(key,self.deccol)*u.deg))
                    print(f'pointing {self.getentry(key,self.pointingcol)} RA/Dec=({self.getentry(key,self.racol)*u.deg:.4f},{self.getentry(key,self.deccol)*u.deg:.4f}), Moon separation {sep.deg:.2f}')
                    self.setentry(key,'Moon_sep',sep.deg)
                    self.setentry(key,'Moon_illum',self.illum)
                    
                    (exposure,t) = self.mkjsonscript4object(key,[filt],exptimes,repeatfilters=repeatfilters)
                    tall += t
                    exposures.extend(exposure)
                #keys.reverse()
              
            outfilename = outfile(target,outrootdir=outrootdir,outsuffix=outsuffix,outsubdir=outsubdir,priority=priority)
            print('total time for script (exposure time and %f seconds readout):\n %.0f seconds, %.1f minutes, %.2f hours' % (self.readouttime_sec,tall,tall/60.0,tall/3600.0))       
            rmfile(outfilename)
            print('Saving ', outfilename)
            open(outfilename,'w').writelines('[\n' + ',\n'.join(exposures)+ '\n]\n')
        
        if sortbyra: keysused = self.sortkeysbycols(keysused,'RA')
        self.printtxttable(keys=keysused)
        print('!!! setting propID: %s !!!' % self.options.propid)
        
        
        
if __name__=='__main__':
    
    if 1==0:
        
        print(ctio,ctio.location)
        
        YYMMDD = '211028'
        date = f'20{YYMMDD[:2]}-{YYMMDD[2:4]}-{YYMMDD[4:6]}T06:00:00'
        t = Time(date)
        print(t)
        t+=1
        print('vvv',t)
        m = moon.moon_phase_angle(t)
        illum = moon.moon_illumination(t)
        print(m.to_value('deg'),illum)
        
        moon = get_moon(t,ctio.location)
        print('moon',moon)
        target =  SkyCoord(34*u.deg,-4*u.deg)
        print('target',target)
        sep = moon.separation(target)
        print(sep.deg)
        sys.exit(0)
    
    mkjson = mkjsonclass()

    parser = mkjson.add_options(usage = 'mkjson.py pointingfile')
    mkjson.options,  args = parser.parse_args()

    if parser.has_option('help') and mkjson.options.help:
        sys.exit(0)


    mkjson.racol = mkjson.options.racol
    mkjson.deccol = mkjson.options.deccol
    mkjson.pointingcol = mkjson.options.pointingcol

    outsubdir = None
    outdate = mkjson.options.outdate
    # if outdate is now, get the current date
    if mkjson.options.outdate.lower()=='now':
        t = Time.now()
        m = re.search('^\d\d(\d\d)\-(\d\d)\-(\d\d)T',t.to_value('isot'))
        if m is not None:
            YYMMDD = ''.join(m.groups())
            outdate = YYMMDD
        else:
            raise RuntimeError('Could not parse {} for YYMMDD'.format(t.to_value('isot')))

    #print(outdate)
    #sys.exit(0)
    
    # make sure outdate is formatted correctly
    if re.search('^\d{6}$',outdate) is not None:
        #all good!
        pass
    elif re.search('^\d\d(\d\d)\-(\d\d)\-(\d\d)',outdate) is not None:
        m = re.search('^\d\d(\d\d)\-(\d\d)\-(\d\d)',outdate)
        if m is not None:
            YYMMDD = ''.join(m.groups())
            outdate = YYMMDD
        else:
            raise RuntimeError('Could not parse {} for YYYY-MM-DD'.format(outdate))
    else:
        if outdate!='':
            raise RuntimeError('Could not parse {} for YYMMDD'.format(outdate))
    
    # just a sanity test for years, months, and days!
    if outdate!='':
        if not(outdate[:2] in ['21','22','23','24','25']):
            raise RuntimeError(f'year {outdate[:2]} is not in years 21-24!')
        if not(int(outdate[2:4])<=12):
            raise RuntimeError(f'month {outdate[2:4]} is not <=12!')
        if not(int(outdate[4:6])<=31):
            raise RuntimeError(f'day {outdate[2:4]} is not <=31!')
        outsubdir = outdate
        
        # Take the current North/South America date, add 1 day, and 6am, 
        # to get the appropriate UTC date form the middle of the night or slightly later
        date = f'20{outdate[:2]}-{outdate[2:4]}-{outdate[4:6]}T06:00:00'
        t = Time(date)+1
        print(f'Setting Moon to UT date {t}')
        mkjson.moon = get_moon(t,ctio.location)
        mkjson.illum = moon.moon_illumination(t)
        print(f'MOON ILLUMINATION: {mkjson.illum:.2f}')
        mkjson.twi={}
        for horizon in mkjson.horizons:
            mkjson.twi[horizon] =  ctio.tonight(t-0.5,horizon=-horizon*u.deg)
            print('UT %d deg twilight: %s %s' % (-horizon,mkjson.twi[horizon][0].to_value('isot'),mkjson.twi[horizon][1].to_value('isot')))
        
        ctiotz = timezone('America/Santiago')
        for horizon in mkjson.horizons:
            dt0 = mkjson.twi[horizon][0].to_datetime(timezone=ctiotz)
            dt1 = mkjson.twi[horizon][1].to_datetime(timezone=ctiotz)
            print('LOCAL %d deg twilight: %s   %s' % (-horizon,dt0.strftime("%m/%d/%Y,   %H:%M:%S"),dt1.strftime("%H:%M:%S")))
            #print(dt.strftime("%m/%d/%Y, %H:%M:%S"))
            #loc_dt = eastern.localize()
    else:
        print('WARNING: no date specified, cannot calculate the Moon parameters')
 
    if mkjson.options.outsubdir is not None:
        outsubdir = mkjson.options.outsubdir
        
    print('outdate:',outdate)
    print('outsubdir:',outsubdir)
    
    mkjson.updatedither(1,mkjson.options.d1)
    mkjson.updatedither(2,mkjson.options.d2)
    mkjson.updatedither(3,mkjson.options.d3)
    mkjson.updatedither(4,mkjson.options.d4)
    mkjson.updatedither(5,mkjson.options.d5)
    mkjson.updatedither(6,mkjson.options.d6)
    mkjson.updatedither(7,mkjson.options.d7)
    mkjson.updatedither(8,mkjson.options.d8)
    mkjson.updatedither(9,mkjson.options.d9)

    if mkjson.options.rectangledither!=None:
        mkjson.setrectangledither(delta=mkjson.options.rectangledither)

    pointingfile = mkjson.options.pointingfile
    
    repeatfilters=[]

    mkjson.loadfile(pointingfile)
    mkjson.mkjsonscript4all(args,
                            outrootdir=mkjson.options.outrootdir,
                            outsuffix=mkjson.options.outsuffix,
                            outsubdir=outsubdir,
                            priority=mkjson.options.priority,
                            repeatfilters=repeatfilters,
                            sortbyra=mkjson.options.sortbyra)
    
    if outdate!='':
        print(f'MOON ILLUMINATION: {mkjson.illum:.2f}')

        for horizon in mkjson.horizons:
            print('UT %d deg twilight: %s %s' % (-horizon,mkjson.twi[horizon][0].to_value('isot'),mkjson.twi[horizon][1].to_value('isot')))
            #if horizon==14:
            #    dt_14_h=(mkjson.twi[horizon][1]-mkjson.twi[horizon][0]).to_value('hour')
            #    print(f'time between {-horizon} deg twilight: {dt_14_h:.2f}h')
            #    sys.exit(0)

        ctiotz = timezone('America/Santiago')
        for horizon in mkjson.horizons:
            dt0 = mkjson.twi[horizon][0].to_datetime(timezone=ctiotz)
            dt1 = mkjson.twi[horizon][1].to_datetime(timezone=ctiotz)
            print('LOCAL %d deg twilight: %s   %s' % (-horizon,dt0.strftime("%m/%d/%Y,   %H:%M:%S"),dt1.strftime("%H:%M:%S")))
        night_midpoint = mkjson.twi[mkjson.horizons[0]][0]+0.5*(mkjson.twi[mkjson.horizons[0]][1]-mkjson.twi[mkjson.horizons[0]][0])
        dt = night_midpoint.to_datetime(timezone=ctiotz)
        print('LOCAL night MIDPOINT:',dt.strftime("%m/%d/%Y,   %H:%M:%S"))

