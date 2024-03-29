#!/usr/bin/env python
#from __future__ import unicode_literals

import re,sys,string,math,os,types,time
if 'PIPE_PYTHONSCRIPTS' in os.environ:
    sys.path.append(os.environ['PIPE_PYTHONSCRIPTS']+'/tools')
#import tools
import numpy


### Converts sexigesimal notation to decimal degrees or decimal hours (if option 'ra=True' invoked)
def sex2deg(sexigecimal, ra=False):
    ### 2005/12/02 - AR: make sure it is in sexagesimal format!
    # is it a string? if not check if it is None
    if not (type(sexigecimal) is str): #types.StringType):
        if type(sexigecimal) == None:
            raise RuntimeError("ERROR: sex2deg cannot handle 'None'!")
        return sexigecimal
    # Does it have ':' or ' '? If not, it must be a float in string format, just convert it and return
    if re.search('\:|\s',sexigecimal) == None:
        return(float(sexigecimal))

    s1, s2, s3 = list(map(float, re.split('[DHMSdhms:\s]', sexigecimal.strip())[0:3]))
    # Get the sign
    if re.search('-', sexigecimal):
        sign = -1
    else:
        sign = 1

    deg = abs(s1) + s2 / 60. + s3 / 3600.

    deg *= sign

    if ra:
        deg *= 15.

    return deg

def hex2int(val):
    if type(val) is str:#types.StringType:
        val = int(eval(val))
    return(val)

# some little defs to compare lists and tuples
def unique(seq):
    if seq == None or seq == []: return []
    d = {}
    for x in seq:
        d[x] = 1
    return list(d.keys())

def makepath(path,raiseError=1):
    if path == '':
        return(0)
    if not os.path.isdir(path):
        os.makedirs(path)
        if not os.path.isdir(path):
            if raiseError == 1:
                raise RuntimeError('ERROR: Cannot create directory %s' % path)
            else:
                return(1)
    return(0)

def makepath4file(filename,raiseError=1):
    path = os.path.dirname(filename)
    if not os.path.isdir(path):
        return(makepath(path,raiseError=raiseError))
    else:
        return(0)

# class to test for skipping
# most important class parameters:
#   self.usemask, self.usemaskOP, self.skipmask, self.skipmaskOP
#
# most important routine: test4use(self,value,maskargument='default'):
# if maskargument is not None, then no mask is applied
# if maskargument=number, then (skipmask,skipmaskOP,usemask,usemaskOP)=(number,'or',None,None)
# if maskargument=tuple, then (skipmask,skipmaskOP,usemask,usemaskOP)=(tuple), where the
# tuple is filled up with None if necessary
#
# algorithm:
# test4use(value,maskarguments):
# as default: use it (usemask and skipmask == None)
# if usemask!=None: use only if (applymask(value,usemask,usemaskOP)
# return 1 if (applymask(value,usemask,usemaskOP) and not applymask(value,skipmask,skipmaskOP)
# where applymask(value,mask,maskOP):
# if mask!=None:
#    V = value & mask
#    return 1 if maskOP='or' and V>0 (Use if at least one bit is set in both, value and mask)
#    return 1 if maskOP='and' and V==value (Use only if all bits in value are also set in mask)
#    return 1 if maskOP='atleast' and V==mask (Use only if all bits in mask are also set in value)
#    return 1 if maskOP='exact' and value==mask
#    else return 0
# else return 1
#
#
# Examples:
# * usemask=0x5,         usemaskOP='or':      any value that has bit 0x1 and/or 0x4 set is used, e.g.,
#                                             0x1,0x4,0x5,0x9 etc (0x1+X, 0x4+X, X any number)
# * usemask=0x5,         usemaskOP='and':     ONLY values=0x1,0x4,0x5 are used
# * usemask=0x5,         usemaskOP='atleast': ONLY values=0x5+X are used, where X can be any number
# * usemask=0x5,         usemaskOP='exact':   ONLY value=0x5 is used
#
# * usemask=0x5,         usemaskOP='or',
# * skipmask=0x9,        skipmaskOP='exact':  any value that has bit 0x1 and/or 0x4 set is used, e.g.,
#                                             0x1,0x4,0x5,0x7 etc (0x1+X, 0x4+X, X any number).
#                                             Exception: 0x9
# Some special cases:
# * skpimask=0xffffffff, usemaskOP='or':       only value==0 is not used
# * usemask=0x0,         usemaskOP='atleast':  only value==0 is used
#
class maskclass:
    def __init__(self,maskargument='default'):
        (self.usemask,self.usemaskOP,self.skipmask,self.skipmaskOP)=(None,None,None,None)
        self.definemask(maskargument)
    def applymask(self,value,mask,maskOP):
        if mask!=None:
            if value == None:
                return(0)

            if (maskOP=='exact'):
                if (value==mask):return(1)

                else:
                    if (type(value) is str) or (type(mask) is str):
                        return(0)
            V = value & mask
            if (maskOP=='or'):
                if (V>0):return(1)
            elif (maskOP=='and'):
                if (V==value):return(1)
            elif (maskOP=='atleast'):
                if (V==mask):return(1)
            elif (maskOP=='exact'):
                if (value==mask):return(1)
            else:
                raise RuntimeError('ERROR: wrong maskOP:').with_traceback(maskOP)
        else:
            return(1)


    def test4use(self,value,maskargument='default'):
        #get the masks
        (usemask,usemaskOP,skipmask,skipmaskOP)=self.parsemaskargument(maskargument)
        # make a quick test
        if (usemask==None and skipmask==None):
            return(1)

        #default
        useflag=1
        skipflag=0
        if (usemask!=None):
            useflag=self.applymask(value,usemask,usemaskOP)
        if (skipmask!=None):
            skipflag=self.applymask(value,skipmask,skipmaskOP)
        return(useflag and (not skipflag))


    def definemask(self,maskargument):
        (self.usemask,self.usemaskOP,self.skipmask,self.skipmaskOP)=self.parsemaskargument(maskargument)

    def parsemaskargument(self,maskargument):

        # use default values if wanted
        if maskargument=='default':
            return(self.usemask,self.usemaskOP,self.skipmask,self.skipmaskOP)

        # don't do any masking!
        if maskargument==None:
            return(None,None,None,None)

        # parse the argument
        if type(maskargument) is int:
            return(maskargument,'or',None,None)
        elif type(maskargument) in [tuple,list]:
            if len(maskargument)==1:
                # standard usemask, operator is 'or'
                return(maskargument[0],'or',None,None)
            elif len(maskargument)==2:
                return(maskargument[0],maskargument[1],None,None)
            elif len(maskargument)==3:
                return(maskargument[0],maskargument[1],maskargument[2],'or')
            elif len(maskargument)==4:
                return(maskargument[0],maskargument[1],maskargument[2],maskargument[3])
            else:
                raise RuntimeError('bad arguments for parsemaskargument:').with_traceback(maskargument)

class txttableclass:
    def __init__(self, **keys):
        # used in combination with self.data[key][self.col4mask] to select data
        self.mask=maskclass()
        # default: col4mask='_mask'. This defines which col is used when masking is applied
        self.col4mask='_mask'

        # clear
        self.clear()
        self.clearcoldefs()

        self.verbose=0

        # Debugging setting
        self.debug = False

    def cleartable(self):
        self.internalrowid=0
        self.data={}
        self.allrowkeys=[]
        self.filename=''
        self.header=[]
        self.comments=[]
        #self.cols=[]
        #self.clearcoldefs

    def clearcoldefs(self):
        self.cols=[]
        self.colinfo={}

        self.header4parsing=None

        self.colmap={}
        self.cols2map=[]

    def clearproperties(self):
        self.properties={}


    def cleardefs(self):
        self.errorflag=0

        # define how to find and treat the header
        self.indicator4headerpattern='^#'
        # get rid of header indicator, e.g. leading '#'
        self.deleteindicator4header=1
        self.deleteleadingwhitespace=1
        self.deletetrailingwhitespace=1

        # define how the line are read in
        self.inputseparator='\s+'
        # define what pattern indicates undefined or missing data (only for numerical columns!)
        self.setpattern4undefined('^-$')
        # what are comments?
        self.indicator4comments='^#'
        # skip empty lines?
        self.skipemptylines=1
        # require that there is an entry, i.e. the number of entries in the line has to match the number of
        # entries in the header
        self.requireentry=1
        # if there are tabs in an entry, expand them to whitespace (Note, that happens AFTER parsing!)
        self.expandtabs=0

        # define how the table should be printed
        self.outputseparator=' '
        self.outputundefined='-'
        self.outputadd2header='#'
        self.outputadd2line=' '
        self.outputadd2lineend=''

        # THIS SHOULD BE '_id'! if it is different, then it will use the values (string format) of the specified
        # column as rowkeys. This is somewhat dangerous, e.g. keys are not guaranteed to be unique.
        # CHANGE THIS PARAMETERS AT YOUR OWN RISK!
        self.col4rowkey='_id'

        # align strings to the left
        self.stringleftalign = False

        self.phantomfillchar = '0'

    def clear(self):
        self.cleartable()
        self.clearcoldefs()
        self.cleardefs()
        self.clearproperties()

    def setcolmap(self,col2map,col):
        self.colmap[col2map]=col
        self.cols2map=list(self.colmap.keys())

    def setcol2visible(self,cols,visible):
        if type(cols) is str:
            cols=[cols,]
        for col in cols:
            if not visible==None:
                self.colinfo[col]['visible']=visible
            else:
                if 'visible' not in self.colinfo[col]:
                    self.colinfo[col]['visible']=0

    def setcol2latexphantom(self,cols,latexphantomflag=True):
        if type(cols) is str:
            cols=[cols,]
        for col in cols:
            if col in self.colinfo:
                if not latexphantomflag==None:
                    self.colinfo[col]['latexphantomflag']=latexphantomflag
                else:
                    if 'latexphantomflag' not in self.colinfo[col]:
                        self.colinfo[col]['latexphantomflag']=False

    def setpattern4undefined(self,pattern):
        if pattern==None:
            self.pattern4undefined=None
        else:
            self.pattern4undefined=re.compile(pattern)

    def entrydefined(self,key,col):
        if self.getentry(key,col)==None:
            return(0)
        else:
            return(1)

    def getcol2use(self,col):
        if col in self.cols2map:
            return self.colmap[col]
        else:
            return col

    def getentry(self,key,col, verbose=False, exit4Error=True):
        """
        Return the value from the given row (key) and column.
        Formats known column values to specified type.
        See 'txttableclass.configcols'.
        Only accepts single values of 'key' and 'col'.

        If 'key' not in 'self.data' returns None.
        """
        col=self.getcol2use(col)
        #print 'NNNNN:',col,self.data[key].has_key(col),self.data[key].keys()
        #if not self.data.has_key(key):  return None
        try:
#            key, col, self.data[key][col]
#            print self.data

            #if not self.data[key].has_key(col):
            #    return None
            if self.data[key][col] is None:
                return(None)
            if self.colinfo[col]['type']=='o':
                return(self.data[key][col])
            if self.colinfo[col]['type'] in ['d','x']:
                return(int(self.data[key][col]))
            if self.colinfo[col]['type']=='f':
                return(float(self.data[key][col]))
            if self.colinfo[col]['type'] in ['ld','lx']:
                return(int(self.data[key][col]))
            #if self.colinfo[col]['type']=='s':
                # AR 05/17/2018: we need to add unicode.
            if sys.version_info >= (3, 0):
                if type(self.data[key][col]) in [str,bytes]:#,unicode]:
                    if type(self.data[key][col]) in [str,bytes]:
                        return(self.data[key][col])
                    if type(self.data[key][col])==float:
                        return("%f" % self.data[key][col])
                    if type(self.data[key][col])==int:
                        return("%d" % self.data[key][col])
                    if type(self.data[key][col])==int:
                        return("%ld" % self.data[key][col])
                    if type(self.data[key][col]) in [int,float,str]:
                        return(int(self.data[key][col]))
                    raise RuntimeError('ERROR: Don\'t know to convert '+self.data[key][col]+' of type '+type(self.data[key][col])+' into a string!')
                else:
                    raise RuntimeError('ERROR: wrong type %s for column %s' % (self.colinfo[col]['type'],col))
            else:
                if type(self.data[key][col]) in [str,bytes]:
                    if type(self.data[key][col]) in [str,bytes]:
                        return(self.data[key][col])
                    if type(self.data[key][col])==float:
                        return("%f" % self.data[key][col])
                    if type(self.data[key][col])==int:
                        return("%d" % self.data[key][col])
                    if type(self.data[key][col])==int:
                        return("%ld" % self.data[key][col])
                    if type(self.data[key][col]) in [int,float,str]:
                        return(int(self.data[key][col]))
                    raise RuntimeError('ERROR: Don\'t know to convert '+self.data[key][col]+' of type '+type(self.data[key][col])+' into a string!')
                else:
                    raise RuntimeError('ERROR: wrong type %s for column %s' % (self.colinfo[col]['type'],col))

        except Exception as e:
            if verbose:
                print('Error (texttable.getentry): returning unformatted value')
                print(e)
                print(('key:',key,'col:',col, 'val:', self.data[key][col], 'datatype:', type(self.data[key][col]), 'coltype', self.colinfo[col]['type']))
                if not exit4Error:
                    return self.data[key][col]

    def inc(self, key, col, incval=1, verbose=True, exit4Error=True ):
        v = self.getentry(key,col,verbose=verbose, exit4Error=exit4Error)
        if v==None:
            if exit4Error:
                print(("Cannot increment None value for key ",key," and col ",col))
                raise RuntimeError('None value not allowed')
            else:
                v=0
        v+=incval
        self.setentry(key,col,v)

    def setentry(self,key,col,value, verbose=True):
        col=self.getcol2use(col)
        try:
            if value!=None:
                if self.colinfo[col]['type']=='o':
                    self.data[key][col]=value
                elif (type(value) is str) and (self.pattern4undefined!=None) and self.pattern4undefined.search(value):
                    self.data[key][col]=None
                elif self.colinfo[col]['type']=='f':
                    self.data[key][col]=float(value)
                elif self.colinfo[col]['type']=='d':
                    self.data[key][col]=int(float(value))
                elif self.colinfo[col]['type']=='x':
                    if type(value) is str:
                        self.data[key][col]=int(value,16)
                    else:
                        self.data[key][col]=int(value)
                elif self.colinfo[col]['type']=='ld':
                    self.data[key][col]=int(value)
                elif self.colinfo[col]['type']=='lx':
                    if type(value) is str:
                        self.data[key][col]=int(value,16)
                    else:
                        self.data[key][col]=int(value)
                elif self.colinfo[col]['type']=='s':
                    self.data[key][col]=value
                elif self.colinfo[col]['type']=='deg':
                    # convert to decimal degrees
                    self.data[key][col]=sex2deg(value)
                else:
                    raise RuntimeError('ERROR: wrong type %s for column %s' % (self.colinfo[col]['type'],col))
            else:
                self.data[key][col]=None
        except ValueError as e:
            if verbose:
                print("setentry error")
                print(e)
                print(('key:',key,'col:',col))
                print(('value:',self.data[key][col]))


    def setentries(self,key,cols,values):
        if type(cols) is str:
            cols=[cols,]
        if type(values) is str:
            values=[values,]
        if len(cols) != len(values):
            raise RuntimeError('ERROR: need some number of cols (%d) than values (%d)' % (len(cols),len(values)))
        irange = list(range(len(cols)))
        for i in irange:
            self.setentry(key,cols[i],values[i])

    def getentries(self,keys,col,unique=0):
        if keys == None:
            return None
        if (unique):
            return(unique([(self.getentry(key,col)) for key in keys]))
        else:
            return([(self.getentry(key,col)) for key in keys])

    def setcol2value(self,col,value,mask='default',col4mask='default',keys=None):
        if keys== None:
            keys = self.rowkeys(mask=mask,col4mask=col4mask)
        for key in keys:
            self.setentry(key,col,value)

    def setcols2value(self,cols,value,mask='default',col4mask='default',keys=None):
        if keys == None:
            keys = self.rowkeys(mask=mask,col4mask=col4mask)
        for key in keys:
            for col in cols:
                self.setentry(key,col,value)

    def clearcol(self, col, keys=None):
        if keys == None:
            keys = self.allrowkeys
        for key in keys:
            self.data[key][col] = None

    def redotypecasting(self,col):
        for key in self.allrowkeys:
            self.setentry(key,self.getcol2use(col),self.data[key][col])

    def search4entry(self,col,value,mask='default',col4mask='default'):
        keys = self.rowkeys(mask=mask,col4mask=col4mask)
        for key in keys:
            if self.getentry(key,col)==value:
                return(key)
        return None

    def search4entries(self,col,value,mask='default',col4mask='default'):
        keys = self.rowkeys(mask=mask,col4mask=col4mask)
        matchkeys = []
        for key in keys:
            if self.getentry(key,col)==value:
                matchkeys.append(key)
        return matchkeys

    def setcoltype(self,col,coltype=None,colformat=None):
        # coltype indicates the type of the column
        # (f=float,x=hexadec.,s=string,d=integer). This is used when
        # the table is read to typecast the values!
        # colformat is used to format the col when written out

        # don't do anything if coltype=None
        if coltype==None:
            if col not in self.colinfo:
                # error if col doesn't exist
                print(('ERROR: coltype must be defined for non-existing column: ',col))
                exit(0)
        else:
            self.colinfo[col]['type']=coltype
            # check for error
            if not coltype in ('f','s','d','x','o','ld','lx','deg'):
                raise RuntimeError('ERROR: unknown column type: ').with_traceback(coltype)

        # if there is no colformat info: set the default colformat
        if colformat==None and self.colinfo[col]['format']==None:
            colformat = 'default'

        # now actually do set the format accordingly:
        if colformat == 'default':
            #default format
            if self.colinfo[col]['type'] in ('f','s','d'):
                self.colinfo[col]['format']='%'+self.colinfo[col]['type']
            elif self.colinfo[col]['type']=='x':
                self.colinfo[col]['format']='0x%x'
            elif self.colinfo[col]['type']=='o':
                self.colinfo[col]['format']='%s'
            elif self.colinfo[col]['type']=='deg':
                self.colinfo[col]['format']='%f'
            else:
                print(('ERROR: unknown column type: ',coltype))
                exit(0)
        # Don't change format if colformat == None, otherwise assign new format
        elif not colformat == None:
            # set the format to the format given
            self.colinfo[col]['format']=colformat

    def configcols(self,cols,coltype=None,colformat=None,visible=None,latexphantomflag=False,defaultvalue=None):
        """
        Set the type of the value returned for a given column.
        E.g. from 'getentry(key,col)'.
        """

        if type(cols) is str:
            cols=[cols,]
        for col in cols:            # test if the column already exist
            if col == None: continue
            if col in self.colinfo:
                newcolflag=0
                oldcoltype=self.colinfo[col]['type'] # save the previous setting
                # if the coltype is changed, and no new format given: give it default format
                if (not (oldcoltype == coltype)) and colformat==None and (not (coltype==None)):
                    colformat='default'
            else:
                newcolflag=1
                self.cols.append(col)
                self.colinfo[col]={} #initialize colinfo
                oldcoltype=''
                # as default: columns are type string
                if coltype==None:
                    coltype='s'
                # new col: give it default format if none given
                if colformat==None:
                    colformat='default'
            # set the type self.colinfo[col]['type'] and self.colinfo[col]['format']
            self.setcoltype(col,coltype,colformat=colformat)

            # set if the column is visible, i.e. if it is printed by default by printtexttable
            self.setcol2visible(col,visible)

            # latex table: set if instead of spaces you want to use phantom{0}
            self.setcol2latexphantom(col,latexphantomflag)

            # set column to the defaultvalue if necessary
            if newcolflag or (not defaultvalue==None):
                self.setcol2value(col,defaultvalue)
                self.colinfo[col]['autoformat']='%s'
            else:
                # redo typecasting if necessary
                if (not newcolflag) and (not coltype==oldcoltype) :
                    self.redotypecasting(col)

    def colexist(self,col):
        return(col in self.cols)

    def colsexist(self,cols,verbose=0):
        for col in cols:
            if not self.colexist(col):
                if verbose: print(('ERROR: col %s does not exist!' % col))
                return(0)
        return(1)

    def checktxttable4error(self):
        # some error checking...
        #if len(self.allrowkeys)==0:
        #    print 'ERROR: no measurements found in '+self.filename
        #    return(1)
        if not self.colsexist(self.cols):
            print(('ERROR: something is wrong with the col definition in ' +self.filename))
            return(2)

        return(0)

    def getrow_as_hash(self,key,cols=None,raiseError=True):
        if cols==None:
            cols = self.__cols2use__()
        if key in self.allrowkeys:
            h = {}
            for col in cols:
                h[col]=self.getentry(key,col)
            return(h)
        else:
            if raiseError:
                raise RuntimeError('ERROR: key %s out of range!' % (key))
            else:
                return None

    def getrow(self,key):
        pass
        if key in self.allrowkeys:
            return(self.data[key])
        else:
            raise RuntimeError('ERROR: key %s out of range!' % (key))

    def add2row(self,rowkey,dict2add,forceExit=False):
        for col in dict2add:
            if self.colexist(col):
                self.setentry(rowkey,col,dict2add[col])
            else:
                print(('ERROR (texttable.txttableclass.add2row): You are trying to add a column, %s, that does not exist!' % col))
                if forceExit:
                    raise RuntimeError("Exiting...")
                else:
                    print('I am going to keep going without adding the entry for this column, I hope you know what you are doing.\n')

    def newrow(self,dict2add,newcolsokflag=0,key=None):
        if not (type(dict2add) is dict):
            print('ERROR: wrong type passed to setrow, dictionary expected!')
            return(1)

        if key == None:
            # get the rowkey, in general  self.internalrowid
            if self.col4rowkey!='_id':
                if self.col4rowkey not in dict2add:
                    raise RuntimeError('ERROR: trying to add new row to table, but col key %s doesnt exist!' % self.col4rowkey)
                rowkey = dict2add[self.col4rowkey]
            else:
                rowkey=self.internalrowid
        else:
            rowkey = key

        self.data[rowkey]={}

        # add new cols if necessary
        if newcolsokflag:
            for col in list(dict2add.keys()):
                if not self.colexist(col):
                    self.configcols([col],'s')


        # make sure all cols exist
        for col in self.cols:
            self.data[rowkey][col]=None

        # set some special columns
        #self.data[rowkey]['_line']=None
        self.data[rowkey]['_id']=self.internalrowid
        self.data[rowkey]['_mask']=0

        # set the rest of the columns
        for col in dict2add:
            if self.colexist(col):
                self.setentry(rowkey,col,dict2add[col])

        self.internalrowid+=1

        # get the primary keys
        self.allrowkeys.append(rowkey)
        #self.allrowkeys=self.data.keys()

        return rowkey

    def newrow_quick(self,dict2add, cols, key = None):
        """ This is a quicker newrow, but with less error checking: The cols are passed and have to exist! """
        if not (type(dict2add) is dict):
            print('ERROR: wrong type passed to setrow, dictionary expected!')
            return(1)

        if key == None:
            # get the rowkey, in general  self.internalrowid
            if self.col4rowkey!='_id':
                if self.col4rowkey not in dict2add:
                    raise RuntimeError('ERROR: trying to add new row to table, but col key %s doesnt exist!' % self.col4rowkey)
                rowkey = dict2add[self.col4rowkey]
            else:
                rowkey=self.internalrowid
        else:
            rowkey = key

        self.data[rowkey]={}

        # make sure all cols exist
        for col in self.cols:
            self.data[rowkey][col]=None

        # set some special columns
        #self.data[rowkey]['_line']=None
        self.data[rowkey]['_id']=self.internalrowid
        self.data[rowkey]['_mask']=0

        # set the rest of the columns
        for col in cols:
            self.setentry(rowkey,col,dict2add[col])

        self.internalrowid+=1

        # get the primary keys
        self.allrowkeys.append(rowkey)
        #self.allrowkeys=self.data.keys()

        return rowkey

    def getcol(self,col,keys=None,mask='default',col4mask='default'):
        if keys == None: keys = self.rowkeys(mask=mask,col4mask=col4mask)
#        return [self.data.getentry(c, col) for c in keys]
        return [self.data[c][col] for c in keys]

    def delrow(self,rowkey):
        del(self.data[rowkey])
        self.allrowkeys=list(self.data.keys())

    def parsetable(self,headerstring,lines,startindex=0,stopindex=None,addflag=0,skipheaderline=1,skipcomments=1,onlyparse2nextheader=0,maxsplit=0, onlyparsecols=None,noheader=False,comments=None,saveoriginalline=False):
        if not addflag:
            self.cleartable()

        if comments != None:
            self.comments.extend(comments)

        #This is a very risky thing to do and generally it's not recommended but sometimes you receive
        #data tables that have no header and it would be useful to be able to parse them.
        if noheader:
            headerstring = lines[startindex]
            colindex = 1

        # first parse the header to get the position of the primary cols
        self.header4parsing = headerstring  # save the original header

        if self.deleteleadingwhitespace:
            headerstring = headerstring.lstrip()
        if self.deletetrailingwhitespace:
            headerstring = headerstring.rstrip()

        # header=headerstring.split(self.inputseparator)
        headerinfo=re.split(self.inputseparator,headerstring,maxsplit=maxsplit)

        # if the col is not set yet, set it to a string
        headerinfoindices = list(range(len(headerinfo)))
        if onlyparsecols:
            newheaderinfoindices = []
            for col in onlyparsecols:
                for c in headerinfoindices:
                    if headerinfo[c] == col:
                        newheaderinfoindices.append(c)
                        break
            headerinfoindices = newheaderinfoindices
        for c in headerinfoindices:
            if noheader:
                col = 'col_'+str(colindex)
                headerinfo[c] =col
                if self.expandtabs:
                    headerinfo[c]=headerinfo[c].expandtabs()
                self.configcols(col, visible=1)
                colindex+=1
            else:
                if self.expandtabs:
                    headerinfo[c]=headerinfo[c].expandtabs()
                if not headerinfo[c] in self.cols:
                    self.configcols(headerinfo[c],visible=1)

        # special case: should rarely be used!
        col4rowkeyindex=None
        if self.col4rowkey!='_id':
            for c in headerinfoindices:
                if self.col4rowkey==headerinfo[c]:
                    col4rowkeyindex=c
                    break
            if col4rowkeyindex==None:
                raise RuntimeError('Could not find column %s!' % (self.col4rowkey))

        # add internal cols if not yet added
        if not '_id' in self.cols:
            self.configcols('_id','d','%d',visible=0)
        if saveoriginalline and (not '_line' in self.cols):
            self.configcols('_line','s','%s',visible=0)
        if (not '_mask' in self.cols):
            self.configcols('_mask','x','0x%04x',visible=0)

        # save the starting id
        internalrowidstart=self.internalrowid

        # precompile the search pattern: speeds it up!
        headersearchpattern=re.compile(self.indicator4headerpattern)
        commentsearchpattern=re.compile(self.indicator4comments)

        ##### now parse the table
        # set the startindex
        imin=startindex
        if (len(lines)<=0):
            #print "empty table"
            return(imin)

        # set the last index
        if stopindex!=None:
            imax=stopindex
        else:
            imax=len(lines)

        # We have to set this to something because we're returning it at the end
        i = imin
        irange = list(range(imin,imax))
        for i in irange:
            if skipheaderline or onlyparse2nextheader:
                # only do the parsing test after the first real line is read
                headerflag=headersearchpattern.search(lines[i])
                #if self.debug:
                #    print 'TEST111111111111111111111111111111',i,headerflag,onlyparse2nextheader,self.internalrowid,internalrowidstart,lines[i]
                # stop if next header is there
                if onlyparse2nextheader and self.internalrowid>internalrowidstart and headerflag:
                    break
                # skip header
                if skipheaderline and headerflag:
                    i+=1
                    continue


            if skipcomments and commentsearchpattern.search(lines[i]):
                self.comments.append(lines[i])
                i+=1
                continue

            if self.deleteleadingwhitespace:
                lines[i] = lines[i].lstrip()
            if self.deletetrailingwhitespace:
                lines[i] = lines[i].rstrip()
            if self.skipemptylines and lines[i]=='':
                i+=1
                continue

            lineinfo= re.split(self.inputseparator,lines[i],maxsplit=maxsplit)
            #lineinfo= re.split(self.inputseparator,lines[i])
            if len(lineinfo) != len(headerinfo):
                if (self.requireentry):
                    print((lineinfo,maxsplit))
#                    print (headerinfo[8], lineinfo[8])
#                    print (headerinfo[18], lineinfo[18])
#                    print (headerinfo[38], lineinfo[38])
#                    print lineinfo[38]
                    print((len(headerinfo), len(lineinfo)))
                    raise RuntimeError('While parsing table: header and line have unequal number of elements: %d != %d !!!! (%s != %s)' % (len(headerinfo), len(lineinfo),headerinfo,lineinfo))

            if col4rowkeyindex == None:
                rowkey=self.internalrowid
            else:
                rowkey=lineinfo[col4rowkeyindex]


            self.data[rowkey]={}
            # safe the original line if wanted
            if saveoriginalline:
                self.data[rowkey]['_line']=lines[i]
            # count the line number (excludes comments and header)
            self.data[rowkey]['_id']=self.internalrowid
            # set mask to 0
            self.data[rowkey]['_mask']=0

            for c in headerinfoindices:
                if c<len(lineinfo):
                    if self.expandtabs:
                        lineinfo[c]=lineinfo[c].expandtabs()
                    #print 'BBBBBBBBBBBBBBBBBBBBB',c,headerinfo[c],lineinfo[c]
                    self.setentry(rowkey,headerinfo[c],lineinfo[c])
                else:
                    self.setentry(rowkey,headerinfo[c],None)

            self.internalrowid+=1
            i+=1

        # get the primary keys
        self.allrowkeys=list(self.data.keys())
        # if _mask does not exist, add it to the cols
        if (not '_mask' in self.cols):
            self.configcols('_mask','x','0x%04x',defaultvalue=0)

        return(i)

    def search4nextheader(self,lines,startindex=0,stopindex=None,takelastheader=0,header1stline=False):
        headerstring=""
        leftoverlines=[]
        headersearchpattern=re.compile(self.indicator4headerpattern)
        commentsearchpattern=re.compile(self.indicator4comments)

        # set the last index
        if stopindex!=None:
            imax=stopindex+1
        else:
            imax=len(lines)

        i=startindex
        irange = list(range(startindex,imax))
        for i in irange:
        #for i in range(startindex,imax):
            # header?
            if headersearchpattern.search(lines[i]):
                self.header.append(lines[i])
                if self.deleteindicator4header:
                    headerstring=headersearchpattern.sub('',lines[i],1)
                else:
                    headerstring=lines[i]
                # only break if the next line is not a header as well
                if takelastheader:
                    if (i+1<imax) and (not headersearchpattern.search(lines[i+1])):
                        break
                    else:
                        i+=1
                        #continue
                else:
                    break
            # comment?
            if i>=imax:
                continue
            if commentsearchpattern.search(lines[i]):
                self.comments.append(lines[i])
                i+=1
            else:
                if header1stline:
                    i=0
                    headerstring=lines[0]
                else:
                    raise RuntimeError("Unaccountable lines when trying to find next header: %s" % lines[0])
        # some more error checking
        if headerstring=='':
            raise RuntimeError('could not find the header!')
        return headerstring,i

    def getfromdb(self,db,columns,command,keys=None,columnconversion=None):

        if columnconversion!=None:
            convcols = list(columnconversion.keys())
            for col in columns:
                if not (col in convcols):
                   columnconversion[col]=col
        #print columnconversion

        newkeys = []
        c = db.cursor()
        db.text_factory=str
        #db.row_factory=dict
        c.execute(command)
        i = 0
        for row in c:
            #print 'xxx',row
            if keys !=None:
                key = keys[i]
            else:
                key = None
            key = self.parsedboutputtuple(columns,row,key=key,columnconversion=columnconversion)
            newkeys.append(key)
            i+=1
        return(newkeys)

    def parsedboutputtuple(self,columns,dboutput,key=None,columnconversion=None):
        if key == None:
            key = self.newrow({})
        if len(columns) != len(dboutput):
            print(('ERROR: # of columns (%d) unequal # of db output (%d)' % (len(columns),len(dboutput))))
            sys.exit(0)
        for i in range(len(columns)):
            if columnconversion!=None:
                col = columnconversion[columns[i]]
                self.setentry(key,col,dboutput[i])
            else:
                self.setentry(key,columns[i],dboutput[i])
        return(key)

    def loadfile(self,filename,addflag=0,maxsplit=0,takelastheader=0, onlyparsecols=None,noheader=False,headerstring=None, loadstartline=None,headerindex=None,header1stline=False):
        self.errorflag=0

        if not os.path.isfile(filename):
            print(('Warning: The file %s does not exist!' % (filename)))
            self.errorflag=1
            return(1)

        # read the files
        lines=open(filename).readlines()

        # get the header
        if noheader:
            headerstring = None
            i = 0
        elif headerstring!=None:
            i = 0
        elif headerindex !=None:
            (headerstring,i) = self.search4nextheader(lines,startindex=headerindex,takelastheader=takelastheader,header1stline=header1stline)
        else:
            (headerstring,i)=self.search4nextheader(lines,startindex=0,takelastheader=takelastheader,header1stline=header1stline)

        if loadstartline:
            i = loadstartline
        # parse the table and fill self.data

        self.parsetable(headerstring,lines,startindex=i,addflag=addflag,maxsplit=maxsplit, onlyparsecols=onlyparsecols,noheader=noheader,comments=self.comments)
        self.filename=filename
        self.errorflag=self.checktxttable4error()
        return(self.errorflag)

    def loadfitsfile(self,fitsfilename,addflag=0,headerstring=None,defaultheaderstring=True,loadstartline=1,getzptmag=True):
        if defaultheaderstring:
            if headerstring:
                print(('WARNING: headerstring specified and defaultheaderstring=True.\n I think you want to use the specified headerstring, so that is what I am going to do:\n %s' % headerstring))
            else:
                dheaderstring = 'Xpos\tYpos\tM\tdM\tflux\tdflux\ttype\tpeakflux\tsigx\tsigxy\tsigy\tsky\tchisqrt\tclass\tFWHM1\tFWHM2\tFWHM\tangle\textendedness\tflag\tmask\tNmask'
                print(("Using default headerstring for fits file:\n%s" % dheaderstring))
                headerstring = dheaderstring
        errorflag = self.loadfile(fitsfilename, headerstring=headerstring,loadstartline=loadstartline)
        if getzptmag:
            import subprocess
            zptmag = subprocess.getoutput('gethead %s ZPTMAG' % fitsfilename)
            if len(zptmag.split()) !=1 :
                print(('Something went wrong with getting ZPTMAG: %s\nSetting it to None' % zptmag))
                self.zptmag = None
            self.zptmag = float(zptmag)
        return errorflag

    def configcmpcols(self,cols2config=['Xpos','Ypos','M','flux','dflux','type']):
        for col in cols2config:
            if   col == 'Xpos': self.configcols('Xpos', 'f', '%.3f', visible=1)
            elif col == 'Ypos': self.configcols('Ypos', 'f', '%.3f', visible=1)
            elif col == 'M': self.configcols('M', 'f', '%.4f', visible=1)
            elif col == 'flux': self.configcols('flux', 'f', '%.5f', visible=1)
            elif col == 'dflux': self.configcols('dflux', 'f', '%.6f', visible=1)
            elif col == 'type': self.configcols('type', 'x', visible=1)
            else:
                print(('Unknown col for configcmpcols %s\nYou need to configure this col on your own!' % col))

    def parsefitsheaderline(self,fitsheaderline,keypattern,stringvalpattern,valpattern,commentpattern):
        # get the keyword
        r=keypattern.match(fitsheaderline)
        if (r==None) or (len(r.groups())!=2):
            return('','','','')
        (key,restofline)=r.groups()

        # get the value ...
        val=''
        comment=''
        type=''
        r=stringvalpattern.search(restofline)
        # first, test if it is a string
        if r:
            (val,comment)=r.groups()
            type='string'
        else:
            # if not it's a number
            r=valpattern.search(restofline)
            if r:
                (val,comment)=r.groups()
                type='numeric'

        # strip unnecessary stuff
        val     = val.strip()
        comment = commentpattern.sub('',comment).strip()
        #print 'key:%s  val:<%s>  type:<%s>   comment:<%s>' % (key,val,type,comment)

        return(key,val,type,comment)

    def parsefitsheader(self,hdrAsString):
        hdr = {}
        # how many lines are in the header?
        N=int((len(hdrAsString)-1)/80.0)
        if N!=(len(hdrAsString)-1)/80.0:
            raise RuntimeError('ERROR: Something is wrong with this fitsheader, number of chars (%d) is not a multiple of 80 ' % (len(hdrAsString)-1))

        # compile the patterns. makes it faster
        keypattern = re.compile('(\S+)\s*=\s*(.*)')
        stringvalpattern = re.compile('\'(.*)\'(.*)')
        valpattern = re.compile('(\S*)(.*)')
        commentpattern = re.compile('\s*/')

        #loop through the header
        for c in range(N):
            hdrline = hdrAsString[c*80:(c+1)*80]
            #print hdrline
            (fitskey,val,type,comment)=self.parsefitsheaderline(hdrline,keypattern,stringvalpattern,valpattern,commentpattern)
            if fitskey != '':
                hdr[fitskey]={}
                hdr[fitskey]['val']=val
                hdr[fitskey]['type']=type
                hdr[fitskey]['comment']=comment
                hdr[fitskey]['#']=c
                hdr[fitskey]['line']=hdrline
        return(hdr)

    def loadcmpfile(self,filename,addflag=0,maxsplit=0, onlyparsecols=None):
        self.errorflag=0
        self.filename=filename
        if not os.path.isfile(filename):
            self.cmpheader = {}
            print(('Warning: The file %s does not exist!' % (filename)))
            self.errorflag=1
            return(1)

        import time
        # read the files
        lines=open(self.filename).readlines()
        # read the fitsheader
        self.cmpheader = self.parsefitsheader(lines[0])

        if 'NCOLTBL' not in self.cmpheader:
            print(('Warning: The file %s seems not to be a cmpfile! NCOLTBL is missing from the header' % (filename)))
            self.errorflag=1
            return(1)

        # create the headerstring
        headerstring = ''
        for c in range(1,int(self.cmpheader['NCOLTBL']['val'])+1):
            headerstring = headerstring + " " + self.cmpheader['COLTBL%d' % c]['val']

        # parse the table and fill self.data
        self.parsetable(headerstring,lines,startindex=1,addflag=addflag,maxsplit=maxsplit, onlyparsecols=onlyparsecols)

        self.errorflag=self.checktxttable4error()
        self.filename=filename
        return(self.errorflag)

    def loadcmpheader(self,filename,addflag=0,maxsplit=0):
        self.errorflag=0
        self.filename=filename

        if not os.path.isfile(filename):
            print(('Warning: The file %s does not exist!' % (filename)))
            self.errorflag=1
            return(1)

        # read the files
        lines=open(self.filename).readlines()
        # read the fitsheader
        self.cmpheader = self.parsefitsheader(lines[0])

        if 'NCOLTBL' not in self.cmpheader:
            print(('Warning: The file %s seems not to be a cmpfile! NCOLTBL is missing from the header' % (filename)))
            self.errorflag=1
            return(1)

#         # create the headerstring
#         headerstring = ''
#         for c in range(1,int(self.cmpheader['NCOLTBL']['val'])+1):
#             headerstring = headerstring + " " + self.cmpheader['COLTBL%d' % c]['val']

#         # parse the table and fill self.data
#         self.parsetable(headerstring,lines,startindex=1,addflag=addflag,maxsplit=maxsplit)

        self.errorflag=0#self.checktxttable4error()
        return(self.errorflag)


    def loadsexcatfile(self,filename,addflag=0):
        self.errorflag=0



    def save2file(self,filename,**args):
        if 'verbose' in args:
            print(('Saving ',filename))
        f = open(filename,'w')
        self.printtxttable(file=f,**args)
        f.close()
        self.savedfilename = filename


    def append2file(self,filename,**args):
        f = open(filename,'a')
        self.printtxttable(file=f,**args)
        f.close()
        self.savedfilename = filename

    def save2texfile(self,filename,cols=None,colnamehash={},colalignhash={},
                     onlydata2file=False,twocolumnflag=False,tablelabelstring='tab:',**args):
        outputseparator_bkp = self.outputseparator
        outputadd2lineend_bkp = self.outputadd2lineend
        self.outputseparator=' & '
        self.outputadd2lineend='\\\\'
        #self.outputundefined=outputundefined
        colline=''
        colnumberline=''
        formatstring = ''
        cols = self.__cols2use__(cols=cols)
        for i in range(len(cols)):
            if colline != '':
                colline += ' & '
                colnumberline += ' & '
            if cols[i] in colnamehash:
                colname=colnamehash[cols[i]]
            else:
                colname=cols[i]
            if cols[i] in colalignhash:
               formatstring +=colalignhash[cols[i]]
            else:
                if self.colinfo[cols[i]]['type']=='s':
                    formatstring += 'l'
                else:
                    formatstring += 'r'

            colline += '\colhead{%s}' % colname
            colnumberline += '\colhead{(%d)}' % i
        f = open(filename,'w')
        if onlydata2file:
            fheader = open(filename+'.header','w')
            ftail = open(filename+'.tail','w')
        else:
            fheader=f
            ftail=f
        if twocolumnflag:
            fheader.write('\\begin{deluxetable*}{%s}\n' % formatstring)
        else:
            fheader.write('\\begin{deluxetable}{%s}\n' % formatstring)
        fheader.write('\\tabletypesize{\\scriptsize}\n')
        fheader.write('\\tablecaption{\n')
        fheader.write('\\label{%s}}\n' % tablelabelstring)
        fheader.write('\\tablehead{\n%s\\\\\n%s\n}\n\startdata\n' % (colline,colnumberline))
        self.printtxttable(file=f,printheader=0,cols=cols,**args)
        ftail.write('\\enddata\n')
        ftail.write('\\tablecomments{\n')
        if twocolumnflag:
            ftail.write('}\n\end{deluxetable*}\n')
        else:
            ftail.write('}\n\end{deluxetable}\n')
        f.close()
        if onlydata2file:
            fheader.close()
            ftail.close()
        self.outputseparator = outputseparator_bkp
        self.outputadd2lineend = outputadd2lineend_bkp

    def save2plaintextable(self,filename,cols=None,colnamehash={},colalignhash={},
                               onlydata2file=False,twocolumnflag=False,tablelabelstring='tab:',**args):
        outputseparator_bkp = self.outputseparator
        outputadd2lineend_bkp = self.outputadd2lineend
        self.outputseparator=' & '
        self.outputadd2lineend='\\\\'
        #self.outputundefined=outputundefined
        colline=''
        colnumberline=''
        formatstring = ''
        cols = self.__cols2use__(cols=cols)
        for i in range(len(cols)):
            if colline != '':
                colline += ' & '
                colnumberline += ' & '
            if cols[i] in colnamehash:
                colname=colnamehash[cols[i]]
            else:
                colname=cols[i]
            if cols[i] in colalignhash:
               formatstring +=colalignhash[cols[i]]
            else:
                if self.colinfo[cols[i]]['type']=='s':
                    formatstring += 'l'
                else:
                    formatstring += 'r'

            colline += '%s' % colname
            colnumberline += '(%d)' % i
        f = open(filename,'w')
        if onlydata2file:
            fheader = open(filename+'.header','w')
            ftail = open(filename+'.tail','w')
        else:
            fheader=f
            ftail=f
        if twocolumnflag:
            fheader.write('\\begin{table*}\n')
        else:
            fheader.write('\\begin{table}\n')
        fheader.write('\\tabletypesize{\\scriptsize}\n')
        fheader.write('\\centering\n')
        fheader.write('\\begin{tabular}{%s}\n' % formatstring)
        fheader.write('\\hline\n')
        fheader.write(colline+'\\\\\n')
        fheader.write(colnumberline+'\\\\\n')
        fheader.write('\\hline\n')
        self.printtxttable(file=f,printheader=0,cols=cols,**args)
        ftail.write('\\hline\n')
        ftail.write('\\end{tabular}\n')
        ftail.write('\\caption{\\label{%s}}\n' % tablelabelstring)
        if twocolumnflag:
            ftail.write('\\end{table*}\n')
        else:
            ftail.write('\\end{table}\n')
        
        f.close()
        if onlydata2file:
            fheader.close()
            ftail.close()
        self.outputseparator = outputseparator_bkp
        self.outputadd2lineend = outputadd2lineend_bkp

    def save2db(self,keys=None):
        """
        Currently does nothing
        """
        pass


    def colvalue2string(self,key,col,verbose=True):
        if self.getentry(key,col, verbose=verbose)==None:
            return('%s' % self.outputundefined)
        else:
            if col in self.colinfo:
                return(self.colinfo[col]['format'] % (self.getentry(key,col)))
            else:
                # That's a hack! remove it asap
                if type(self.getentry(key,col)) in [int,int]:
                    return('%ld' % (self.getentry(key,col,verbose=verbose)))
                elif type(self.getentry(key,col)) in [float]:
                    return('%f' % (self.getentry(key,col,verbose=verbose)))
                elif type(self.getentry(key,col)) in [str]:
                    return('%s' % (self.getentry(key,col,verbose=verbose)))
                else:
                    raise RuntimeError('ERROR: unknown type of ').with_traceback(self.getentry(key,col))


    def colheader2string(self,col):
        return('%s' % col)

    def __cols2use__(self,cols=None,showallcols=0):
        if cols!=None:
            if type(cols) is str:
                cols=[cols,]
            return(cols)
        collist=[]
        for col in self.cols:
            if self.colinfo[col]['visible'] or (showallcols and (not col.startswith('_'))):
                collist.append(col)
        return(collist)

    def emptyline(self,cols=None,showallcols=0,autoformat=None, verbose=True):
        s=''
        cols = self.__cols2use__(cols=cols,showallcols=showallcols)
        for col in cols:
            if s=='':
                s=self.outputadd2line
            else:
                s+=self.outputseparator
            if autoformat:
                if (col not in self.colinfo) and (col in self.cols2map):
                    col = self.colmap[col]
                if (col not in self.colinfo):
                    s+='%s' % self.outputundefined
                else:
                    s+=self.colinfo[col]['autoformat'] % self.outputundefined
            else:
                s+='%s' % self.outputundefined
        s+=self.outputadd2lineend
        return(s)

    def outputline(self,key,cols=None,showallcols=0,autoformat=None, verbose=True):
        s=''
        cols = self.__cols2use__(cols=cols,showallcols=showallcols)
        for col in cols:
            if s=='':
                s=self.outputadd2line
            else:
                s+=self.outputseparator
            if autoformat:
                if (col not in self.colinfo) and (col in self.cols2map):
                    col = self.colmap[col]
                if (col not in self.colinfo):
                    s+=self.colvalue2string(key,col, verbose=verbose)
                else:
                    if self.colinfo[col]['latexphantomflag']:
                        stemp = self.colinfo[col]['autoformat'] % self.colvalue2string(key,col, verbose=verbose)
                        if re.search('^\s+$',stemp) or (stemp == self.outputundefined) or self.getentry(key,col)==None:
                            pass
                        else:
                            stemp = re.sub(' ','\phantom{%s}' % self.phantomfillchar,stemp)
                        s+=stemp
                    else:
                        s+=self.colinfo[col]['autoformat'] % self.colvalue2string(key,col, verbose=verbose)
            else:
                s+=self.colvalue2string(key,col)
        s+=self.outputadd2lineend
        return(s)

    def outputheader(self,cols=None,showallcols=0,autoformat=None):
        s=''
        cols = self.__cols2use__(cols=cols,showallcols=showallcols)
        for col in cols:
            colname = col
            if s=='':
                s=self.outputadd2header
            else:
                s+=self.outputseparator
            if autoformat:
#                if (not self.colinfo.has_key(col)) and (col in self.cols2map):
#                    col = self.colmap[col]
#                s+=self.colinfo[col]['autoformat'] % self.colheader2string(colname)
                if (col not in self.colinfo) and (col in self.cols2map):
                    col = self.colmap[col]
#                 try:
#                     s+=self.colinfo[col]['autoformat'] % self.colheader2string(colname)
#                 except: #allows for older formats (with perhaps extra columns) of files where the columns are configured
#                     print "WARNING (texttable:outputheaer): This column, %s,  is not in colinfo so it will be converted to a string.  If this is unacceptable, check your code\n" % col
#                     s+=self.colheader2string(col)
                if (col not in self.colinfo):
                    print(("WARNING (texttable.outputheaer): This column, %s,  is not in colinfo so it will be converted to a string.  If this is unacceptable, check your code\n" % col))
                    s+=self.colheader2string(colname)
                else:
                    s+=self.colinfo[col]['autoformat'] % self.colheader2string(colname)
            else:
                s+=self.colheader2string(col)
        s+=self.outputadd2lineend
        return(s)

    def applyautoformat(self,keys=None,cols=None,showallcols=0,showallrows=0,mask='default',col4mask='default'):
        cols = self.__cols2use__(cols=cols,showallcols=showallcols)
        if keys == None:
            keys = self.rowkeys(keys=keys,mask=mask,col4mask=col4mask,showallrows=showallrows)
        for col in cols:
            colname = col
            if col not in self.colinfo:
                if col in self.cols2map:
                    col = self.colmap[col]
                else:
                    continue
            maxlen=0
            if not self.colinfo[col]['latexphantomflag']:
                maxlen=max(maxlen,len(self.colheader2string(col)))
            for key in keys:
                # skip math expressions
                if re.search('^\$.*\$$',self.colvalue2string(key,col)):
                    continue
                maxlen=max(maxlen,len(self.colvalue2string(key,col)))
            if self.stringleftalign and self.colinfo[col]['type'] == 's':
                self.colinfo[col]['autoformat']='%'+('-%d' % maxlen)+'s'
            else:
                self.colinfo[col]['autoformat']='%'+('%d' % maxlen)+'s'

    def printtxttable(self,file=None,keys=None,cols=None,showallcols=0,showallrows=0,
                      mask='default',col4mask='default',autoformat=1,printheader=20,linelist=None,printheader4nodata=1,extralinesbeforedata=None,verbose=True):

        if extralinesbeforedata!=None:
            for line in extralinesbeforedata:
                if line[0]!='#':
                    line = '#'+line
                if line[-1]!='\n':
                    line +='\n'
                if file!=None:
                    file.write(line)
                if linelist != None:
                    linelist.append(line)


        keys = self.rowkeys(keys=keys,mask=mask,col4mask=col4mask,showallrows=showallrows)

        if printheader4nodata == 0:
            if keys == None or len(keys)==0:
                return 0

        if autoformat!=None:
            self.applyautoformat(keys=keys,cols=cols,showallcols=showallcols,showallrows=showallrows,mask=mask,col4mask=col4mask)

        if linelist == None:
            if file==None:
                file=sys.stdout

        counter=1

        # only print the header if user wants it
        if printheader>0:
            header=self.outputheader(cols=cols,showallcols=showallcols,autoformat=autoformat)+'\n'
            if file!=None:
                file.write(header)
            if linelist != None:
                linelist.append(header)


        for key in keys:
            line=self.outputline(key,cols=cols,showallcols=showallcols,autoformat=autoformat, verbose=verbose)+'\n'
            if file!=None:
                file.write(line)
            if linelist != None:
                linelist.append(line)
            counter+=1
            # for STDOUT: reprint header every X lines if X>1, X=printheader
            if file==sys.stdout and printheader>1:
                if math.fmod(counter,printheader)==0.0:
                    file.write(header)
                    counter=1
        return(len(keys))

    def countNskipped(self,mask='default',col4mask='default'):
        return(self.allrowkeys-len(self.rowkeys(mask=mask,col4mask=col4mask)))
    def countNused(self,mask='default',col4mask='default'):
        return(len(self.rowkeys(mask=mask,col4mask=col4mask)))

    def definemask(self,mask):
        self.mask.definemask(mask)

    def rowkeys(self,keys=None,mask='default',col4mask='default',showallrows=0):
        """
        Returns all row keys matching keyword 'mask' and 'col4mask' conditions.
        If no mask or col4mask conditions are set or showallrows is set to True
        then returns all keys.
        If 'keys' is set then just returns those 'keys'.
        """
        # If keys are given the it overrides all the other options
        ## 2004/10/30 - MWV: I don't understand the purpose of this.
        ##   The analogous option for 'selectkeys' makes complete sense.
        ##    This makes none.
        ##   When would you use 'self.rowkeys(keys=keylist)'?
        ##   It just returns 'keylist' which you had in the first place.
        if keys!=None:
            if not (type(keys) is list):
                keys=[keys,]
            return(keys)

        # make a quick test if all rows are returned
        if showallrows:
            return(self.allrowkeys)
        if self.mask.parsemaskargument(mask)==(None,None,None,None):
            return(self.allrowkeys)
        # go through all rows and decide which rows to use
        keylist=[]
        for key in self.allrowkeys:
            if self.test4use(key,mask=mask,col4mask=col4mask):
                keylist.append(key)
        return(keylist)

    def selectkeys(self,keys=None,mask='default',col4mask='default'):
        # If keys are given the it overrides all the other options
        if keys==None:
            keys = self.allrowkeys
        # make a quick test if all rows are returned
        if self.mask.parsemaskargument(mask)==(None,None,None,None):
            return(keys)
        # go through all rows and decide which rows to use
        keylist=[]
        for key in keys:
            if self.test4use(key,mask=mask,col4mask=col4mask):
                keylist.append(key)
        return(keylist)

    def test4use(self,key,mask='default',col4mask='default'):
        if col4mask=='default':
            col4mask=self.col4mask
        #self.data[key][col4mask] = self.getentry(key, col4mask)
        #return(self.mask.test4use(self.data[key][col4mask],maskargument=mask))
        return(self.mask.test4use(self.getentry(key,col4mask),maskargument=mask))

    def sortkeysbycols(self,keys,cols,asstring=1,reverse=0):
        if keys == None:
            return None
        if asstring:
            self.applyautoformat(cols=cols,keys=keys)
            temp = [(self.outputline(key,cols=cols,showallcols=1,autoformat=1),key) for key in keys]
            temp.sort()
            if reverse: temp.reverse()
            return([x[1] for x in temp])
        else:
            if not (type(cols) is str):
                raise RuntimeError("sorting can only be done with one col if asstring=0")
            if not (self.colinfo[cols]['type'] in ['d','f']):
                print(('ERROR: cannot sort this type (%s) as number!' % (self.colinfo[cols]['type'])))
                sys.exit(0)
            temp = [(self.getentry(key,cols),key) for key in keys]
            temp.sort()
            if reverse: temp.reverse()
            return([x[1] for x in temp])

    def sortbycols(self,cols,asstring=1,reverse=0):
        if asstring:
            self.applyautoformat(cols=cols,showallrows=1)
            temp = [(self.outputline(key,cols=cols,showallcols=1,autoformat=1),key) for key in self.allrowkeys]
            temp.sort()
            if reverse: temp.reverse()
            self.allrowkeys=[x[1] for x in temp]
        else:
            if not (type(cols) is str):
                raise RuntimeError("sorting can only be done with one col if asstring=0")
            temp = [(self.getentry(key,cols),key) for key in self.allrowkeys]
            temp.sort()
            if reverse: temp.reverse()
            self.allrowkeys=[x[1] for x in temp]

    def sortkeysbycols2(self,keys=None,cols='all',reverse=0):
        if keys == None:
            return None
        if keys == 'all':
            keys = list(self.data.keys())
        if keys == 'allunmasked':
            keys = self.rowkeys()
 #       print 'FFFFFFFFFF', type(cols), n_cols
        if type(cols) is str:
            cols = [cols]
#        print 'GGGGGGGGG', type(cols)
        n_cols = len(cols)
        temp=[]
        self.applyautoformat(cols=cols,keys=keys)

        for key in keys:
            keytemp=[]
            for col in cols:
                keytemp.append(self.data[key][col])
            keytemp.append(key)
            temp.append(keytemp)
        temp.sort()
        if reverse: temp.reverse()
        return([x[n_cols] for x in temp])

    def searchcols(self,fctn,cols,*args):
        """
        This function passes ([getentry(key,cols[]),*args) to fctn, which must return either 0 or 1
        e.g. searchcols(somefctn,['ra','dec'],ra0,dec0) passes
        (getentry(key,'ra'),getentry(key,'dec'),ra0,dec0) to somefctn, with key in allrowkeys
        """
        goodkeys=[]
        for key in self.allrowkeys:
            temp=[]
            for c in cols:
                temp.append(self.getentry(key,c))
            for i in range(len(args)):
                temp.append(args[i])

            if fctn(*tuple(temp)):
                goodkeys.append(key)
        return(goodkeys)


    def searchcolkeys(self,fctn,cols,keylist,*args):
        """
        This function passes ([getentry(key,cols[]),*args) to fctn, which must return either 0 or 1
        e.g. searchcols(somefctn,['ra','dec'],ra0,dec0) passes
        (getentry(key,'ra'),getentry(key,'dec'),ra0,dec0) to somefctn, with key in allrowkeys
        """
        goodkeys=[]
        for key in keylist:
            temp=[]
            for c in cols:
                temp.append(self.getentry(key,c))
            for i in range(len(args)):
                temp.append(args[i])

            if fctn(*tuple(temp)):
                goodkeys.append(key)
        return(goodkeys)


    def key4minval(self,col,keys = None):
        if keys == None:
            keys = self.rowkeys()
        if len(keys)<1:
            return None
        keymin   = keys[0]
        minval   = self.getentry(keymin,col)
        for key in keys:
            if self.getentry(key,col)!=None:
                if minval == None or self.getentry(key,col)<minval:
                    minval = self.getentry(key,col)
                    keymin = key
        return keymin

    def key4maxval(self,col,keys = None):
        if keys == None:
            keys = self.rowkeys()
        if len(keys)<1:
            return None
        keymax   = keys[0]
        maxval = self.getentry(keymax,col)
        for key in keys:
            if self.getentry(key,col)!=None:
                if maxval == None or self.getentry(key,col)>maxval:
                    maxval = self.getentry(key,col)
                    keymax = key
        return keymax

    def minentry(self,col,keys = None):
        key = self.key4minval(col,keys = keys)
        if key == None:
            return None
        else:
            return self.getentry(key,col)

    def maxentry(self,col,keys = None):
        key = self.key4maxval(col,keys = keys)
        if key == None:
            return None
        else:
            return self.getentry(key,col)

    def median(self,col,keys=None,default = None,fraction=0.5):
        if keys == None:
            keys = self.rowkeys()
        if keys == None or len(keys)<1:
            return None,None
        keys = self.sortkeysbycols(keys,col,asstring=0)
        medianindex = int(len(keys)*fraction)
        key = keys[medianindex]
        if medianindex<len(keys)-1:
            w1 = 1.0-(len(keys)*fraction-int(len(keys)*fraction))
            w2 = 1.0 - w1
            val1 = self.getentry(keys[medianindex],col)
            val2 = self.getentry(keys[medianindex+1],col)
            if w2==0.0:
                median = val1
            elif w1==0.0:
                median = val2
            else:
                median = (val1*w1+val2*w2)/(w1+w2)
            #print len(keys)*fraction,int(len(keys)*fraction),w1,w2,val1,val2,median
        else:
            median = self.getentry(key,col)
        #self.printtxttable(cols=['flux_c','dflux_c'],keys=keys)
        #print median
        return median,key

    def calcaverage(self,keys,col,default = None):
        sum = 0.0
        N = 0
        for key in keys:
            if self.getentry(key,col) == None:
                continue
            sum += self.getentry(key,col)
            N += 1
        average = default
        if (N>0):
            average = sum/N
        return average

    def calcaverage_sigmacut(self,col,mu,sigma,keys=None,Nsigma=0.0,verbose=0,returnusekeys=False):
        if keys == None:
            keys = self.rowkeys()
        if not ('__temp4calcaverage_skipflag' in self.cols):
            self.configcols(['__temp4calcaverage_skipflag'],'d',visible=0)
            for key in keys:
                self.setentry(key,'__temp4calcaverage_skipflag',1)
        Nchanged = 0
        Nskipped=Nused=0
        sum=0.0
        for key in keys:
            skipflag = 0
            if sigma>0.0 and Nsigma>0.0:
                delta = mu - self.getentry(key,col)
                #print 'FFFF', delta,Nsigma,sigma
                if delta * delta > Nsigma * Nsigma * sigma * sigma:
                    skipflag = 1
                if skipflag != self.getentry(key,'__temp4calcaverage_skipflag'):
                    Nchanged += 1
            self.setentry(key,'__temp4calcaverage_skipflag',skipflag)
            if skipflag:
                Nskipped += 1
                continue
            sum += self.getentry(key,col)
            Nused+=1
        if Nused<1:
            if verbose:
                print('WARNING: no data to calculate average!')
            if returnusekeys:
                return(None,None,None,None,None,[])
            else:
                return(None,None,None,None,None)
        mu = sum/Nused
        sigma = None
        sum2=0.0
        if returnusekeys:
            usekeys=[]
        if Nused>1:
            Nused2=0
            for key in keys:
                if self.getentry(key,'__temp4calcaverage_skipflag'):
                    continue
                if returnusekeys:
                    usekeys.append(key)
                diff = (mu - self.getentry(key,col))
                sum2 += diff*diff
                Nused2 += 1
            if Nused2 != Nused:
                raise RuntimeError('ERROR: BUG!!!!!! Nused2 != Nused (%d!=%d) in calcaverage_sigmacut' % (Nused2,Nused))
            sigma = math.sqrt(1.0/(Nused-1)*sum2)
        if returnusekeys:
            return(mu,sigma,Nchanged,Nused,Nskipped,usekeys)
        else:
            return(mu,sigma,Nchanged,Nused,Nskipped)

    def calcaverage_sigmacutloop(self,col,keys=None,Nsigma=0.0,Nitmax=10,verbose=0,returnusekeys=False):
        self.configcols(['__temp4calcaverage_skipflag'],'d',visible=0)
        if keys == None:
            keys = self.rowkeys()
        for key in keys:
            self.setentry(key,'__temp4calcaverage_skipflag',1)
        sigma = 0.0
        mu = 0
        Nchanged=Nused=Nskipped=0
        i = 0
        converged = 0
        while ((i<Nitmax) or (Nitmax==0)) and (not converged):
            if returnusekeys:
                (mu,sigma,Nchanged,Nused,Nskipped,usekeys) = self.calcaverage_sigmacut(col,mu,sigma,keys=keys,Nsigma=Nsigma,verbose=verbose,returnusekeys=returnusekeys)
            else:
                (mu,sigma,Nchanged,Nused,Nskipped) = self.calcaverage_sigmacut(col,mu,sigma,keys=keys,Nsigma=Nsigma,verbose=verbose,returnusekeys=returnusekeys)
            # Only do a sigma cut if wanted
            if verbose>=2:
                if mu != None and sigma != None:
                    print(("i:%d  mu:%f  sigma:%f Nchanged:%d Nused:%d Nskipped:%d" % (i,mu,sigma,Nchanged,Nused,Nskipped)))
                else:
                    print(("i:",i," mu:",mu," sigma:",sigma," Nchanged:",Nchanged," Nused:",Nused," Nskipped:",Nskipped))
                    break
            if Nsigma == 0.0:
                break
            if (i>0) and (Nchanged==0):
                converged=1
                break
            if sigma==None or sigma==0.0:
                break
            i+=1
        if returnusekeys:
            return(mu,sigma,Nchanged,Nused,Nskipped,usekeys)
        else:
            return(mu,sigma,Nchanged,Nused,Nskipped)

    def calcaverage_errorcut(self,keys,col,mu,mu_err,sigma,dcol = None,dcol2 = None,Nsigma=0.0,medianflag=False,verbose=0,skipcol='__temp4calcaverage_skipflag'):
        if not (skipcol in self.cols):
            self.configcols([skipcol],'d',visible=0)
            self.setcol2value(skipcol,1,keys=keys)
        #print skipcol
        #sys.exit(0)
        Nchanged = 0
        Nskipped=Nused=0
        C1=C2=0.0
        mu_err2=mu_err*mu_err
        if dcol2 == None :
            for key in keys:
                skipflag = 0
                if sigma>0.0 and Nsigma>0.0:
                    delta = mu - self.getentry(key,col)
                    if delta * delta > Nsigma * Nsigma * (self.getentry(key,dcol) * self.getentry(key,dcol)):
                        if verbose>2: print(("I am getting skipped", abs(delta), Nsigma*self.getentry(key,dcol), (abs(delta)> Nsigma*self.getentry(key,dcol))))
                        if verbose>2:
                            print(("I am getting skipped0", mu,self.getentry(key,col),self.getentry(key,dcol)))
                            print(("I am getting skipped1", abs(delta), Nsigma*self.getentry(key,dcol), (abs(delta)> Nsigma*self.getentry(key,dcol))))
                        skipflag = 1
                    else:
                        if verbose>2:
                            print(("I DONT get skipped0", mu,self.getentry(key,col),self.getentry(key,dcol)))
                            print(("I DONT get skipped1", abs(delta), Nsigma*self.getentry(key,dcol), (abs(delta)> Nsigma*self.getentry(key,dcol))))
                    if skipflag != self.getentry(key,skipcol):
                        Nchanged += 1
                self.setentry(key,skipcol,skipflag)
                if skipflag:
                    Nskipped += 1
                    continue
                inverror2 = 1.0/self.getentry(key,dcol)
                inverror2 = inverror2 * inverror2
                C1 += self.getentry(key,col)*inverror2
                C2 += inverror2
                Nused+=1
            if Nused<1:
                if verbose:
                    print(('WARNING: no data to calculate average of col %s!' % (col)))
                    print(('mu = %.2f, Nsigma=%.0f' % (mu,Nsigma)))
                    self.printtxttable(keys=keys,cols=[col,dcol])
                return (None,None,None,None,None,None,None)
            if medianflag and len(keys)>2:
                vals = self.col_as_list(col,keys=keys)
                vals.sort()
                index = int(len(vals)*0.5)
                mu = vals[index]
            else:
                mu = C1/C2                                                    
            mu_err = math.sqrt(1.0/C2)
            sigma = None
            C1=0.0
            X2norm=0.0
            if Nused>1:
                Nused2=0
                for key in keys:
                    if self.getentry(key,skipcol):
                        continue
                    diff = (mu - self.getentry(key,col))
                    Cdummy = C1
                    C1 += diff*diff
                    #print 'NNN',mu,self.getentry(key,col),diff,self.getentry(key,dcol),diff*diff/(self.getentry(key,dcol)*self.getentry(key,dcol))
                    X2norm += diff*diff/(self.getentry(key,dcol)*self.getentry(key,dcol))
                    Nused2 += 1
                if Nused2 != Nused:
                    raise RuntimeError('ERROR: BUG!!!!!! Nused2 != Nused (%d!=%d) in calcaverage_errorcut' % (Nused2,Nused))
                sigma = math.sqrt(1.0/(Nused-1)*C1)
                X2norm = X2norm/(Nused2-1)

        elif dcol2 != None :    #quicker if no verbose
            for key in keys:
                skipflag = 0
                if sigma>0.0 and Nsigma>0.0:
                    delta = mu - self.getentry(key,col)
                    if delta * delta > Nsigma * Nsigma * self.getentry(key,dcol2):
                        if verbose>2: print(("I am getting skipped", abs(delta), Nsigma*math.sqrt(self.getentry(key,dcol2)), (abs(delta)> Nsigma*math.sqrt(self.getentry(key,dcol2)))))
                        if verbose>2:
                            print(("I am getting skipped0", mu,self.getentry(key,col),math.sqrt(self.getentry(key,dcol2))))
                            print(("I am getting skipped1", abs(delta), Nsigma*math.sqrt(self.getentry(key,dcol2)), (abs(delta)> Nsigma*math.sqrt(self.getentry(key,dcol2)))))
                        skipflag = 1
                    else:
                        if verbose>2:
                            print(("I DONT get skipped0", mu,self.getentry(key,col),math.sqrt(self.getentry(key,dcol2))))
                            print(("I DONT get skipped1", abs(delta), Nsigma*math.sqrt(self.getentry(key,dcol2)), (abs(delta)> Nsigma*math.sqrt(self.getentry(key,dcol2)))))
                    if skipflag != self.getentry(key,skipcol):
                        Nchanged += 1
                self.setentry(key,skipcol,skipflag)
                if skipflag:
                    Nskipped += 1
                    continue
                inverror2 = 1.0/self.getentry(key,dcol2)
                C1 += self.getentry(key,col)*inverror2
                C2 += inverror2
                Nused+=1
            if Nused<1:
                if verbose:
                    print(('WARNING: no data to calculate average of col %s!' % (col)))
                    print(('mu = %.2f, Nsigma=%.0f' % (mu,Nsigma)))
                    self.printtxttable(keys=keys,cols=[col,dcol2])
                return (None,None,None,None,None,None,None)

            mu_err = math.sqrt(1.0/C2)
            sigma = None
            C1=0.0
            X2norm=0.0
            if Nused>1:
                Nused2=0
                for key in keys:
                    if self.getentry(key,skipcol):
                        continue
                    diff = (mu - self.getentry(key,col))
                    Cdummy = C1
                    C1 += diff*diff
                    X2norm += diff*diff/(self.getentry(key,dcol2))
                    Nused2 += 1
                if Nused2 != Nused:
                    raise RuntimeError('ERROR: BUG!!!!!! Nused2 != Nused (%d!=%d) in calcaverage_errorcut' % (Nused2,Nused))
                sigma = math.sqrt(1.0/(Nused-1)*C1)
                X2norm = X2norm/(Nused2-1)

        return(mu,mu_err,sigma,X2norm,Nchanged,Nused,Nskipped)

    def calcaverage_errorcutloop(self,keys,col,dcol=None,dcol2 = None,Nsigma=0.0,Nitmax=10,medianflag4firstiteration=False,returnusedkeys=False,skipcol='__temp4calcaverage_skipflag',verbose=0):
        if dcol == None and dcol2 == None :
            print('WARNING! GIVE AT LEAST ONE ERROR COLUMN IN CALCAVERAGE_ERRORCUTLOOP, SQUARED OR NOT')
            return False
        if (not (skipcol in self.cols)) or self.colinfo[skipcol]['type']!='d':
            self.configcols([skipcol],'d',visible=0)
        self.setcol2value(skipcol,1,keys=keys)
        #self.configcols(['__temp4calcaverage_skipflag'],'d',visible=0)
        #for key in keys:
        #    self.setentry(key,'__temp4calcaverage_skipflag',1) 
        sigma = 0.0
        mu = mu_err = 0.0
        Nchanged=Nused=Nskipped=0
        i = 0
        converged = 0
        while ((i<Nitmax) or (Nitmax==0)) and (not converged):
            if i==0:
                Nsigma_temp = 2*Nsigma
                medianflag = medianflag4firstiteration                
            else:
                Nsigma_temp = Nsigma
                medianflag = False
            (mu,mu_err,sigma,X2norm,Nchanged,Nused,Nskipped) = self.calcaverage_errorcut(keys,col,mu,mu_err,sigma,dcol,dcol2,Nsigma=Nsigma_temp,medianflag=medianflag,skipcol=skipcol,verbose=verbose)
            # Only do a sigma cut if wanted
            if verbose>=1 and mu!=None and sigma!=None:
                print(("mu:%f mu_err:%f sigma:%f X2norm:%f Nchanged:%d Nused:%d Nskipped:%d" % (mu,mu_err,sigma,X2norm,Nchanged,Nused,Nskipped)))
            if Nsigma == 0.0:
                break
            if (i>0) and (Nchanged==0):
                converged=1
                break
            if sigma==None or sigma==0.0:
                break
            #if verbose:
            #    print 'Nused round',i, Nused, len(keys)
            i+=1
        if returnusedkeys:
            usedkeys = self.CUT_inrange(skipcol,0,0,keys=keys)
            return(mu,mu_err,sigma,X2norm,Nchanged,Nused,Nskipped,usedkeys)
        else:
            return(mu,mu_err,sigma,X2norm,Nchanged,Nused,Nskipped)

    def calcweightedaverage(self,keys,col,dcol, verbose=False):
        if len(keys)<1:
            print('No Data!')
            return(None,None,None)
        sum = sumweight = 0.0

        for key in keys:
            if self.getentry(key,dcol) == 0.0:
                if verbose:
                    print(("Skipping key %d b/c weight, %s, have 0.0 value\n" % (key, dcol)))
                continue
            weight = 1.0/(self.getentry(key,dcol)*self.getentry(key,dcol))
            sum += weight * self.getentry(key,col)
            sumweight+= weight
        mean = sum/sumweight
        error = math.sqrt(1.0/sumweight)
        X2 = 0.0
        for key in keys:
            if self.getentry(key,dcol) == 0.0:
                continue
            diffnorm = (mean-self.getentry(key,col))/self.getentry(key,dcol)
            X2 += diffnorm*diffnorm
        X2norm = X2/len(keys)
        return(mean,error,X2norm)

    # This routine doesn't really belong into texttable, it belongs into
    # lctableclass. Since we need it for some eventsstatsclass, which uses
    # only texttables and not lctableclass, we have to make a hack right now
    # and put it here
    def selectlckeys(self,keys=None,photcodes=None,dophottypes=None,onlyforced=0,onlynonforced=0,MJDmin=None,MJDmax=None,skipbadim=0,onlybadim=0,Nmaskmax=None,skipifpeakmasked=None):
        if keys == None:
            keys = self.rowkeys()
        FORCEDPHOTFLAG = 0x10
        # check the passed photcodes: make it a list, and make sure it is in hexadecimal format
        if type(photcodes) is int:
            photcodes = [photcodes,]
        if photcodes!=None:
            for i in range(len(photcodes)):
                photcodes[i] = hex2int(photcodes[i]) & 0xffff

        # check the passed dophot types: make it a list
        if type(dophottypes) is int:
            dophottypes = [dophottypes,]

        self.configcols(['type'],'x')

        newkeys = []
        for lckey in keys:
            # skip bad images if wanted
            if skipbadim:
                if self.asint(lckey,'badim'):
                    continue

            # skip good images if wanted
            if onlybadim:
                if self.asint(lckey,'badim')==0:
                    continue


            # check if it is forced
            if onlyforced:
                if not (self.getentry(lckey,'type') & FORCEDPHOTFLAG):
                    continue

            # check if it is unforced
            if onlynonforced:
                if (self.getentry(lckey,'type') & FORCEDPHOTFLAG):
                    continue

            # check if it is the correct dophot type
            if dophottypes != None:
                foundflag = 0
                lct = self.getentry(lckey,'type') & 0xf
                for t in dophottypes:
                    if t == lct:
                        foundflag = 1
                        break

                if not foundflag:
                    # bad type
                    continue

            # check if it is the correct photcode
            if photcodes != None:
                lcpc = self.getentry(lckey,'photcode')
                if lcpc == None:
                    continue
                lcpc = hex2int(lcpc) & 0xffff
                foundflag = 0
                for pc in photcodes:
                    if pc == lcpc:
                        foundflag = 1
                        break
                if not foundflag:
                    # bad photcode
                    continue
            # check for MJD
            if MJDmin is not None or MJDmax is not None:
                mjd = self.getentry(lckey,'MJD')
                if mjd is not None:
                    mjd = float(mjd)
                    if MJDmin is not None and mjd<MJDmin:
                        continue
                    if MJDmax is not None and mjd>MJDmax:
                        continue
            if Nmaskmax != None:
                Nmask = self.getentry(lckey,'Nmask')
                if Nmask == None:
                    continue
                Nmask = int(Nmask)
                if Nmask>int(Nmaskmax):
                    continue
            if skipifpeakmasked !=None and skipifpeakmasked!=0:
                flag = self.getentry(lckey,'Nmask')
                if flag == None:
                    continue
                if hex2int(flag) & 0x8:
                    continue



            newkeys.append(lckey)
        return newkeys

    def col_as_list(self,col,keys=None):
        if keys == None: keys = self.rowkeys()
        list = [self.getentry(key,col) for key in keys]
        return list

    def col_asint_list(self,col,keys=None):
        if keys == None: keys = self.rowkeys()
        intlist = [self.asint(key,col) for key in keys]
        return intlist

    def col_asfloat_list(self,col,keys=None):
        if keys == None: keys = self.rowkeys()
        list = [self.asfloat(key,col) for key in keys]
        return list

    def col_ashex_list(self,col,keys=None):
        if keys == None: keys = self.rowkeys()
        list = [self.ashex(key,col) for key in keys]
        return list

    def col_asstring_list(self,col,keys=None):
        if keys == None: keys = self.rowkeys()
        list = [self.asstring(key,col) for key in keys]
        return list

    def cols_as_dict(self,key,cols):
        if type(cols) is str:
            cols = [cols,]
        dict = {}
        for col in cols:
            dict[col]=self.getentry(key,col)
        return dict

    def row_as_list(self,key,cols=None):
        if cols == None:
            cols = self.cols
        else:
            if type(cols) is str:
                cols = [cols,]
        list = [self.getentry(key,col) for col in cols]
        return list

    def row_as_dict(self,key,cols=None):
        if cols == None:
            cols = self.cols
        else:
            if type(cols) is str:
                cols = [cols,]
        dict = {}
        for col in cols:
            dict[col]=self.getentry(key,col)
        return dict

    def CUT_more_than(self, col, maxval, keys = None, CUTlimit=False):
        if keys == None: keys = self.rowkeys()
        if maxval == None:
            return keys
        newkeys = []
        for key in keys:
            val = self.getentry(key,col)
            if (val==None) or (val!=None and val>maxval):
                continue
            if CUTlimit and val==maxval:
                continue
            newkeys.append(key)

            #if val!=None and val<=maxval:
            #    newkeys.append(key)
        return(newkeys)

    def CUT_less_than(self, col, minval, keys = None, CUTlimit=False):
        if keys == None: keys = self.allrowkeys
        if minval == None:
            return keys
        newkeys = []
        for key in keys:
            val = self.getentry(key,col)
            if (val==None) or (val!=None and val<minval):
                continue
            if CUTlimit and val==minval:
                continue
            newkeys.append(key)
            #if val!=None and val>=minval:
            #    newkeys.append(key)
        return(newkeys)

    def CUT_inrange(self, col, minval, maxval, keys=None, excludelowerlimit=False,excludeupperlimit=False):
        if keys == None: keys = self.allrowkeys
        if minval == None and maxval == None:
            return keys
        newkeys = []
        for key in keys:
            val = self.getentry(key,col)
            if val!=None:
                if minval!=None and val<minval:
                    continue
                if maxval!=None and val>maxval:
                    continue
                if excludelowerlimit and val==minval:
                    continue
                if excludeupperlimit and val==maxval:
                    continue
                newkeys.append(key)
        return(newkeys)

    def CUT_outrange(self, col, minval, maxval, keys=None):
        if keys == None: keys = self.allrowkeys
        if minval == None and maxval == None:
            return keys
        newkeys = []
        for key in keys:
            val = self.getentry(key,col)

            if val!=None:
                if minval!=None and val<minval:
                    newkeys.append(key)
                elif maxval!=None and val>maxval:
                    newkeys.append(key)
        return(newkeys)


    def CUT_none_vals(self, cols, keys = None):
        if keys == None: keys = self.allrowkeys
        newkeys = []
        if type(cols) is str:
            for key in keys:
                if self.getentry(key,cols)!=None:
                    newkeys.append(key)
        else:
            for key in keys:
                NoneFlag=False
                for col in cols:
                    NoneFlag |= (self.getentry(key,col)==None)
                if not NoneFlag:
                    newkeys.append(key)
        return(newkeys)

    def CUT_zero_vals(self, col, keys=None):
        if keys == None: keys = self.allrowkeys
        newkeys = []
        for key in keys:
            val = self.getentry(key,col)
            if val != 0.0:
                newkeys.append(key)
        return newkeys

    def get_binned_keys(self, bins, col, keys=None, minkeys4bin=None, binstepsize=0.1,minval=None, maxkeys4bin=None):
        if keys == None: keys = self.rowkeys()
        Nbins = len(bins)
        binnedkeys = list(range(Nbins))
        binds = list(range(Nbins))

        for bin in binds:
            if bin == 0:
                #print bins[bin]
                if minval !=None:
                    bkeys = self.CUT_inrange(col=col, keys=keys, maxval=bins[bin], minval=minval)
                else:
                    bkeys = self.CUT_more_than(col=col, maxval=bins[bin], keys=keys)
                if minkeys4bin:
                    i=0
                    while len(bkeys)>int(minkeys4bin*.1) and len(bkeys) < minkeys4bin:
                        bins[bin:] = [b+binstepsize for b in bins[bin:]]
                        if minval != None:
                            bkeys = self.CUT_inrange(col=col, maxval=bins[bin], minval=minval, keys=keys)
                        else:
                            bkeys = self.CUT_more_than(col=col, maxval=bins[bin], keys=keys)
                        i+=1
                        if i>10:
                            break
                if maxkeys4bin:

                    i=0
                    while len(bkeys) > maxkeys4bin:
                        bins[bin:] = [b-binstepsize for b in bins[bin:]]
                        if minval != None:
                            bkeys = self.CUT_inrange(col=col, maxval=bins[bin], minval=minval, keys=keys)
                        else:
                            bkeys = self.CUT_more_than(col=col, maxval=bins[bin], keys=keys)
                        i+=1
                        if i>10:
                            break

                binnedkeys[bin] = bkeys
                keys = [key for key in keys if key not in bkeys]
            else:
                bkeys = self.CUT_inrange(col=col, maxval=bins[bin], minval=bins[bin-1], keys=keys)

                if minkeys4bin:
                    i=0
                    while len(bkeys) > int(minkeys4bin*.1) and len(bkeys) < minkeys4bin:
                        bins[bin:] = [b+binstepsize for b in bins[bin:]]
                        #print 'last bin', bins[Nbins-1], 'current bin', bins[bin], 'n keys', len(bkeys)
                        bkeys = self.CUT_inrange(col=col, maxval=bins[bin], minval=bins[bin-1], keys=keys)
                        i+=1
                        if i>10:
                            break

                if maxkeys4bin:
                    i=0

                    while len(bkeys) > maxkeys4bin:
                        bins[bin:] = [b-binstepsize for b in bins[bin:]]
                        if minval != None:
                            bkeys = self.CUT_inrange(col=col, maxval=bins[bin], minval=minval, keys=keys)
                        else:
                            bkeys = self.CUT_more_than(col=col, maxval=bins[bin], keys=keys)
                        i+=1
                        if i>10:
                            break

                binnedkeys[bin] = bkeys
                keys = [key for key in keys if key not in bkeys]

        return bins, binnedkeys

    def asfloat(self,key,col,default=None):
        try:
            val = self.getentry(key,col)
            if val == None: return(default)
            if type(val) is float: return(val)
            return(float(val))
        except:
            return default

    def asint(self,key,col,default=None):
        try:
            val = self.getentry(key,col)
            if val == None: return(default)
            if type(val) is int: return(val)
            return(int(eval(val)))
        except:
            return default

    def ashex(self,key,col,default=None):
        try:
            val = self.getentry(key,col)
            if val == None: return(default)
            #if type(val) is types.StringType:
            return(int(val,16))
        except:
            return default

    def asstring(self,key,col,default=None,format='%s'):
        try:
            val = self.getentry(key,col)
            if val == None: return(default)
            #if type(val) is types.StringType:
            return(format % val)
        except:
            return default

    def keys_undefcolval(self,cols,keys = None):
        if not (type(cols) is list): cols=[cols,]
        if keys == None: keys = self.rowkeys()
        undefkeys = []
        for key in keys:
            undef = 0
            for col in cols:
                if self.getentry(key,col)==None:
                    undef = 1
                    break
            if undef == 1:
                undefkeys.append(key)
        return(undefkeys)

    def keys_defcolval(self,cols,keys = None):
        if not (type(cols) is list): cols=[cols,]
        if keys == None: keys = self.rowkeys()
        defkeys = []
        for key in keys:
            undef = 0
            for col in cols:
                if self.getentry(key,col)==None:
                    undef = 1
                    break
            if undef == 0:
                defkeys.append(key)
        return(defkeys)

    def colMin(self, col, keys=None):
        if keys == None: keys = self.rowkeys()
        minval = self.getentry(keys[0], col)
        for k in keys:
            if self.getentry(k, col) < minval:
                minval = self.getentry(k, col)
        return minval

    def colMax(self, col, keys=None):
        if keys == None: keys = self.rowkeys()
        maxval = self.getentry(keys[0], col)
        for k in keys:
            if self.getentry(k, col) > maxval:
                maxval = self.getentry(k, col)
        return maxval


    def mkds9regionlist(self,xcol,ycol,keys=None,color='red',symbol='circle',WCSflag=False,save2file=None):
        """
        symbol: circle, box, diamond, cross, X, boxcircle
        colors: white, black, red, green, blue, cyan, magenta, yellow
        """
        if keys == None: keys = self.allrowkeys
        out = ['global color=%s font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source' % color]
        if WCSflag:
            out.append('fk5')
        for key in keys:
            out.append('point(%s,%s) # point=%s color=%s' % (self.asstring(key,xcol),self.asstring(key,ycol),symbol,color))
        if save2file!=None:
            f = open(save2file,'w')
            for line in out:
                f.write(line+'\n')
            f.close()
        return(out)

    def straightline(self,keys,xcol,ycol,dycol=None,dycol2=None,dyval=1.0):

        if (len(keys)<2):
            return(False,0,0,0,0)

        # Bevington, page 104
        # y(x) = a + bx
        if dycol==None and dycol2==None:
            sumx2=sumx=sumy=sumxy=delta=0.0
            for key in keys:
                sumx2 += self.getentry(key,xcol)*self.getentry(key,xcol)
                sumx  += self.getentry(key,xcol)
                sumy  += self.getentry(key,ycol)
                sumxy += self.getentry(key,xcol)*self.getentry(key,ycol)
            delta = len(keys)*sumx2 - sumx*sumx
            offset_a     = 1.0/delta * (sumx2*sumy - sumx*sumxy)
            slope_b      = 1.0/delta * (len(keys)*sumxy - sumx*sumy)

            if dyval==None:
                counter=0
                diffylist = numpy.zeros((len(keys)))
                for key in keys:
                    diffylist[counter]=self.getentry(key,xcol)*slope_b+offset_a - self.getentry(key,ycol)
                    counter+=1
                dyval = numpy.std(diffylist)

            offset_a_err = math.sqrt(dyval*dyval/delta*sumx2)
            slope_b_err  = math.sqrt(len(keys)*dyval*dyval/delta)


        elif dycol != None and dycol2 == None:
            sumx2=sumx=sumy=sumxy=sumis2=delta=0.0
            for key in keys:
                inversesigma2 = 1.0/(self.getentry(key,dycol)*self.getentry(key,dycol))
                sumis2 += inversesigma2
                sumx2 += inversesigma2 * self.getentry(key,xcol)*self.getentry(key,xcol)
                sumx  += inversesigma2 * self.getentry(key,xcol)
                sumy  += inversesigma2 * self.getentry(key,ycol)
                sumxy += inversesigma2 * self.getentry(key,xcol)*self.getentry(key,ycol)
            delta = sumis2*sumx2 - sumx*sumx
            offset_a     = 1.0/delta * (sumx2*sumy - sumx*sumxy)
            slope_b      = 1.0/delta * (sumis2*sumxy - sumx*sumy)
            offset_a_err = math.sqrt(1.0/delta * sumx2)
            slope_b_err  = math.sqrt(1.0/delta * sumis2)
        elif dycol2 != None:
            sumx2=sumx=sumy=sumxy=sumis2=delta=0.0
            for key in keys:
                inversesigma2 = 1.0/(self.getentry(key,dycol2))
                sumis2 += inversesigma2
                sumx2 += inversesigma2 * self.getentry(key,xcol)*self.getentry(key,xcol)
                sumx  += inversesigma2 * self.getentry(key,xcol)
                sumy  += inversesigma2 * self.getentry(key,ycol)
                sumxy += inversesigma2 * self.getentry(key,xcol)*self.getentry(key,ycol)
            delta = sumis2*sumx2 - sumx*sumx
            offset_a     = 1.0/delta * (sumx2*sumy - sumx*sumxy)
            slope_b      = 1.0/delta * (sumis2*sumxy - sumx*sumy)
            offset_a_err = math.sqrt(1.0/delta * sumx2)
            slope_b_err  = math.sqrt(1.0/delta * sumis2)

        return(slope_b,slope_b_err,offset_a,offset_a_err)



    def straightline_errorcut(self,keys,xcol,ycol,offset,slope,x0=0.0,dycol=None,dycol2=None,Nsigma=3.0,dyval=1.0,verbose=0):

        if (len(keys)<2):
            return(False,0,0,0,0)
        if not ('__temp4calcaverage_skipflag' in self.cols):
            self.configcols(['__temp4calcaverage_skipflag'],'d',visible=0)
            for key in keys:
                self.setentry(key,'__temp4calcaverage_skipflag',1)
        Nchanged = 0
        Nskipped=Nused=0
        C1=C2=0.0

        # Bevington, page 104
        # y(x) = a + bx
        if dycol2 == None:
            sumx2=sumx=sumy=sumxy=sumis2=delta=0.0
            for key in keys:
                skipflag = 0
                if Nsigma>0.0:
                    delta = slope*(self.getentry(key,xcol)-x0)+offset - self.getentry(key,ycol)
                    if delta * delta > Nsigma * Nsigma * (self.getentry(key,dycol) * self.getentry(key,dycol)):
                        if verbose>2: print(("I am getting skipped", abs(delta), Nsigma*self.getentry(key,dycol), (abs(delta)> Nsigma*self.getentry(key,dycol))))
                        if verbose>2:
                            print(("I am getting skipped0", self.getentry(key,xcol),self.getentry(key,ycol),self.getentry(key,dycol)))
                            print(("I am getting skipped1", abs(delta), Nsigma*self.getentry(key,dycol), (abs(delta)> Nsigma*self.getentry(key,dycol))))
                        skipflag = 1
                    else:
                        if verbose>2:
                            print(("I DONT get skipped0", self.getentry(key,xcol),self.getentry(key,ycol),self.getentry(key,dycol)))
                            print(("I DONT get skipped1", abs(delta), Nsigma*self.getentry(key,dycol), (abs(delta)> Nsigma*self.getentry(key,dycol))))
                    if skipflag != self.getentry(key,'__temp4calcaverage_skipflag'):
                        Nchanged += 1
                self.setentry(key,'__temp4calcaverage_skipflag',skipflag)
                if skipflag:
                    Nskipped += 1
                    continue
                inverror2 = 1.0/(self.getentry(key,dycol)*self.getentry(key,dycol))
                C1 += self.getentry(key,ycol)*inverror2
                C2 += inverror2
                Nused+=1
            if Nused<1:
                if verbose:
                    print(('WARNING: no data to calculate average of col %s!' % (ycol)))
                    print(('slope = %.2f, Nsigma=%.0f' % (slope,Nsigma)))
                    self.printtxttable(keys=keys,cols=[ycol,dycol])
                return (None,None,None,None,None,None,None)
            if Nused==1:
                Nused2=0
                for key in keys:
                    if self.getentry(key,'__temp4calcaverage_skipflag'):
                        continue
                    data = self.getentry(key,ycol)
                    data_err = self.getentry(key,dycol)
                    Nused2 += 1
                if Nused2 != Nused:
                    raise RuntimeError('ERROR: BUG!!!!!! Nused2 != Nused (%d!=%d) in calcaverage_errorcut' % (Nused2,Nused))
                offset_a     = data
                slope_b      = 0
                offset_a_err = data_err
                slope_b_err  = 0
            if Nused>1:
                Nused2=0
                for key in keys:
                    if self.getentry(key,'__temp4calcaverage_skipflag'):
                        continue
                    inversesigma2 = 1.0/(self.getentry(key,dycol)*self.getentry(key,dycol))
                    sumis2 += inversesigma2
                    sumx2 += inversesigma2 * (self.getentry(key,xcol)-x0)*(self.getentry(key,xcol)-x0)
                    sumx  += inversesigma2 * (self.getentry(key,xcol)-x0)
                    sumy  += inversesigma2 * self.getentry(key,ycol)
                    sumxy += inversesigma2 * (self.getentry(key,xcol)-x0)*self.getentry(key,ycol)
                    Nused2 += 1
                if Nused2 != Nused:
                    raise RuntimeError('ERROR: BUG!!!!!! Nused2 != Nused (%d!=%d) in calcaverage_errorcut' % (Nused2,Nused))
                delta = sumis2*sumx2 - sumx*sumx
                offset_a     = 1.0/delta * (sumx2*sumy - sumx*sumxy)
                slope_b      = 1.0/delta * (sumis2*sumxy - sumx*sumy)
                offset_a_err = math.sqrt(1.0/delta * sumx2)
                slope_b_err  = math.sqrt(1.0/delta * sumis2)

        else :   #quicker if no verbose
            sumx2=sumx=sumy=sumxy=sumis2=delta=0.0
            for key in keys:
                skipflag = 0
                if Nsigma>0.0:
                    delta = slope*(self.getentry(key,xcol)-x0)+offset - self.getentry(key,ycol)
                    if delta * delta > Nsigma * Nsigma * (self.getentry(key,dycol2)) :
                        #if verbose==1: print "I am getting skipped", abs(delta), Nsigma*math.sqrt(self.getentry(key,dycol2)), (abs(delta)> Nsigma*math.sqrt(self.getentry(key,dycol2)))
                        if verbose>2:
                            print(("I am getting skipped0", self.getentry(key,xcol),self.getentry(key,ycol),math.sqrt(self.getentry(key,dycol2))))
                            print(("I am getting skipped1", abs(delta), Nsigma*math.sqrt(self.getentry(key,dycol2)), (delta*delta > Nsigma*Nsigma*self.getentry(key,dycol2))))
                        skipflag = 1
                    else:
                        if verbose>2:
                            print(("I DONT get skipped0",self.getentry(key,xcol),self.getentry(key,ycol),math.sqrt(self.getentry(key,dycol2))))
                            print(("I DONT get skipped1", abs(delta), Nsigma*math.sqrt(self.getentry(key,dycol2)), (delta*delta > Nsigma*Nsigma*self.getentry(key,dycol2))))
                    if skipflag != self.getentry(key,'__temp4calcaverage_skipflag'):
                        Nchanged += 1
                self.setentry(key,'__temp4calcaverage_skipflag',skipflag)
                if skipflag:
                    Nskipped += 1
                    continue
                inverror2 = 1.0/self.getentry(key,dycol2)
                Nused+=1
            if Nused<1:
                if verbose:
                    print(('WARNING: no data to calculate average of col %s!' % (ycol)))
                    print(('slope = %.2f, Nsigma=%.0f' % (slope,Nsigma)))
                    self.printtxttable(keys=keys,cols=[ycol,dycol2])
                return (None,None,None,None,None,None,None)
            if Nused==1:
                Nused2=0
                for key in keys:
                    if self.getentry(key,'__temp4calcaverage_skipflag'):
                        continue
                    data = self.getentry(key,ycol)
                    data_err = math.sqrt(self.getentry(key,dycol2))
                    Nused2 += 1
                if Nused2 != Nused:
                    raise RuntimeError('ERROR: BUG!!!!!! Nused2 != Nused (%d!=%d) in calcaverage_errorcut' % (Nused2,Nused))
                offset_a     = data
                slope_b      = 0
                offset_a_err = data_err
                slope_b_err  = 0
            if Nused>1:
                Nused2=0
                for key in keys:
                    if self.getentry(key,'__temp4calcaverage_skipflag'):
                        continue
                    inversesigma2 = 1.0/(self.getentry(key,dycol2))
                    sumis2 += inversesigma2

                    sumx2 += inversesigma2 * (self.getentry(key,xcol)-x0)*(self.getentry(key,xcol)-x0)
                    sumx  += inversesigma2 * (self.getentry(key,xcol)-x0)
                    sumy  += inversesigma2 * self.getentry(key,ycol)
                    sumxy += inversesigma2 * (self.getentry(key,xcol)-x0)*self.getentry(key,ycol)

                    Nused2 += 1
                if Nused2 != Nused:
                    raise RuntimeError('ERROR: BUG!!!!!! Nused2 != Nused (%d!=%d) in calcaverage_errorcut' % (Nused2,Nused))
                delta = sumis2*sumx2 - sumx*sumx
                offset_a     = 1.0/delta * (sumx2*sumy - sumx*sumxy)
                slope_b      = 1.0/delta * (sumis2*sumxy - sumx*sumy)
                offset_a_err = math.sqrt(1.0/delta * sumx2)
                slope_b_err  = math.sqrt(1.0/delta * sumis2)

        return(slope_b,slope_b_err,offset_a,offset_a_err,Nchanged,Nused,Nskipped)


    def straightline_errorcutloop(self,keys,xcol,ycol,x0=0.0,dycol=None,dycol2=None,Nsigma=3.0,dyval=1.0,Nitmax=10,verbose=0):
        if dycol == None and dycol2 == None :
            print('WARNING! GIVE AT LEAST ONE ERROR COLUMN IN STRAIGHTLINE_ERRORCUTLOOP, SQUARED OR NOT')
            return False
        self.configcols(['__temp4calcaverage_skipflag'],'d',visible=0)
        for key in keys:
            self.setentry(key,'__temp4calcaverage_skipflag',1)
        slope_b = offset_a = 0.0
        Nchanged=Nused=Nskipped=0
        i = 0
        converged = 0
        while ((i<Nitmax) or (Nitmax==0)) and (not converged):
            if i==0:
                Nsigma_temp = 0.0
            elif i==1:
                Nsigma_temp = 3*Nsigma
            else:
                Nsigma_temp = Nsigma
            (slope_b,slope_b_err,offset_a,offset_a_err,Nchanged,Nused,Nskipped) = self.straightline_errorcut(keys,xcol,ycol,offset_a,slope_b,x0=x0,dycol=dycol,dycol2=dycol2,Nsigma=Nsigma_temp,dyval=dyval,verbose=verbose)
            # Only do a sigma cut if wanted
            if verbose>=1 and slope_b!=None and slope_b != None:
                print(("%02d: slope_b:%f slope_b_err:%f offset_a:%f offset_a_err:%f Nchanged:%d Nused:%d Nskipped:%d" % (i,slope_b,slope_b_err,offset_a,offset_a_err,Nchanged,Nused,Nskipped)))
            if Nsigma == 0.0:
                break
            if (i>1) and (Nchanged==0):
                break
            if slope_b==None :
                break
            #if verbose:
            #    print 'Nused round',i, Nused, len(keys)
            i+=1
        return(slope_b,slope_b_err,offset_a,offset_a_err,Nchanged,Nused,Nskipped)

    def cpcol2col(self,t1,t2,cols1,cols2,keys1,keys2,checkcols=True,raiseErrorFlag=True,icols=None):
        many2many = False
        one2many = False
        if type(keys1)==list:
            many2many = True
            if len(keys1)!=len(keys2):
                if raiseErrorFlag:
                    raise RuntimeError("%d != %1" % (len(keys1),len(keys2)))
                else:
                    return(1)
        elif type(keys2)==list:
            one2many = True

        if checkcols:
            if len(cols1)!=len(cols2):
                if raiseErrorFlag:
                    raise RuntimeError("%d != %1" % (len(cols1),len(cols2)))
                else:
                    return(2)

            if not t2.colsexist(cols2):
                for col in cols2:
                    if not t2.colexist(col):print(('ERROR: %s does not exist in destination table!' % col))
                if raiseErrorFlag:
                    raise RuntimeError("column(s) do not exist in destination table")
                else:
                    return(3)

            if not t1.colsexist(cols1):
                for col in cols1:
                    if not t1.colexist(col):print(('ERROR: %s does not exist in source table!' % col))
                if raiseErrorFlag:
                    raise RuntimeError("column(s) do not exist in source table")
                else:
                    return(4)

        if icols == None: icols = list(range(len(cols1)))
        if many2many or one2many:
            ikeys = list(range(len(keys2)))
            for ikey in ikeys:
                for icol in icols:
                    if one2many:
                        # one2many
                        t2.setentry(keys2[ikey],cols2[icol],t1.getentry(keys1,cols1[icol]))
                    else:
                        # many2many
                        t2.setentry(keys2[ikey],cols2[icol],t1.getentry(keys1[ikey],cols1[icol]))
        else:
            # one2one
            for icol in icols:
                t2.setentry(keys2,cols2[icol],t1.getentry(keys1,cols1[icol]))

        return(0)


    def cpmatchedcols(self,t1,matchcol1,t2,matchcol2,cols1,cols2,keys2=None,raiseErrorFlag=True,skipifnotexist=False,addmultientries=False):
        ### cp from 1 to 2 ###

        if len(cols1)!=len(cols2):
            if raiseErrorFlag:
                raise RuntimeError("%d != %1" % (len(cols1),len(cols2)))
            else:
                return(1)

        if not t2.colsexist(cols2):
            for col in cols2:
                if not t2.colexist(col):print(('ERROR: %s does not exist in destination table!' % col))
            if raiseErrorFlag:
                raise RuntimeError("column(s) do not exist in destination table")
            else:
                return(2)

        if not t1.colsexist(cols1):
            for col in cols1:
                if not t1.colexist(col):print(('ERROR: %s does not exist in source table!' % col))
            if raiseErrorFlag:
                raise RuntimeError("column(s) do not exist in source table")
            else:
                return(3)

        irange = list(range(len(cols1)))
        if keys2==None:
            keys2=t2.rowkeys()
        for key2 in keys2:
            val2 = t2.getentry(key2,matchcol2)
            if val2!=None:
                keys1 = t1.CUT_inrange(matchcol1,val2,val2)
                if len(keys1)==0:
                    if raiseErrorFlag:
                        raise RuntimeError('ERROR! cpfluxcols: Could not find entry %s in column %s!' % (val2,matchcol1))
                    else:
                        if skipifnotexist:
                            continue
                        else:
                            return(4)
                elif len(keys1)==1:
                    for i in irange:
                        t2.setentry(key2,cols2[i],t1.getentry(keys1[0],cols1[i]))
                else:
                    if raiseErrorFlag:
                        raise RuntimeError('ERROR!!! more than one entry matching %s in column %s' % (val2,matchcol1))
                    else:
                        if addmultientries:
                            t2.setentry(key2,cols2[i],','.join(t1.col_asstring_list(cols1[i],keys=keys1)))
                        if skipifnotexist:
                            continue
                        else:
                            return(5)
            else:
                if raiseErrorFlag:
                    raise RuntimeError('ERROR!!! None value in %s' % (matchcol2))
                else:
                    return(6)
        return(0)

    def cpcols(self,t1,matchcol1,t2,matchcol2,cols1,cols2,keys2=None,keys1all=None,
               raiseErrorFlag=True,skipifnotexist=False,addmultientries=False,requireallkeys1matched=True,returnunmatchedkeys=False,skipifNone=False):
        ### cp from 1 to 2 ###

        if len(cols1)!=len(cols2):
            if raiseErrorFlag:
                raise RuntimeError("%d != %1" % (len(cols1),len(cols2)))
            else:
                return(1)

        if not t2.colsexist(cols2):
            for col in cols2:
                if not t2.colexist(col):print(('ERROR: %s does not exist in destination table!' % col))
            if raiseErrorFlag:
                raise RuntimeError("column(s) do not exist in destination table")
            else:
                return(2)

        if not t1.colsexist(cols1):
            for col in cols1:
                if not t1.colexist(col):print(('ERROR: %s does not exist in source table!' % col))
            if raiseErrorFlag:
                raise RuntimeError("column(s) do not exist in source table")
            else:
                return(3)

        irange = list(range(len(cols1)))
        if keys2==None:
            keys2=t2.rowkeys()

        keys1matched=[]

        for key2 in keys2:
            val2 = t2.getentry(key2,matchcol2)
            if val2!=None:
                keys1 = t1.CUT_inrange(matchcol1,val2,val2,keys=keys1all)
                #print 'VVVVVVVV',val2,keys1
                #if val2=='2004dt':
                #    sys.exit(0)
                keys1matched.extend(keys1)
                if len(keys1)==0:
                    if skipifnotexist:
                        continue
                    else:
                        if raiseErrorFlag:
                            t1.printtxttable()
                            t1.printtxttable(cols=[matchcol1])
                            raise RuntimeError('ERROR! cpfluxcols: Could not find entry %s in column %s!' % (val2,matchcol1))
                        else:
                            return(4)
                elif len(keys1)==1:
                    for i in irange:
                        t2.setentry(key2,cols2[i],t1.getentry(keys1[0],cols1[i]))
                else:
                    if addmultientries:
                        for i in irange:
                            t2.setentry(key2,cols2[i],','.join(t1.col_asstring_list(cols1[i],keys=keys1)))
                    else:
                        if raiseErrorFlag:
                            raise RuntimeError('ERROR!!! more than one entry matching %s in column %s' % (val2,matchcol1))
                        else:
                            if skipifNone:
                                continue
                            else:
                                return(5)
            else:
                if raiseErrorFlag:
                    raise RuntimeError('ERROR!!! None value in %s' % (matchcol2))
                else:
                    if skipifnotexist:
                        continue
                    else:
                        return(6)
        if keys1all==None:
            keys1notmatched=tools.AnotB(t1.allrowkeys,tools.unique(keys1matched))
        else:
            keys1notmatched=tools.AnotB(keys1all,tools.unique(keys1matched))


        if requireallkeys1matched:
            if len(keys1notmatched)>0:
                if raiseErrorFlag:
                    t1.printtxttable(cols=[matchcol1].extend(cols1),keys=keys1notmatched)
                    raise RuntimeError('ERROR: not all entries in column %s of t1 got matched to column %s of t2' % (matchcol1,matchcol2))
                print(('ERROR: not all entries in column %s of t1 got matched to column %s of t2' % (matchcol1,matchcol2)))
                t1.printtxttable(cols=[matchcol1].extend(cols1),keys=keys1notmatched)
                if returnunmatchedkeys:
                    return(7,keys1notmatched)
                return(7)

        if returnunmatchedkeys:
            return(0,keys1notmatched)
        return(0)
