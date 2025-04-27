#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 23 13:04:40 2025

@author: arest
"""

import sys, os, re
from obsplan_baseclass import obsplan_baseclass
import numpy as np



class modify_obsplanclass(obsplan_baseclass):
    def __init__(self):
        obsplan_baseclass.__init__(self)
    
    def define_optional_arguments(self,parser=None,usage=None,conflict_handler='resolve'):
        parser = obsplan_baseclass.define_optional_arguments(self,parser=parser,usage=usage,conflict_handler=conflict_handler)
        parser.add_argument('-a','--add', type=str, nargs="+",  default=None, help='list of jsonID\'s to add to the obsplan')
        parser.add_argument('-r','--remove', type=int, nargs="+",  default=None, help='list of jsonID\'s to remove to the obsplan')
        return(parser)
    
    def add_targets(self,addlist):
        if (addlist is None) or (addlist is []):
            if self.verbose>1: print('No targets to add....')
            return(0)
        
        ixs = []
        for s in addlist:
            if re.search('^\d+$',s) is not None:
                jsonID = int(s)
                ix_found = self.jsontable.ix_equal('jsonID',jsonID)
                if len(ix_found)==0:
                    raise RuntimeError(f'Could not find jsonID {jsonID}')
                elif len(ix_found)>1:
                    self.jsontable.write(indices=ix_found)
                    raise RuntimeError(f'More than one entry for jsonID {jsonID}')
            else:
                ix_found = self.jsontable.ix_equal('json_short',s)
                if len(ix_found)==0:
                    raise RuntimeError(f'Could not find json_short={s}!')
                elif len(ix_found)>1:
                    self.jsontable.write(indices=ix_found)
                    raise RuntimeError(f'More than one entry for json_short {s}')
            ixs.append(ix_found[0])
            
        if self.verbose>2:
            print(f'Adding {len(ixs)} json files to the obsplan:')
            self.jsontable.write(indices=ixs)

        # make sure they are not anymore in the ordered list
        self.jsontable.t.loc[ixs,'order']=-1

        # now at the targets to the ordered list
        self.add_unordered_targets(ixs)

        # recalculate the slew times
        self.calc_slew()
        
    def remove_targets(self,removelist):
        if (removelist is None) or (removelist is []):
            if self.verbose>1: print('No targets to remove....')
            return(0)
        
        ixs = []
        for jsonID in removelist:
            print(f'jsonID: {jsonID}')
            ix = self.jsontable.ix_equal('jsonID',jsonID)
            if len(ix)==0:
                raise RuntimeError(f'Could not find jsonID {jsonID}')
            elif len(ix)>1:
                self.jsontable.write(indices=ix)
                raise RuntimeError(f'More than one entry for jsonID {jsonID}')
            
            ixs.append(ix[0])
        
        # set the order to -1, and remove all the times
        self.jsontable.t.loc[ixs,'order']=-1
        self.jsontable.t.loc[ixs,'UT']=np.nan
        self.jsontable.t.loc[ixs,'CLT']=np.nan
        self.jsontable.t.loc[ixs,'t_slew[m]']=np.nan
        self.jsontable.t.loc[ixs,'t_tot[m]']=np.nan
        self.jsontable.write(indices=ixs)
        
        # recalculate the slew times
        self.calc_slew()
        
        return(0)
        
    
if __name__ == '__main__':
    


    obsplan = modify_obsplanclass()
    parser = obsplan.define_optional_arguments()
    args = parser.parse_args()

    # the arguments are saved in query.params
    obsplan.get_arguments(args)
    obsplan.initialize4date(obsplan.params['date'])
    
    #obsplan.init_airmass()
    #sys.exit(0)
    obsplan.load_airmass()
    
    if not os.path.isfile(obsplan.obsplan_filename):
        raise RuntimeError(f'Obsplan {obsplan.obsplan_filename} does not exist!')
        #obsplan.find_new_jsonfiles()
        
    obsplan.loadobsplan()
    obsplan.find_new_jsonfiles()
    obsplan.update_tables()
    
    obsplan.calc_slew()
    
    obsplan.remove_targets(args.remove)
    obsplan.add_targets(args.add)
    
    obsplan.calc_obsUTtime(firsthalf=args.firsthalf,secondhalf=args.secondhalf,add2starttime_minutes=args.add2starttime_minutes)

    obsplan.calc_airmass()
    obsplan.calc_moon_distances()
    
    if args.makeplots or args.showplots:
        obsplan.airmass_plot()
        obsplan.moon_plot()

    obsplan.wrap_up_and_summary()
    obsplan.saveobsplan()
    obsplan.wrap_up_and_summary()

    if args.showplots:
        obsplan.show_plots()
