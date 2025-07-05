#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 23 13:04:40 2025

@author: arest
"""

import sys, os
from obsplan_baseclass import obsplan_baseclass

class new_obsplanclass(obsplan_baseclass):
    def __init__(self):
        obsplan_baseclass.__init__(self)
    


if __name__ == '__main__':

    obsplan = new_obsplanclass()
    parser = obsplan.define_optional_arguments()
    args = parser.parse_args()

    # the arguments are saved in query.params
    obsplan.get_arguments(args)
    obsplan.initialize4date(obsplan.params['date'])
    
    #obsplan.init_airmass()
    obsplan.load_airmass()
    #if os.path.isfile(obsplan.obsplan_filename):
    #    raise RuntimeError(f'Obsplan {obsplan.obsplan_filename} already exists!')
    
    obsplan.find_new_jsonfiles()
    obsplan.update_tables()
    
    (ixs_ordered,ixs_notordered) = obsplan.generic_target_ordering()

    obsplan.calc_slew(ixs_ordered=ixs_ordered)
    obsplan.add_unordered_targets(ixs_notordered,ixs_ordered=ixs_ordered)
    obsplan.calc_obsUTtime(firsthalf=args.firsthalf,secondhalf=args.secondhalf,add2starttime_minutes=args.add2starttime_minutes)

    obsplan.calc_airmass()
    obsplan.calc_moon_distances()
    if args.makeplots or args.showplots:
        obsplan.airmass_plot()
        obsplan.moon_plot()

    obsplan.saveobsplan()
    obsplan.wrap_up_and_summary()

    if args.showplots:
        obsplan.show_plots()
