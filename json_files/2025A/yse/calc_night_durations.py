#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 12:31:14 2024

@author: arest
"""
from astroplan import Observer,moon
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import (EarthLocation, SkyCoord, AltAz, get_sun,
                                 get_moon, Angle, Longitude)

if __name__=='__main__':
    YYYYMMDDs = [
        '20250215',
        '20250315',
        '20250415',
        '20250515',
        '20250615',
        '20250715',
        ]
    ctio = Observer.at_site("CTIO",timezone='America/Santiago')

    horizon = 14
    for YYYYMMDD in YYYYMMDDs:
        date = f'{YYYYMMDD[:4]}-{YYYYMMDD[4:6]}-{YYYYMMDD[6:8]}T06:00:00'
        t = Time(date)
        # add 0.5 days since t is in UT. Otherwise the values are not correct
        twi =  ctio.tonight(t+0.5,horizon=-horizon*u.deg)
        print('UT %d deg twilight: %s %s' % (-horizon,twi[0].to_value('isot'),twi[1].to_value('isot')))
        dt = twi[1] - twi[0]
        print(f'night duration: {dt.to_value("h"):.2f} hours')