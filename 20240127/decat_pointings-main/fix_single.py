#!/usr/bin/env python
import sys
import os
import numpy as np
import get_silicon
import astropy.table as at


if __name__=='__main__':
    # TODO - argparse...
    try:
        snid, ra, dec, yseid = sys.argv[1:]
    except Exception as e:
        usage = f'{sys.argv[0]} SNID RA Dec YSE_field'
        print(usage)
        sys.exit(-1)

    # check if we already have an entry for this object
    # if yes, just print it and end
    snid = str(snid)
    yse = at.Table.read('yse/fieldmaps.txt', format='ascii')
    if snid in yse['SNID']:
        ind = yse['SNID'] == snid
        print(yse[ind][0])
        sys.exit(0)

    ra  = [float(ra), ]
    dec = [float(dec),]

    # do not use the center two CCDs
    ccds = list(sorted(set(np.arange(1, 63, 1)) - set([28, 35])))
    thisccd = np.random.choice(ccds, size=1)

    field_ra, field_dec = get_silicon.get_field_center_for_target_on_specific_ccd(ra, dec, thisccd)
    rows = [(yseid, field_ra[0], field_dec[0], ra[0], dec[0], snid, thisccd[0], 1),]
    names=('YSEID', 'RA', 'Dec',  'RACand', 'DecCand', 'SNID', 'candCCD', 'priority')
    dtype = ('S50', 'f8', 'f8', 'f8', 'f8', 'S10', 'i4', 'i4')
    out = at.Table(rows=rows, names=names, dtype=dtype)
    formats={'RACand':'%.6f', 'DecCand':'%.6f', 'RA':'%.6f', 'Dec':'%.6f' }
    for col, form in formats.items():
        out[col].format = form
    print(out)


