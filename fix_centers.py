#!/usr/bin/env python
import sys
import os
import get_silicon
from astropy.table import Table


if __name__=='__main__':
    night = '20210318'
    prop_file = f'2021A/{night}/decat_YSE_proposed_{night}.txt'
    prop = Table.read(prop_file, format='ascii')

    field_ra, field_dec = get_silicon.get_field_center_for_target_on_specific_ccd(prop['RACand'], prop['DecCand'], prop['CandCCD'])
    prop.rename_column('RACand', 'Cand_RA')
    prop.rename_column('DecCand', 'Cand_Dec')
    formats={'Cand_RA':'%.6f', 'Cand_Dec':'%+.6f', 'Field_RA':'%.6f', 'Field_Dec':'%+.6f' }
    prop['Field_RA'] = field_ra
    prop['Field_Dec'] = field_dec
    out_file = f'2021A/{night}/decat_YSE_list_{night}.txt'
    include_names=('ID', 'Field_RA', 'Field_Dec', 'Cand_RA', 'Cand_Dec', 'SNID', 'CandCCD')
    out_cols = [prop[x] for x in include_names]
    out = Table(out_cols, names=include_names)
    out.write(out_file, format='ascii.fixed_width',\
            formats=formats, names=include_names, delimiter='   ',\
            overwrite=True)


