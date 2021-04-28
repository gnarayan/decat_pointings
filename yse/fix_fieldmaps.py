import astropy.table as at
f = at.Table.read('fieldmaps.txt', format='ascii')
f['Dec'].format='%.6f'
f['DecCand'].format='%.6f'
f['RACand'].format='%.6f'
f['RA'].format='%.6f'
f.sort(['RA', 'Dec', 'SNID'])
f.write('new_fieldmaps.txt', format='ascii.fixed_width', delimiter='  ', overwrite=True)
