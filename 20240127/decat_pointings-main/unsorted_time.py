from glob import glob
import numpy as np
import slewtimes
import sys
folder = sys.argv[1]

for opdir in glob(folder):
    date = opdir.split('/')[-1]
    jsons = glob(opdir+'/*.json')
    t = slewtimes.total_time_from_jsons(np.sort(jsons),sort=False)
    print('%s\t%.1f minutes'%(date,t/60.))
