from glob import glob
import numpy as np
import slewtimes
import sys
folder = sys.argv[1]

for opdir in glob(folder):
    date = opdir.split('/')[-1]
    jsons = glob(opdir+'/*.json')
    t,t_p1 = slewtimes.total_time_from_jsons(np.sort(jsons))
    print('%s\t%.1f minutes'%(date,t/60.))
    print('%s_P1\t%.1f minutes'%(date,t_p1/60.))
