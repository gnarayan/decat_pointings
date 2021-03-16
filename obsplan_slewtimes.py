from glob import glob
import numpy as np
import slewtimes

for opdir in glob('jsons/obsplans/*'):
    date = opdir.split('/')[-1]
    jsons = glob(opdir+'/*.json')
    for json in jsons:
        t = slewtimes.time_for_single_json(json)
        print(json,'%s\t%.1f minutes'%(date,t/60.))
    t = slewtimes.total_time_from_jsons(np.sort(jsons))
    print('%s\t%.1f minutes'%(date,t/60.))
