import numpy as np
import pandas as pd
import gspread
from oauth2client.service_account import ServiceAccountCredentials
import readobslogs as ro
import makeobservabilityplot as mop
from datetime import datetime
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import matplotlib.pyplot as plt
import os
import editjson as ej

def split(word):
    return [char for char in word]

scope = ['https://spreadsheets.google.com/feeds',
         'https://www.googleapis.com/auth/drive']

creds = ServiceAccountCredentials.from_json_keyfile_name('.gclient.json', scope)
client = gspread.authorize(creds)

sheet = client.open("DEBASS Sample").sheet1

list_of_lists = sheet.get_all_values()

goodcollist = ['snid','Following?','TNS class','RA','DEC']
colnames = list_of_lists[0]
dict = {}
for i,col in enumerate(colnames):
    if col in goodcollist:
        dict[col] = sheet.col_values(i+1)[1:]

df = pd.DataFrame.from_dict(dict)

from datetime import datetime
currdate=datetime.today().strftime('%Y-%m-%d')
datestr = input(f'Please Enter Observing Date in the following format (YYYY-MM-DD)\ndefault {currdate}\n')
if datestr == '': datestr=currdate
obsdict = ro.run(verbose=False)


skiprows = 60
for i,row in df.iterrows():
    if i < skiprows: continue
    if 'FINISHED' in row['Following?']: continue
    if 'ABANDON' in row['Following?']:	continue
    if 'LOST' in row['Following?']:  continue
    if '91' in row['TNS class']:  continue
    if ('ia' in row['TNS class'].lower()) | ('?' in row['TNS class']):
        os.system('clear')
        print(obsdict[row['snid']])
        mop.doplot(datestr,ra=float(row['RA']),dec=float(row['DEC']),name=row['snid'],block=False)
        priority = input('Please enter a priority for this object (ie 1 2 3)\n')
        filters = split(input('Please enter filters for this object (ie griz)\n'))
        exptimes = []
        for f in filters:
            exptime = input(f'Enter Exptime {f} (default 15)\n')
            if exptime == '': exptime = '15'
            exptimes.append(exptime)
        plt.clf()
        if os.path.exists('jsons/2020B-0053_DEBASS_Brout/BASE/%s.json'%row['snid']):
            os.system('rm jsons/2020B-0053_DEBASS_Brout/EVERYTHING/%s_P*.json'%(row['snid']))
            ej.edit('jsons/2020B-0053_DEBASS_Brout/BASE/%s.json'%row['snid'],priority,filters,exptimes,
                 'jsons/2020B-0053_DEBASS_Brout/EVERYTHING/%s_P%s.json'%(row['snid'],priority))
        else:
            print('WARNING: jsons/2020B-0053_DEBASS_Brout/BASE/%s.json'%row['snid']+'\nDoes Not Exist')
            input('press enter to continue')
