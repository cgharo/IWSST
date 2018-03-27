import matplotlib.pyplot as plt
import matplotlib.colors as color
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import xarray as xr
#from cmocean import cm
import os
import fnmatch as fn
import pandas as pd
from glob import iglob
import matplotlib.dates as mdates
import datetime
import locale

locale.setlocale(locale.LC_ALL,'en_US.utf8')

datadir = '/home/cercache/users/cgonzale/data/test'#root path for data
graphdir = '/home/cercache/users/cgonzale/snapshots/'#root path for snapshots of SST


#Parameters that will change to define region and period of interest
Region = 'Mozamb4'

region_folder = datadir + os.sep + Region + os.sep

#Find available files in the region
rootdir_glob = region_folder+ '**/*'
file_list = [f for f in iglob(rootdir_glob, recursive=True) if os.path.isfile(f)]

#check how many available  sensors for the given region
fsensor = [f.rsplit(os.sep,4)[1] for f in file_list]

sensor = list(set(fsensor))

# Read date from filenames and convert the string into datetime
date = [f.rsplit(os.sep,1)[1][0:14] for f in file_list if 'modis' not in f]
#Take into account different file name for modis
date_modis = [f.rsplit(os.sep,1)[1][0:8] for f in file_list if 'modis'  in f ]
#add 000000 to account for HHMMSS
date_modis_h = [d+'000000' for d in date_modis]
#concatenate both lists
date_merge = date + date_modis_h



format='%Y%m%d%H%M%S'
date_time = pd.to_datetime(date_merge, format=format)




#Find dates with more than one SST image

data_dict = {'date' : date_time, 'sensor':fsensor,'path':file_list}
df = pd.DataFrame(data_dict)

df = df.sort_values('date')

#Append new column for day_index
df['day_index'] = df['date'].apply(lambda x:str(x.year)+str(x.month)+str(x.day))
df = df[df.day_index.duplicated(keep=False)]
a = df.groupby('day_index').apply(lambda x: list(x.path))
list_day_images = a.values


#print results per date and path to data file
data_file_fname = pd.DataFrame([])

for days in range(len(list_day_images)):
    date = [list_day_images[days][0].rsplit(os.sep,1)[1][0:8]]
    print('###############')
    print(date)
    print('    files:')
    for f in list_day_images[days]:
        data_file_fname = data_file_fname.append(pd.DataFrame({'date': date,'path':f}, index=[0]))
        print(f)


#save csv file
filename = datadir+ os.sep + 'Available_images_LOG' +os.sep + Region+'_available_SST_images.csv'
data_file_fname.to_csv(filename, sep='\t')


#Figure with available data
#inizialitation of variable to plot
plot_sensor = np.zeros([len(date_time)])

years = mdates.YearLocator()   # every year
months = mdates.MonthLocator()  # every month
days = mdates.DayLocator()  # every day
yearsFmt = mdates.DateFormatter('%Y')

symb_sensor = ['o','^','*']


locale.setlocale(locale.LC_ALL,'en_US.utf8')
fig, ax = plt.subplots()

for s in range(len(sensor)):
    isensor = [True if sensor[s] in i else False for i in file_list]
    ax.plot(date_time[isensor], plot_sensor[isensor]+s+0.5,symb_sensor[s],label = sensor[s])


datemin = datetime.date(date_time.min().year, 1, 1)
datemax = datetime.date(date_time.date.max().year + 1, 1, 1)
ax.set_xlim(datemin, datemax)
ax.set_ylim(0, len(sensor))

ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
ax.xaxis.set_major_locator(months)
ax.xaxis.set_minor_locator(days)
ax.xaxis.set_ticks_position('bottom')


ax.legend()
ax2 = ax.twiny()

#offset the twin axes below the host
ax2.spines["bottom"].set_position(("axes", -0.2))

#turn on the frame for the twin axis, but then hide all
ax2.set_frame_on(True)
ax2.patch.set_visible(False)
#for sp in ax2.spines.intervalues():
#    sp.set_visible(False)
ax2.spines["bottom"].set_visible(True)

ax2.xaxis.set_ticks_position('bottom')
ax2.xaxis.set_label_position('bottom')

# format the ticks
ax2.xaxis.set_major_locator(years)
ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax2.set_xlabel('Time')
#ax2.xaxis.set_minor_locator(months)

# rotates and right aligns the x labels, and moves the bottom of the
# axes up to make room for them
fig.autofmt_xdate()

plt.show(block = False)




