import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import xarray as xr
from netCDF4 import Dataset
from cmocean import cm
import pandas as pd
import os
import fnmatch as fn
from glob import iglob
import sys
import argparse
from utils import *
from palette import *
palette = '/home3/homedir7/perso/cgonzale/IWAVE/script/palette/medspiration.rgb'
csst = getColorMap( rgbFile = palette )

######################
def date_to_nth_day(date, format='%Y%m%d'):
    date = pd.to_datetime(date, format=format)
    new_year_day = pd.Timestamp(year=date.year, month=1, day=1)
    return (date - new_year_day).days + 1

def safe_make_folder(i):
    '''Makes a folder if not present'''
    try:  
        os.makedirs(i)
        print('creating: ', i)
    except:
        pass

########################


def plot_snapshots_SST(Region , datadir = '/home/cercache/users/cgonzale/data/test' ,graphdir=  '/home/cercache/users/cgonzale/snapshots'):
  


    region_folder = datadir+os.sep + Region+os.sep
    rootdir_glob = region_folder+ '**/*'
    #find the path of available data
    file_list = [f for f in iglob(rootdir_glob, recursive=True) if os.path.isfile(f)]


    for file in file_list:  
        path_graph =  file.rsplit(os.sep,1)[0].split(os.sep,7)[-1]
        #Create path if it does not exist
        safe_make_folder(graphdir + os.sep+ path_graph)
    
        #may check if the figure already exists before open file and do the plot
        #ds = xr.open_dataset(file)
        #print(file)
        ncfile=Dataset(file,'r')
        lat = ncfile.variables['lat'][:]
        lon = ncfile.variables['lon'][:]
        sst = ncfile.variables['sea_surface_temperature'][:]
        ncfile.close()

        sst = np.squeeze(sst)
        lon2d , lat2d = np.meshgrid(lon,lat)

        fig1=plt.figure(figsize=(9, 5))
        ax = plt.axes(projection=ccrs.PlateCarree())
        #toplt = ds['sea_surface_temperature']
        #hdl = ax.imshow(toplt, origin='lower', transform=ccrs.PlateCarree(), cmap =csst)
        hdl = ax.pcolormesh(lon2d,lat2d,sst, \
	                transform = ccrs.PlateCarree(),cmap=csst)
        ax.coastlines(resolution='110m', color ='k' )
        ax.add_feature(cfeature.LAND, facecolor = '0.75')
        ax.gridlines(draw_labels = True)
        cb = plt.colorbar(hdl,ax=ax)
        cb.set_label('SST [K]')

        name_fig = file.rsplit(os.sep,1)[1].rsplit('.',2)[0]
        plt.savefig(graphdir+os.sep+path_graph+os.sep+name_fig+'.png')
        plt.show(block = False)
        plt.close(fig1)
        print('saving:',graphdir+os.sep+path_graph+os.sep+name_fig+'.png')
        plt.close(fig1)

def main():
    #Region = sys.argv[1]
    parser = argparse.ArgumentParser()
    parser.add_argument('Region', help = 'Name of the region to inspec', type = str)
    parser.add_argument('-o', '--output', help = 'Root directory were plots will be saved')
    parser.add_argument('-i', '--input', help = 'Root directory were NetCDF are stored')

    args = parser.parse_args()
    print('arg.out:', args.output)

    if args.output:
        graphdir = args.output
    else:
        graphdir= '/home/cercache/users/cgonzale/snapshots'

    if args.input:
        datadir = args.input
    else:
        datadir= '/home/cercache/users/cgonzale/data/test'

    
    plot_snapshots_SST(args.Region, datadir = datadir , graphdir =  graphdir)

if __name__ == '__main__':
   main()


#python3 plot_snapshots_SST.py Mozamb1 -o /home/cercache/users/cgonzale/graph_dir -i /home/cercache/users/cgonzale/datain


