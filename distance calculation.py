import xarray as xr
import pandas as pd
import os
import geopandas as gpd
from shapely.geometry import Point
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import interp1d
import pickle
from poly_funcs import geo_plots
#%%
'''file to calculate location flags'''

#define satellite id number and path to the folder with data for one month
satellite_number = 4 #0-cryosat, 1-jason1, 2-jason2, 3-jason3, 4-saral, 5- sentinel3a, 6-envisat
month_folder='D:/thesis/ceda_data/2014/01'
month=month_folder[25:27]
year=month_folder[20:24]

reef_folder="D:/thesis/reefs"
reef_file=os.listdir(reef_folder)
reefs=[]
for file in reef_file:
    full_name=os.path.join(reef_folder, file)
    reefs.append(gpd.read_file(full_name))
reefs = [reef.set_crs(crs=4326) for reef in reefs]    
reefs_meters = [reef.to_crs(crs=3857) for reef in reefs]

# %%
def flags_all_reefs(month_folder,satellite_number):
    """Function to create dictinary with flags for all reefs"""
 
    all_files=[]
    pass_set = set()
    days=os.listdir(month_folder)
    days = days
    for day in days:
        file='D:/thesis/ceda_data/{}/{}/{}/ESACCI-SEASTATE-L3-SWH-MULTI_1D-{}{}{}-fv01.nc'.format(year, month, day, year, month, day)
        all_files.append(file)
    month_02='D:/thesis/ceda_data/2007/02'
    # add the folowing part of code if satellite repeat cycle is more than 30 days
    days=os.listdir(month_02)
    days=days[0:5]
    for day in days:
        file='D:/thesis/ceda_data/2007/02/{}/ESACCI-SEASTATE-L3-SWH-MULTI_1D-200702{}-fv01.nc'.format(day, day)
        all_files.append(file)

    for file in all_files:
        ds=xr.open_dataset(file)
        df=ds.to_dataframe()
        pass_data=df[df['satellite']== satellite_number]
        #making a list of pass numbers
        pass_set.update(pass_data['relative_pass_number'].unique())
    pass_nums = sorted(list(pass_set))

    pass_flag_dict = {}
    for file in all_files:
        ds=xr.open_dataset(file)
        df=ds.to_dataframe()
        sat_data=df[df['satellite']== satellite_number]
        for num in pass_nums:
            pass_data = sat_data[sat_data['relative_pass_number'] == num]
            lon = pass_data['lon']
            lat= pass_data['lat']
            dist_2_coast=pass_data['distance_to_coast']
            #remove irrelevant data
            close_indices = (dist_2_coast < 500000) & (lat < 30) & (lat > -30)
            lon_close = lon[close_indices]
            lat_close = lat[close_indices]
            dist_close=dist_2_coast[close_indices]
            if not len(lon_close) ==0:
                points = gpd.GeoSeries([Point(lon, lat) for lon, lat in zip(lon_close, lat_close)])
                points = points.set_crs(crs=4326)
                points_meters=points.to_crs(crs=3857)

                flag=[]
                #0 - empty zone; 1 - reef; 2 - buffer zone
                for point in points_meters:
                    point_flag = 0
                    for reef_multipolygon in reefs_meters:
                        one_reef_meters=reef_multipolygon.geometry.iloc[0]
                        #one_reef_meters_list=[reef for reef in one_reef_meters]
                        #for reef in one_reef_meters: #one_reef_meters_list:
                        distance = one_reef_meters.distance(point)# this doednt really work bc it calculates one distance to a multipolygon
                        if distance == 0:
                                point_flag = 1
                                break
                        elif 3500 <= distance <= 11000:
                                point_flag = 2
                                break
                    flag.append(point_flag)
                if num not in pass_flag_dict:
                    pass_flag_dict[num] = pd.DataFrame({'lon': lon_close,
                                                         'lat': lat_close, 'flag': flag})
                else: 
                    pass_flag_dict[num] = pd.concat([pass_flag_dict[num], pd.DataFrame({'lon': lon_close,
                                                                                        'lat': lat_close,
                                                                                        'flag': flag})])

           
    return pass_flag_dict

#%% 
#calculate flags
pass_flag_dict = flags_all_reefs(month_folder, satellite_number)
#check flags
fig, ax_c, map_c = geo_plots(reefs,120, 160, -30, 0)
for num in pass_flag_dict:
    lons=pass_flag_dict[num]['lon']
    lats=pass_flag_dict[num]['lat']
    flages=pass_flag_dict[num]['flag']    
    map_c.scatter(lons, lats, s=40, c=flages, cmap=plt.cm.jet, vmin=0, vmax=4, ax=ax_c)

#save flags
with open('sat_{}_flags_final.pickle'.format(satellite_number), 'wb') as handle:
        pickle.dump(pass_flag_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
# %%
