#file to analyse SWH differences on global scale for altimetry and ERA5 on a grid and interpolated on altimetry tracks
import os
import geopandas as gpd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import datetime
import pickle
import time
from common_settings import lonmin, lonmax, latmin, latmax
from common_settings import satellites, sat_flags
from common_settings import months, years
from poly_funcs import sort_by_pass, filter_by_flags_seastate, geo_plots, colorbars , mean_swh
from model_functions import modeled_track_new, filter_by_flags_seastate_mt, model_filter_flags, model_filter_reef_buff
#%%
#set reqiered map to True
altimetry_map = False
ERA5_grid_map = True
ERA5_tracks_map = False

#reef data
reef_folder="D:/thesis/reefs"
reef_file=os.listdir(reef_folder)
reefs=[]
for file in reef_file:
    full_name=os.path.join(reef_folder, file)
    reefs.append(gpd.read_file(full_name))
reefs = [reef.set_crs(crs=4326) for reef in reefs]    
reefs_meters = [reef.to_crs(crs=3857) for reef in reefs]

#%%
start_time=time.time()

differences_all = []
att_all = [] 
differences_all_mt = []
att_all_mt = []
differences_all_mt_energy = []
att_all_mt_energy = []
modelled_pass_dict = {}
if altimetry_map == True:
    fig_reef, ax_reef, map_reef = geo_plots(reefs, lonmin, lonmax, latmin, latmax) #map of SWH on reefs
    fig_buf, ax_buf, map_buf = geo_plots(reefs, lonmin, lonmax, latmin, latmax) #map of SWH before reefs
    fig_diff, ax_diff, map_diff = geo_plots(reefs, lonmin, lonmax, latmin, latmax) #map of SWH differences
    
if ERA5_tracks_map == True:
    fig_reef_mt, ax_reef_mt, map_reef_mt = geo_plots(reefs, lonmin, lonmax, latmin, latmax) #map of SWH on reefs
    fig_buf_mt, ax_buf_mt, map_buf_mt = geo_plots(reefs, lonmin, lonmax, latmin, latmax) #map of SWH before reefs
    fig_diff_mt, ax_diff_mt, map_diff_mt = geo_plots(reefs, lonmin, lonmax, latmin, latmax) #map of SWH differences

    fig_reef_energy, ax_reef_energy, map_reef_energy = geo_plots(reefs, lonmin, lonmax, latmin, latmax) #map of SWH on reefs
    fig_buf_energy, ax_buf_energy, map_buf_energy = geo_plots(reefs, lonmin, lonmax, latmin, latmax) #map of SWH before reefs
    fig_diff_energy, ax_diff_energy, map_diff_energy = geo_plots(reefs, lonmin, lonmax, latmin, latmax) #map of SWH differences

if ERA5_grid_map == True:
    fig_reef_grid, ax_reef_grid, map_reef_grid = geo_plots(reefs, lonmin, lonmax, latmin, latmax) #map of SWH on reefs
    fig_buf_grid, ax_buf_grid, map_buf_grid = geo_plots(reefs, lonmin, lonmax, latmin, latmax) #map of SWH before reefs
    fig_diff_grid, ax_diff_grid, map_diff_grid = geo_plots(reefs, lonmin, lonmax, latmin, latmax) #map of SWH differences

#altimetry and ERA5 interpolated on satellite tracks
for year in years:
    for month in months:
        for i in range(len(satellites)):
            print('Satellite id =',i)
            month_folder = 'D:/thesis/ceda_data/{}/{}'.format(year, month)
            model = 'D:/thesis/model_data/model_{}_{}.nc'.format(year,month)#it is better to separate model data by moth for faster interpolation
            model_reference_time = datetime.datetime(int(year), int(month), 1, 0, 0, 0)
            #sort satellite data 
            pass_nums, pass_data_dict = sort_by_pass(month_folder, satellites[i], lonmin, lonmax, latmin, latmax)
            
            if altimetry_map == True:
                buffer_filtered, reef_filtered, difference_swh = filter_by_flags_seastate(pass_data_dict, sat_flags, pass_nums, i, ssl=1, ssh=100)
                #calculate mean per location for plots
                mean_reef = mean_swh(reef_filtered)
                mean_buffer = mean_swh(buffer_filtered)
                mean_differeces  = mean_swh(difference_swh)
                for num in reef_filtered:
                    lons_reef = mean_reef[num]['lon'].tolist() 
                    lats_reef = mean_reef[num]['lat'].tolist() 
                    swhs_reef = mean_reef[num]['swh'].tolist()
                
                    lons_buffer = mean_buffer[num]['lon'].tolist() 
                    lats_buffer = mean_buffer[num]['lat'].tolist() 
                    swhs_buffer = mean_buffer[num]['swh'].tolist()

                    lons_diff = mean_differeces[num]['lon'].tolist()
                    lats_diff = mean_differeces[num]['lat'].tolist()
                    swh_diff_mean = mean_differeces[num]['swh'].tolist()

                    swhs_diff = difference_swh[num]['swh'].tolist()
                    att_percent = difference_swh[num]['attenuation_percent'].tolist()

                    differences_all.extend(swhs_diff)
                    att_all.extend(att_percent)
                    #plot all data on maps
                    map_reef.scatter(lons_reef, lats_reef, s=30, c=swhs_reef, cmap=plt.cm.jet, vmin=0, vmax=3, ax=ax_reef)
                    map_buf.scatter(lons_buffer, lats_buffer, s=30, c=swhs_buffer, cmap=plt.cm.jet, vmin=0, vmax=3, ax=ax_buf)
                    map_diff.scatter(lons_diff, lats_diff, s=30, c=swh_diff_mean, cmap=plt.cm.bwr, vmin=-1, vmax=1, ax=ax_diff)
            
            #model era5 on satelite tracks
            if ERA5_tracks_map == True:    
                for num in pass_data_dict:
                    modeled_track_df, rmse_reef, corr_reef = modeled_track_new(model, model_reference_time, pass_data_dict[num])
                    modelled_pass_dict[num]=modeled_track_df
                print('interpolation done')
                buffer_filtered_mt, reef_filtered_mt, difference_swh_mt = filter_by_flags_seastate_mt(modelled_pass_dict, sat_flags, pass_nums, i, ssl=1, ssh=100)

                for num in reef_filtered_mt:
                    #plot model data on sat tracks
                    lons_reef_mt = reef_filtered_mt[num]['lon'].tolist() 
                    lats_reef_mt = reef_filtered_mt[num]['lat'].tolist() 
                    swhs_reef_mt = reef_filtered_mt[num]['swh'].tolist()
                
                    lons_buffer_mt = buffer_filtered_mt[num]['lon'].tolist() 
                    lats_buffer_mt = buffer_filtered_mt[num]['lat'].tolist() 
                    swhs_buffer_mt = buffer_filtered_mt[num]['swh'].tolist()

                    lons_diff_mt = difference_swh_mt[num]['lon'].tolist()
                    lats_diff_mt = difference_swh_mt[num]['lat'].tolist()
                    swhs_diff_mt = difference_swh_mt[num]['diff_swh'].tolist()
                    att_percent_mt = difference_swh_mt[num]['att_swh'].tolist()
                    energy_diff_mt = difference_swh_mt[num]['diff_energy'].tolist()
                    energy_att_percent_mt = difference_swh_mt[num]['att_energy'].tolist()

                    differences_all_mt.extend(swhs_diff_mt)
                    att_all_mt.extend(att_percent_mt)
                    #plot all data on maps
                    map_reef_mt.scatter(lons_reef_mt, lats_reef_mt, s=30, c=swhs_reef_mt, cmap=plt.cm.jet, vmin=0, vmax=4, ax=ax_reef_mt)
                    map_buf_mt.scatter(lons_buffer_mt, lats_buffer_mt, s=30, c=swhs_buffer_mt, cmap=plt.cm.jet, vmin=0, vmax=4, ax=ax_buf_mt)
                    map_diff_mt.scatter(lons_diff_mt, lats_diff_mt, s=30, c=swhs_diff_mt, cmap=plt.cm.bwr, vmin=-1, vmax=1, ax=ax_diff_mt)
                            
#colorbars
if altimetry_map==True:
    clb, clb_buf, clb_diff = colorbars(fig_reef, fig_buf, fig_diff, ax_reef, ax_buf, ax_diff, hmin = 0, hmax = 3)
if ERA5_tracks_map==True:
    clb, clb_buf, clb_diff = colorbars(fig_reef_mt, fig_buf_mt, fig_diff_mt, ax_reef_mt, ax_buf_mt, ax_diff_mt, hmin = 0, hmax = 4)

#ERA5 on a grid
if ERA5_grid_map == True:
    differences_model = []
    attenuation_model = []
    reef_indx, buf_indx, lon_flat,lat_flat, new_time, swh_flat = model_filter_flags(model,reefs_meters)
    reef_filtered_swh, reef_filtered_swh, difference_swh =model_filter_reef_buff(reef_indx, buf_indx, lon_flat, lat_flat,0.5, 100, swh_flat)
    
    #mean per point for maps
    mean_reef = mean_swh(reef_filtered_swh)
    mean_buffer = mean_swh(reef_filtered_swh)
    mean_differences = mean_swh(difference_swh)
    for time_ind in difference_swh:
        diff_model = difference_swh[time_ind]['swh']
        att_model = difference_swh[time_ind]['attenuation_percent']
        differences_model.extend(diff_model)
        attenuation_model.extend(att_model)
        
    reef_swh_model = mean_reef['swh']
    reef_lon_model = mean_reef['lon']
    reef_lat_model = mean_reef['lat']

    buf_swh_model = mean_buffer['swh']
    buf_lon_model = mean_buffer['lon']
    buf_lat_model = mean_buffer['lat']

    diff_swh_model_mean = mean_differences['swh']
    diff_lon_model = mean_differences['lon']
    diff_lat_model = mean_differences['lat']

    #plot all data on maps
    map_reef_grid.scatter(reef_lon_model, reef_lat_model, s=30, c=reef_swh_model, cmap=plt.cm.jet, vmin=0, vmax=4, ax=ax_reef_grid)
    map_buf_grid.scatter(buf_lon_model, buf_lat_model, s=30, c=buf_swh_model, cmap=plt.cm.jet, vmin=0, vmax=4, ax=ax_buf_grid)
    map_diff_grid.scatter(diff_lon_model, diff_lat_model, s=30, c=diff_swh_model_mean, cmap=plt.cm.bwr, vmin=-1, vmax=1, ax=ax_diff_grid)
                       

    clb, clb_buf, clb_diff = colorbars(fig_reef_grid, fig_buf_grid, fig_diff_grid, ax_reef_grid, ax_buf_grid, ax_diff_grid, hmin = 0, hmax = 4)


end_time=time.time()
print("Execution time = ", (end_time-start_time)/60 , 'minutes')
#%%  
#histogram of swh attenuations
plt.hist(differences_all, bins=100)
plt.hist(differences_all_mt, bins=100)
plt.xlabel('difference [m]')
plt.ylabel('Frequency')
plt.title('Histogram of SWH differences ')
plt.legend(['Altimetry', 'ERA5 on tracks'])
plt.grid(True)
plt.xlim(right=2.5) 
plt.xlim(left=-2.5)
plt.show()

#statistics
#count=0
#positive=[]
differences_all = np.array(differences_all)
pos_alt = [y>=0 for y in differences_all] # why so many nans
pos_array =  differences_all[pos_alt]
pos_count = sum(pos_alt)
print(pos_count)
#for n,i in enumerate(differences_all):
#    if i>=0:
#        count=count+1
#        positive.append(att_all[n])
notnan_diff=np.count_nonzero(~np.isnan(differences_all))
percents=(100/notnan_diff)*pos_count
std_alt=np.nanstd(differences_all)
mean_att_alt=np.nanmean(differences_all)

print("Positive differences ", percents, "%")
print('mean att', mean_att_alt, 'm')
print("Standard deviation of SWH attenuation ", std_alt, "m")

# %%
