import pandas as pd
import os
import geopandas as gpd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pickle
import time
from netCDF4 import Dataset
import datetime
from shapely.geometry import Point
from scipy.spatial.distance import cdist
from scipy.stats import pearsonr
import math
from scipy.interpolate import griddata
from scipy.interpolate import RegularGridInterpolator
import numpy.ma as ma

from poly_funcs import geo_plots


reef_folder="D:/thesis/reefs"
reef_file=os.listdir(reef_folder)
reefs=[]
for file in reef_file:
    full_name=os.path.join(reef_folder, file)
    reefs.append(gpd.read_file(full_name))
reefs = [reef.set_crs(crs=4326) for reef in reefs]    
reefs_meters = [reef.to_crs(crs=3857) for reef in reefs]
#%%
#model='D:/thesis/model_data/adaptor.mars.internal-1687950556.6451876-18738-2-a0e6ccf0-3abd-40b4-8f2b-77ec43ed5713.nc'
#reference_time = datetime.datetime(2009, 3, 1, 0, 0, 0)
def mean_per_point(model, reference_time, start_day, number_of_days):
    '''function to calculate mean'''
    fid=Dataset(model)
    lon=fid.variables['longitude'][:]
    lat=fid.variables['latitude'][:]
    swh=fid.variables['swh'][:]
    mwp=fid.variables['mwp'][:]
    time_model=fid.variables['time'][:]
    new_time=[]
    for i in range(len(time_model)):
        time_delta = datetime.timedelta(hours=i)
        result_time = reference_time + time_delta 
        new_time.append(result_time)
    print(new_time)
    for n in range(len(new_time)):
        
        if new_time[n].day == start_day.day:
            start_index=n
            end_index = n+number_of_days*24
            break
    # Extract the relevant portion of swh based on the desired indices
    swh_slice = swh[start_index:end_index]
    mwp_slice = mwp[start_index:end_index]
    lon_mesh, lat_mesh = np.meshgrid(lon, lat)
    # Flatten the arrays 
    lon_flat = lon_mesh.flatten()
    lat_flat = lat_mesh.flatten()

    location_wave_heights={}
    for n in range(len(swh_slice)):
        swh_flat = swh_slice[n].flatten()
        mwp_flat = mwp_slice[n].flatten()
        for i in range(len(lon_flat)):
            if (lon_flat[i],lat_flat[i]) in location_wave_heights:
                                location_wave_heights[(lon_flat[i], lat_flat[i])]['sum_swh'] += swh_flat[i]
                                location_wave_heights[(lon_flat[i], lat_flat[i])]['sum_mwp'] += mwp_flat[i]
                                location_wave_heights[(lon_flat[i], lat_flat[i])]['count'] += 1
            else:
                # If the location hasn't been seen, add it to the dictionary with its wave height
                location_wave_heights[(lon_flat[i], lat_flat[i])] = {
                                    'sum_swh': swh_flat[i], 'sum_mwp': mwp_flat[i],  'count': 1}
                                

        mean_swh=[]
        mean_lon=[]
        mean_lat=[]
        mean_mwp=[]
        # Iterate through each unique location in the pass
        for unique_location, location_data in location_wave_heights.items():
            # Calculate the mean wave height for the location
            mean_wave_height = location_data['sum_swh'] / location_data['count']
            mean_wave_period = location_data['sum_mwp'] / location_data['count']
            mean_lon.append(unique_location[0])
            mean_lat.append(unique_location[1])
            mean_swh.append(mean_wave_height)
            mean_mwp.append(mean_wave_period)
    mean_parameters = pd.DataFrame({'lons': mean_lon, 'lats':mean_lat, 'swh': mean_swh, 'mwp':mean_mwp })
    return mean_parameters

# %% filter out reef points

def filter_reef_buffer(mean_parameters, reefs):
    lon_not_filtered = mean_parameters['lons']
    lat_not_filtered= mean_parameters['lats']
    swh_not_filtered = mean_parameters['swh']
    mwp_not_filtered = mean_parameters['mwp']  

                    
    points = gpd.GeoSeries([Point(lon, lat) for lon, lat in zip(lon_not_filtered, lat_not_filtered)])
    points = points.set_crs(crs=4326)
    points_meters=points.to_crs(crs=3857)

    loc_reef=[]
    swh_reef=[]
    mwp_reef=[]
    loc_buf=[]
    swh_buf=[]
    mwp_buf=[]
    for point in range(len(points)):
        for reef_multipolygon in reefs:
                distance = reef_multipolygon.distance(points_meters[point])
                distance=distance.values
                if distance <=200:
                    loc_reef.append(points[point])
                    swh_reef.append(swh_not_filtered[point])
                    mwp_reef.append(mwp_not_filtered[point])
                    break
                elif 1000 <= distance <= 40000:#justify the range by the grid size of era5
                    loc_buf.append(points[point])
                    swh_buf.append(swh_not_filtered[point])
                    mwp_buf.append(mwp_not_filtered[point])
                    break

    lons_reef = [point.x for point in loc_reef]  
    lats_reef = [point.y for point in loc_reef]

    lons_buf = [point.x for point in loc_buf]  
    lats_buf = [point.y for point in loc_buf]

    reef_param=pd.DataFrame({'lons':lons_reef,'lats':lats_reef, 'swh':swh_reef, 'mwp':mwp_reef})
    buf_param=pd.DataFrame({'lons':lons_buf,'lats':lats_buf, 'swh':swh_buf, 'mwp':mwp_buf})
    return reef_param, buf_param

#%% energy calculations
def energy_flux(param_dict, swh_reef):
    density = 1025  # Water density in kg/m³
    gravity = 9.81  # Acceleration due to gravity in m/s²
    lons=param_dict['lon']
    lats=param_dict['lat']
    swhs=param_dict['swh']
    mwps=param_dict['mwp']
    energy=[]
    for value in range(len(swh_reef)):
        wave_energy =((gravity**2)/(64*np.pi)) * density * swh_reef[value]**2 * mwps[value]
        wave_energy=wave_energy/1000 #to convert from W/m to KW/m
        energy.append(wave_energy)
    
    energy_dict=pd.DataFrame({'lons':lons, 'lats':lats, 'energy_flux': energy})

    return(energy_dict)

#%%tests
#write a function for model one day - one time to create a mask two masks - one for buffer, one for mask. Measurments from every hour multiply w wthis masks. find closes buffer and reef points and substract them 
'''mean_par = mean_per_point(model, reference_time, 2, 3)
reef_param, buffer_param = filter_reef_buffer(mean_par,reefs_meters)
energy_dick = energy_flux(reef_param)
mean_lons=energy_dick['lons']
mean_lats = energy_dick['lats']
mean_swhs = energy_dick['energy_flux']

fig, ax, map_base = geo_plots(reefs)
map_base.scatter(mean_lons, mean_lats, s=30, c=mean_swhs, cmap=plt.cm.viridis, vmin=0, vmax=3, ax=ax)

model='D:/thesis/model_data/adaptor.mars.internal-1687950556.6451876-18738-2-a0e6ccf0-3abd-40b4-8f2b-77ec43ed5713.nc'
fid=Dataset(model)
lon=fid.variables['longitude'][:]
lat=fid.variables['latitude'][:]
swh=fid.variables['swh'][:]
time_model=fid.variables['time'][:]
year = 2009
month = 3
day = 1
reference_time = datetime.datetime(year, month, day, 0, 0, 0)
new_time=[]
for i in range(len(time_model)):
    time_delta = datetime.timedelta(hours=i)
    result_time = reference_time + time_delta
    #print(result_time)
    new_time.append(result_time)
print(new_time[0])'''

# time[0] = 1.03.2009 00:00
# time [25] = 2.032009 00:00
#Australia UTC +10
#data is already cut to my region by copernicus
# make plots of enegy flux for needed dates (date format?)
#%%
def mean_of_flags_seastate_model(lons,lats,swh,flags,sea_state_low, sea_state_high):
    mean_swhs_reef = []
    mean_lat_reef =[]
    mean_lon_reef=[]
    sum_swh_reef = []
    sum_lat_reef = []
    sum_lon_reef = []
    mean_swhs_buf = []
    mean_lat_buf =[]
    mean_lon_buf=[]
    sum_swh_buf = []
    sum_lat_buf = []
    sum_lon_buf = []
    
    counter_buf=0
    counter_reef=0
    is_sequence_open = False
    for i in range(0, (len(flags))):
        if flags[i] == 2 and sea_state_low <= swh[i] <=sea_state_high:
            sum_swh_buf.append(swh[i])
            sum_lon_buf.append(lons[i])
            sum_lat_buf.append(lats[i])
            counter_buf+=1
            is_sequence_open = True
        elif counter_buf>0:
                mean_swhs_buf.append(np.nanmean(sum_swh_buf))
                mean_lon_buf.append(sum_lon_buf[0])
                mean_lat_buf.append(sum_lat_buf[0])
                counter_buf = 0
                sum_swh_buf = []
                sum_lat_buf = []
                sum_lon_buf = []
        if flags[i] == 1: 
            if is_sequence_open == True:
                counter_reef += 1
                sum_swh_reef.append(swh[i])
                sum_lat_reef.append(lats[i])
                sum_lon_reef.append(lons[i])
            
        elif counter_reef > 0:
                mean_swhs_reef.append(np.nanmean(sum_swh_reef[1:(len(sum_swh_reef)-1)]))
                mean_lon_reef.append(np.nanmean(sum_lon_reef))
                mean_lat_reef.append(np.nanmean(sum_lat_reef))
                counter_reef = 0
                sum_swh_reef = []
                sum_lat_reef = []
                sum_lon_reef = []
                is_sequence_open = False
    for i in range(len(flags)-1, -1, -1):
        if flags[i] == 2 and swh[i]>=sea_state_low and swh[i]<=sea_state_high:
            sum_swh_buf.append(swh[i])
            sum_lon_buf.append(lons[i])
            sum_lat_buf.append(lats[i])
            counter_buf+=1
            is_sequence_open = True
        elif counter_buf>0:
            mean_swhs_buf.append(np.nanmean(sum_swh_buf))
            mean_lon_buf.append(sum_lon_buf[0])
            mean_lat_buf.append(sum_lat_buf[0])
            counter_buf = 0
            sum_swh_buf = []
            sum_lat_buf = []
            sum_lon_buf = []
        if flags[i] == 1: 
            if is_sequence_open == True:
                counter_reef += 1
                sum_swh_reef.append(swh[i])
                sum_lat_reef.append(lats[i])
                sum_lon_reef.append(lons[i])
        elif counter_reef > 0:
                mean_swhs_reef.append(np.nanmean(sum_swh_reef[1:(len(sum_swh_reef)-1)]))
                mean_lon_reef.append(np.mean(sum_lon_reef))
                mean_lat_reef.append(np.mean(sum_lat_reef))
                counter_reef = 0
                sum_swh_reef = []
                sum_lat_reef = []
                sum_lon_reef = []
                is_sequence_open = False
    return mean_lon_reef, mean_lat_reef, mean_swhs_reef, mean_lon_buf, mean_lat_buf, mean_swhs_buf
#%%
def model_filter_flags(model,reefs_meters):
    fid=model
    fid=Dataset(model)
    lon=fid.variables['longitude'][:]
    lat=fid.variables['latitude'][:]
    swh=fid.variables['swh'][:]
    time_model=fid.variables['time'][:]
    year = 2017
    month = 1
    day = 1
    reference_time = datetime.datetime(year, month, day, 0, 0, 0)
    new_time=[]
    for i in range(len(time_model)):
        time_delta = datetime.timedelta(hours=i)
        result_time = reference_time + time_delta
        #print(result_time)
        new_time.append(result_time)
    lon_mesh, lat_mesh = np.meshgrid(lon, lat)
    # Flatten the arrays for vectorized scatter plot
    lon_flat = lon_mesh.flatten()
    lat_flat = lat_mesh.flatten()

    points = gpd.GeoSeries([Point(lon, lat) for lon, lat in zip(lon_flat, lat_flat)])
    points = points.set_crs(crs=4326)
    points_meters=points.to_crs(crs=3857)

    flag=[] #0 - empty zone; 1 - reef; 2 - buffer zone
    mask_buffer = []
    mask_reef = []
    for point in range(len(points)):
        point_flag_buf = float("nan")
        point_flag_reef = float("nan")
        for reef_multipolygon in reefs_meters:
            distance = reef_multipolygon.distance(points_meters[point])
            distance=distance.values
            if distance == 0:
                point_flag_reef = 1
                break
            elif 1000 < distance <= 10000:#justify the range by the grid size of era5
                point_flag_buf = 1
                break
        mask_buffer.append(point_flag_buf)
        mask_reef.append(point_flag_reef)
    reef_indices = np.where(~np.isnan(mask_reef))[0]
    buf_indices = np.where(~np.isnan(mask_buffer))[0]
    return reef_indices, buf_indices, lon_flat,lat_flat, new_time, swh

def model_filter_reef_buff(reef_indices, buf_indices, lon_flat, lat_flat, ssl,ssh, swh):
    #THIS DOESNT WORK, METHOD WITH USE OF THE SAME FUNCTION AS FOR ALTIMETRY WAS BETTER(old comment?)
    reef_filtered_swh={}
    buffer_filtered_swh={}
    difference_swh = {}
    for time_model in range(len(swh)):
        swh_flat = swh[time_model].flatten()
        #replace masked array values in model
        swh_flat = ma.filled(swh_flat, float('nan'))

        lon_reef=lon_flat[reef_indices]
        lat_reef = lat_flat[reef_indices]
        swh_reef = swh_flat[reef_indices]
        lon_buf=lon_flat[buf_indices]
        lat_buf=lat_flat[buf_indices]
        swh_buf=swh_flat[buf_indices]

        swhs_reef = np.array(swh_reef)
        reef_locations = np.column_stack((lon_reef, lat_reef))
        swhs_buffer = np.array(swh_buf)
        buff_locations = np.column_stack((lon_buf, lat_buf))
        differences = np.empty(0)

        diff_all=np.empty(0)
        att_all=np.empty(0)
        diff_lon=np.empty(0)
        diff_lat=np.empty(0)


        for n in range(len(reef_locations)):
            distance = cdist([reef_locations[n]], buff_locations, metric='euclidean')
            #print(type(distance))
            min_distance_index = np.argmin(distance[0])
            closest_buffer_swh = swhs_buffer[min_distance_index]
                
            if distance[0][min_distance_index]<=1 and ssl <= closest_buffer_swh <= ssh:

                differences = closest_buffer_swh-swhs_reef[n]
                attenuation_percent=(100/closest_buffer_swh)*differences
                diff_all = np.append(diff_all, differences)
                att_all = np.append(att_all, attenuation_percent)
                diff_lon = np.append(diff_lon, lon_reef[n])
                diff_lat = np.append(diff_lat, lat_reef[n])
        
        reef_filtered_swh[time_model] = pd.DataFrame({'lon': lon_reef,
                                                    'lat': lat_reef, 
                                                    'swh': swh_reef})
                                
        buffer_filtered_swh[time_model] = pd.DataFrame({'lon': lon_buf,
                                                        'lat': lat_buf, 
                                                        'swh': swh_buf})
        difference_swh[time_model] = pd.DataFrame({'lon': diff_lon,
                                        'lat': diff_lat, 
                                        'swh': diff_all, 
                                        'attenuation_percent' : att_all})
    return reef_filtered_swh, buffer_filtered_swh, difference_swh

# %%
def modeled_track_new(model, reference_time, pass_data_dict):
    #model data
    #swh(time, latitude, longitude)
    rmse_reef=0
    corr_reef=[0, 0]
    fid=Dataset(model)
    lon=fid.variables['longitude'][:]
    lat=fid.variables['latitude'][:]
    swh=fid.variables['swh'][:]
    mwp=fid.variables['mwp'][:]
    time_model=fid.variables['time'][:]
    #print(swh)
    #print(swh[1])
    
    #sat data
    lons_sat = pass_data_dict['lon']
    lats_sat = pass_data_dict['lat']
    swhs_sat = pass_data_dict['swh']
    time_sat = pass_data_dict['time']
    #filter out all nans for correlation calculation
    '''condition = ~np.array([np.isnan(x) for x in swhs_sat])
    swhs_sat = swhs_sat[condition]
    lons_sat=lons_sat[condition]
    lats_sat=lats_sat[condition]
    time_sat=time_sat[condition]'''

    # calculate the satellite observation time in seconds
    time_sat = [(timestamp - reference_time).total_seconds() for timestamp in time_sat]
    #model time in seconds
    time_delta =[datetime.timedelta(hours=i) for i in range(len(time_model))]
    new_time_raw = [reference_time + time_delta[n] for n in  range(len(time_delta))]
    new_time=[(timestamp - reference_time).total_seconds() for timestamp in new_time_raw]
    
    sat_array = np.column_stack((time_sat, lats_sat, lons_sat))
    #reverse lat ans swh from model lists 
    reversed_lat = lat[::-1]
    
    #for n, swh_slice in enumerate(swh):
    #     swh[n]=swh_slice[::-1]
    swh= np.flip(swh, axis=1)
    mwp = np.flip(mwp, axis=1)
    #replace masked array values in model
    swh = ma.filled(swh, float('nan'))
    #for n, mwp_slice in enumerate(mwp):
    #     mwp[n]=mwp_slice[::-1]
    #replace masked array values in model
    mwp = ma.filled(mwp, float('nan'))
    interp_swh = RegularGridInterpolator((new_time, reversed_lat, lon), swh, bounds_error=False)
    swh_modeled=interp_swh(sat_array)
    interp_mwp = RegularGridInterpolator((new_time, reversed_lat, lon), mwp, bounds_error=False)
    mwp_modeled=interp_mwp(sat_array)
    #remove nan values from modelled swh to calc correlation 
    condition = ~np.array([np.isnan(x) for x in swh_modeled])
    swh_modeled_nonan = swh_modeled[condition]
    swhs_sat_nonan=swhs_sat[condition]
    if len(swhs_sat)>2:
        rmse_reef=math.sqrt(sum(np.square(swhs_sat-swh_modeled))/len(swhs_sat))
        corr_reef=pearsonr(swhs_sat, swh_modeled)
    #print(corr_reef)
    
    
    modeled_track_df=pd.DataFrame({'lon':lons_sat, 'lat':lats_sat, 'swh':swh_modeled, 'mwp':mwp_modeled, 'time':time_sat})
    
    return modeled_track_df, rmse_reef, corr_reef #, swhs_sat
#%%
def mean_of_flags_seastate_mt(lons, lats, swh, mwp, flags, sea_state_low, sea_state_high, time_sat):
    
    density = 1025  # average sea water density in kg/m³
    gravity = 9.81
    diff_all=[]
    att_all=[]
    diff_energy=[]
    att_energy=[]
    flags=np.array(flags)
    lons=np.array(lons)
    lats=np.array(lats)
    swh=np.array(swh)
    mwp=np.array(mwp)
    time_sat=np.array(time_sat)

    buffer_indx= np.where(flags == 2)[0]
    #extract buffer measurments
    lon_buf=lons[buffer_indx]
    lat_buf=lats[buffer_indx]
    swhs_buf=swh[buffer_indx]
    mwp_buf=mwp[buffer_indx]
    time_buf=time_sat[buffer_indx]

    #average the reef measurments
    mean_swhs_reef = []
    mean_lat_reef =[]
    mean_lon_reef=[]
    mean_mwp_reef=[]
    sum_swh_reef = []
    sum_lat_reef = []
    sum_lon_reef = []
    sum_mwp_reef = []
    time_sat_list=[]
    mean_time=[]
    counter_reef=0
    #calculate mean swh per reef
    for i in range(0, (len(flags))):
        if flags[i] == 1: 
            counter_reef += 1
            sum_swh_reef.append(swh[i])
            sum_lat_reef.append(lats[i])
            sum_lon_reef.append(lons[i])
            sum_mwp_reef.append(mwp[i])
            time_sat_list.append(time_sat[i])
        elif counter_reef > 0:
            mean_swhs_reef.append(np.nanmean(sum_swh_reef[1:(len(sum_swh_reef))-1]))
            mean_mwp_reef.append(np.nanmean(sum_mwp_reef[1:(len(sum_mwp_reef)-1)]))
            mean_lon_reef.append(np.nanmean(sum_lon_reef))
            mean_lat_reef.append(np.nanmean(sum_lat_reef))
            mean_time.append(time_sat_list[0])
            counter_reef = 0
            sum_swh_reef = []
            sum_lat_reef = []
            sum_lon_reef = []
            sum_mwp_reef = []
    if not len(mean_swhs_reef)==0 and not len(swhs_buf)==0:
        reef_locations=np.column_stack((mean_lon_reef, mean_lat_reef))
        buff_locations = np.column_stack((lon_buf, lat_buf))
        for n in range(len(reef_locations)):
            distance = cdist([reef_locations[n]], buff_locations, metric='euclidean')
            min_distance_index = np.argmin(distance[0])
            closest_buffer_swh = swhs_buf[min_distance_index]
            closest_buffer_mwp = mwp_buf[min_distance_index]
            distance[0][min_distance_index]=1000
            seconf_min_indx=np.argmin(distance[0])
            if swhs_buf[seconf_min_indx]>closest_buffer_swh:
                 closest_buffer_swh=swhs_buf[seconf_min_indx]
                 closest_buffer_mwp = mwp_buf[seconf_min_indx]
            #calculate energy for chose buff and reef location
            energy_reef = density * (gravity**2) * (mean_swhs_reef[n]**2) * mean_mwp_reef[n]/(64*np.pi*1000)
            energy_buf = density * (gravity**2) * (closest_buffer_swh**2) * closest_buffer_mwp/(64*np.pi*1000)
            if sea_state_low <= closest_buffer_swh <=sea_state_high:
                differences = closest_buffer_swh-mean_swhs_reef[n]
                differences_energy = energy_buf-energy_reef
                attenuation_percent=(100/closest_buffer_swh)*differences
                att_percent_energy=(100/energy_buf)*differences_energy
                diff_all.append(differences)
                att_all.append(attenuation_percent)
                diff_energy.append(differences_energy)
                att_energy.append(att_percent_energy)

    attenuation_stats= pd.DataFrame({'diff_swh':diff_all,
                                     'att_swh':att_all,
                                     'diff_energy':diff_energy,
                                     'att_energy': att_energy})

    reef_filtered_swh = pd.DataFrame({'lon': mean_lon_reef,
                                      'lat': mean_lat_reef, 
                                      'swh':mean_swhs_reef,
                                      'time':mean_time})
    buffer_filtered_swh = pd.DataFrame({'lon': lon_buf,
                                        'lat': lat_buf, 
                                        'swh': swhs_buf,
                                        'time':time_buf})
    return buffer_filtered_swh, reef_filtered_swh, attenuation_stats

def filter_by_flags_seastate_mt(pass_data_dict, pass_flags_sats, pass_nums, satellite_number, ssl, ssh):
    ''' This function is used to sepearete observations on reef and observations before reef utilizang the flags calculated beforehand
    Input: pass_data_dict - dictionary with all observation for one satellite for one month sorted by passes
    Output:
    '''
    reef_filtered = {}
    buffer_filtered = {}
    difference_swh = {}
    #attenuation_stats = {}
    pass_flag=pass_flags_sats[satellite_number]  #choose flags for specified satellite
    for num in pass_nums:
        if (num in pass_flag) and (num in pass_data_dict) and not pass_data_dict[num].empty and not pass_flag[num].empty:
            #read observation data from sorted by pass number dictionary
            lons = pass_data_dict[num]['lon'].tolist() 
            lats = pass_data_dict[num]['lat'].tolist() 
            swhs = pass_data_dict[num]['swh'].tolist()
            mwp = pass_data_dict[num]['mwp'].tolist()
            time_sat=pass_data_dict[num]['time'].tolist()
            #read flag data from calculated dictionary
            flags= pass_flag[num]['flag']
            lats_flag=pass_flag[num]['lat']
            
            flag_new = griddata(lats_flag, flags, lats, method='linear', fill_value=0)#works perfect
            flag_new=flag_new.round()
            calculations=False
            for i in range(5, len(flag_new)):
                if flag_new[i] == 2 and  (not any(value == 1 for value in flag_new[i-5:i]) or not any(value == 1 for value in flag_new[i+1:i+6])):
                    check=1
                elif flag_new[i]==1:
                    check=2
                    calculations=True
                else:
                    flag_new[i]=0

            

            #if calculations:
            #calculate mean point of each reef and each buffer included in one pass        
            buffer_filtered_num, reef_filtered_num, difference_num = mean_of_flags_seastate_mt(lons, lats, swhs, mwp, flag_new, ssl, ssh, time_sat)
            reef_filtered [num] = reef_filtered_num
            buffer_filtered[num]=buffer_filtered_num
            #attenuation_stats[num] = attenuation_stats_num
            difference_swh[num] = difference_num
            
        
    return buffer_filtered, reef_filtered, difference_swh
#%%
def modeled_track(model, reference_time, target_time, lons_sat, lats_sat):
     
    fid=Dataset(model)
    lon=fid.variables['longitude'][:]
    lat=fid.variables['latitude'][:]
    swh=fid.variables['swh'][:]
    mwp=fid.variables['mwp'][:]
    time_model=fid.variables['time'][:]

    time_delta =[datetime.timedelta(hours=i) for i in range(len(time_model))]

    new_time = [reference_time + time_delta[n] for n in  range(len(time_delta))]
    
    for n in range(len(new_time)):
         if new_time[n].day == target_time.day and new_time[n].hour == target_time.hour:
            start_index = n
            end_index = start_index + 1

    # Extract the relevant portion of swh based on the desired indices
    swh_slice = swh[start_index]
    mwp_slice = mwp[start_index]
    print(len(swh_slice))
    print(len(swh_slice[0]))
    

    sat_array = np.column_stack((lats_sat, lons_sat))
    #reverse lat ans swh from model lists 
    reversed_lat = lat[::-1]
    
    swh_slice= np.flip(swh_slice, axis=0)
    mwp_slice = np.flip(mwp_slice, axis=0)
    #replace masked array values in model
    swh_slice = ma.filled(swh_slice, np.float('nan'))
    #replace masked array values in model
    mwp_slice = ma.filled(mwp_slice, np.float('nan'))
    interp_swh = RegularGridInterpolator((reversed_lat, lon), swh_slice, bounds_error=False, fill_value=None)
    swh_modeled=interp_swh(sat_array)
    interp_mwp = RegularGridInterpolator((reversed_lat, lon), mwp_slice, bounds_error=False, fill_value=None)
    mwp_modeled=interp_mwp(sat_array)

    '''lon_mesh, lat_mesh = np.meshgrid(lon, lat)
    # Flatten the arrays 
    lon_flat = lon_mesh.flatten()
    lat_flat = lat_mesh.flatten()
    swh_flat = swh_slice.flatten()
    mwp_flat = mwp_slice.flatten()
    #filter out masked points
    condition = ~np.array([isinstance(x, np.ma.core.MaskedConstant) for x in swh_flat])
    # Apply the condition to remove masked values
    swh_f = swh_flat[condition]
    mwp_f=mwp_flat[condition]
    lon_f=lon_flat[condition]
    lat_f=lat_flat[condition]
    #model track
    #track_points= np.column_stack((lons_sat, lats_sat))
    #model_points= np.column_stack((lon_flat,lat_flat))
    track_points= list(zip(lons_sat, lats_sat))
    model_points= list(zip(lon_f,lat_f))
    swh_modeled = griddata(model_points, swh_f, track_points, method='linear', fill_value=0)#works perfect
    mwp_modeled = griddata(model_points, mwp_f, track_points, method='linear', fill_value=0)#works perfect
    #interpolation'''
    '''distances = cdist(model_points, track_points, metric='euclidean')
    # Find the index of the closest point for each satellite track point
    closest_indices = np.argmin(distances, axis=0)
    lon_modeled = lon_flat[closest_indices]
    lat_modeled = lat_flat[closest_indices]
    swh_modeled = swh_flat[closest_indices]
    mwp_modeled = mwp_flat[closest_indices]'''
    
    modeled_track_df=pd.DataFrame({'lon':lons_sat, 'lat':lats_sat, 'swh':swh_modeled, 'mwp':mwp_modeled})
    
    return modeled_track_df#,lon_flat, lat_flat, swh_flat #lons_sat, lats_sat, swh_modeled, mwp_modeled
