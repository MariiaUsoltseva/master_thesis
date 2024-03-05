import xarray as xr
import pandas as pd
import os
import geopandas as gpd
from shapely.geometry import Point, Polygon
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import path
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata
from scipy.spatial.distance import cdist
import datetime
from netCDF4 import Dataset
import math
import warnings
#%% list of all pass numbers
#folder for current analisys

def netcdf_to_dataframe(netcdf_file):
    '''Function to convert netCDF file from altimetry dataset to dataframe without losing the acquisition time'''
    
    netds=Dataset(netcdf_file)
    #convert seconds from reference epoch to  data time of data acquisition
    timevar=netds.variables['time'][:]
    seconds = timevar[0]
    reference_epoch = datetime.datetime(1981, 1, 1, 0, 0, 0)
    time_obs = reference_epoch + datetime.timedelta(seconds=seconds)

    lon=netds.variables['lon'][:]
    lat=netds.variables['lat'][:]
    swh_denoised=netds.variables['swh_denoised'][:]
    satellite=netds.variables['satellite'][:]
    pass_num=netds.variables['relative_pass_number'][:]
    dist_2_coast=netds.variables['distance_to_coast']

    df=pd.DataFrame({'lon': lon, 'lat':lat, 'swh_denoised': swh_denoised,
                    "satellite": satellite, "relative_pass_number": pass_num, 
                    'time': time_obs, 'dist_2_coast': dist_2_coast})
    return(df)

def sort_by_pass(month_folder,satellite_number, lonmin, lonmax, latmin, latmax):
    '''the function sorts all data from netcdf files stored in month_folder with days subfolders 
    by number of satellite pass (and number of satellite?) and stores the lon.lat and swh in the dictionary 
    per pass number. Only data between 35 and - 35 degrees latitude and closere than 500 km to cast is saved.
    Input: path to the folder 
            satellite_number - 0 - CryoSat2, 1 - Jason-1, 2 - Jason-2, 3 - Jason-3, 4 - Saral, 5 - Sentinel-3A, 6 - Envisat
            corner coordinates of the region of interest
    Outut: list of relative pass numbers 
            dicionary sorted by pass number '''
    all_files=[]
    pass_set = set()
    days=os.listdir(month_folder)
    month=month_folder[25:27]
    year=month_folder[20:24]
    for day in days:
        file='D:/thesis/ceda_data/{}/{}/{}/ESACCI-SEASTATE-L3-SWH-MULTI_1D-{}{}{}-fv01.nc'.format(year,month,day,year, month, day)
        all_files.append(file)

    for file in all_files:
        df=netcdf_to_dataframe(file)
        pass_data=df[df['satellite']== satellite_number]
        #create a list of pass numbers
        pass_set.update(pass_data['relative_pass_number'].unique())
    pass_nums = sorted(list(pass_set))

    pass_data_dict = {}
    for file in all_files:
        df=netcdf_to_dataframe(file)
        sat_data=df[df['satellite'] == satellite_number]
        for num in pass_nums:
            pass_data = sat_data[sat_data['relative_pass_number'] == num]
            lon_s =np.array(pass_data['lon'])
            lat_s= np.array(pass_data['lat'])
            swh_s=np.array(pass_data['swh_denoised'])
            time_sat_s=np.array(pass_data['time'])
            dist_2_coast_s=np.array(pass_data['dist_2_coast'])
            #condition to exclude all data too far from the coast
            # and not included in the region of interest
            close_indices = (lat_s < latmax) & (lat_s > latmin) & (lon_s < lonmax) & (lon_s > lonmin)
            lon_close = lon_s[close_indices]
            lat_close = lat_s[close_indices]
            swh_close = swh_s[close_indices]
            dist_close = dist_2_coast_s[close_indices]
            time_close = time_sat_s[close_indices]

            if not len(lon_close) ==0:
                if num not in pass_data_dict:
                    pass_data_dict[num] = pd.DataFrame({'lon': lon_close,
                                                        'lat': lat_close, 'swh': swh_close,
                                                          'dist_coast': dist_close,
                                                          'time': time_close})
                else:
                    pass_data_dict[num] = pd.concat([pass_data_dict[num], pd.DataFrame({'lon': lon_close,
                                                                                         'lat': lat_close,
                                                                                         'swh': swh_close,
                                                                                         'dist_coast': dist_close,
                                                                                         'time': time_close})])
    return pass_nums, pass_data_dict

#%%mean swh per coordinate
def mean_swh(pass_reef_dict):
    '''The function calculates mean swh per unique location(lon,lat) for each nominal track
    input: pass_reef_dict - sorted by pass number dictionary, generated in sort_by_pass function
    output: mean_swh_dict - dictionary with mean swh value per location, sorted by relative pass number'''
    mean_swh_dict={}
    for num in pass_reef_dict:
        lons=pass_reef_dict[num]['lon'].tolist()
        lats=pass_reef_dict[num]['lat'].tolist()
        swhs=pass_reef_dict[num]['swh'].tolist()
        #dist=pass_reef_dict[num]['dist_coast'].tolist()
        location_wave_heights={}
        for i in range(len(lons)):
            #if location is already in the dictionary add location to a running sum
            if (lons[i],lats[i]) in location_wave_heights:
                location_wave_heights[(lons[i], lats[i])]['sum_swh'] += swhs[i]
                location_wave_heights[(lons[i], lats[i])]['count'] += 1
            else:
            # If the location isn't in the dictionary, add it with corresponding measurments
                location_wave_heights[(lons[i], lats[i])] = {'sum_swh': swhs[i],
                                                            'count': 1} #'dist_coast': dist[i]}
                    
    mean_swh=[]
    mean_lon=[]
    mean_lat=[]
    # Iterate through each unique location in the pass
    for unique_location, location_data in location_wave_heights.items():
        # Calculate the mean wave height for the location
        mean_wave_height = location_data['sum_swh'] / location_data['count']
        #distance = location_data['dist_coast']
        mean_lon.append(unique_location[0])
        mean_lat.append(unique_location[1])
        mean_swh.append(mean_wave_height)
        
    mean_swh_dict = pd.DataFrame({'lon': mean_lon,
                                            'lat': mean_lat, 
                                            'swh': mean_swh}) 
                                            #'dist_coast': distance})
    return(mean_swh_dict)


#%%plots
def geo_plots(reef_boundary,lonmin, lonmax, latmin, latmax):
    '''Plot map of region of interest specified reef boundaries
    input: reef_boundary - list of reef boundaries in geopandas format
            corner coordinates of region of interest'''
    fig = plt.figure(figsize=[15,15])
    ax = fig.add_subplot(1, 1, 1) 
    map_base = Basemap(projection='cyl',llcrnrlat= latmin,urcrnrlat=latmax,
                llcrnrlon=lonmin, urcrnrlon=lonmax,resolution='c', ax=ax, suppress_ticks=False)
    map_base.drawcoastlines()
    plt.xlabel("Longitude", fontsize=14)
    plt.ylabel("Latitude", fontsize=14)
#   print the reef extends
    for reef in reef_boundary:
        reef.plot(ax=ax, color='dimgray')

    return fig, ax, map_base
#%%
def mean_of_flags_seastate_new(lons, lats, swh, flags, sea_state_low, sea_state_high, time_sat):
    #diff_all=[0,0]
    #att_all=[0,0] does bad thing with array lenghts
    
    flags=np.array(flags)
    lons=np.array(lons)
    lats=np.array(lats)
    swh=np.array(swh)
    time_sat=np.array(time_sat)

    reef_indx=np.where(flags == 1)[0]
    buffer_indx= np.where(flags == 2)[0]

    lon_reef=lons[reef_indx]
    lat_reef = lats[reef_indx]
    swhs_reef = swh[reef_indx]
    time_reef=time_sat[reef_indx]
    lon_buf=lons[buffer_indx]
    lat_buf=lats[buffer_indx]
    swhs_buf=swh[buffer_indx]
    time_buf=time_sat[buffer_indx]
    
    #average the reef measurments
    mean_swhs_reef = []
    mean_lat_reef =[]
    mean_lon_reef=[]
    sum_swh_reef = []
    sum_lat_reef = []
    sum_lon_reef = []
    time_sat_list=[]
    mean_time=[]
    counter_reef=0
    #calculate mean swh per reef
    #print('flags type', type(flags))
    for i in range(0, (len(flags))):
        if flags[i] == 1: 
            counter_reef += 1
            sum_swh_reef.append(swh[i])
            sum_lat_reef.append(lats[i])
            sum_lon_reef.append(lons[i])
            time_sat_list.append(time_sat[i])
        elif counter_reef > 0:
            #if np.any(~np.isnan(np.array(sum_swh_reef))): #condition to get read of runtime warrning (mean of empty slice)
            #suppres empty slice warning
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
            #    foo = np.nanmean(sum_swh_reef)

                mean_swhs_reef.append(np.nanmean(sum_swh_reef)) #[1:(len(sum_swh_reef)-1)]))
                mean_lon_reef.append(np.nanmean(sum_lon_reef))
                mean_lat_reef.append(np.nanmean(sum_lat_reef))
                mean_time.append(time_sat_list[0])

            counter_reef = 0
            sum_swh_reef = []
            sum_lat_reef = []
            sum_lon_reef = []
    #diff_all=np.zeros(len(mean_swhs_reef))
    #att_all=np.zeros(len(mean_swhs_reef))
    diff_all=np.empty(0)
    att_all=np.empty(0)
    diff_lon=np.empty(0)
    diff_lat=np.empty(0)
    if not len(mean_swhs_reef)==0 and not len(swhs_buf)==0:
        reef_locations=np.column_stack((mean_lon_reef, mean_lat_reef))
        buff_locations = np.column_stack((lon_buf, lat_buf))
        for n in range(len(reef_locations)):
            distance = cdist([reef_locations[n]], buff_locations, metric='euclidean')
            #print(type(distance))
            min_distance_index = np.argmin(distance[0])
            closest_buffer_swh = swhs_buf[min_distance_index]
            distance[0][min_distance_index]=1000
            seconf_min_indx=np.argmin(distance[0])
            if swhs_buf[seconf_min_indx]>closest_buffer_swh:
                closest_buffer_swh = swhs_buf[seconf_min_indx]
        
            if sea_state_low <= closest_buffer_swh <=sea_state_high:
                differences = closest_buffer_swh-mean_swhs_reef[n]
                attenuation_percent=(100/closest_buffer_swh)*differences
                diff_all = np.append(diff_all, differences)
                att_all = np.append(att_all, attenuation_percent)
                diff_lon = np.append(diff_lon, mean_lon_reef[n])
                diff_lat = np.append(diff_lat, mean_lat_reef[n])
                #diff_all[n]=differences
                #att_all[n]=attenuation_percent
            #else:
                #diff_all[n]=float('nan')
                #att_all[n]=float('nan')
                
    reef_filtered_swh = pd.DataFrame({'lon': mean_lon_reef,
                                      'lat': mean_lat_reef, 
                                      'swh':mean_swhs_reef,
                                      'time':mean_time})
    buffer_filtered_swh = pd.DataFrame({'lon': lon_buf,
                                        'lat': lat_buf, 
                                        'swh': swhs_buf,#})
                                        'time':time_buf})
    difference_swh = pd.DataFrame({'lon': diff_lon,
                                        'lat': diff_lat, 
                                        'swh': diff_all, 
                                        'attenuation_percent' : att_all})
    return buffer_filtered_swh, reef_filtered_swh, difference_swh

def filter_by_flags_seastate(pass_data_dict, pass_flags_sats, pass_nums, satellite_number, ssl, ssh):
    ''' This function is used to sepearete observations on reef and observations before reef utilizang the flags calculated beforehand
    Input: pass_data_dict - dictionary with all observation for one satellite for one month sorted by passes
    Output:
    '''
    diff_all=[0,0]
    att_all=[0,0]
    reef_filtered = {}
    buffer_filtered = {}
    difference_swh = {}
    pass_flag=pass_flags_sats[satellite_number]  #choose flags for specified satellite
    for num in pass_nums:
        if (num in pass_flag) and (num in pass_data_dict) and not pass_data_dict[num].empty and not pass_flag[num].empty:
            #read observation data from sorted by pass number dictionary
            lons = pass_data_dict[num]['lon'].tolist() 
            lats = pass_data_dict[num]['lat'].tolist() 
            swhs = pass_data_dict[num]['swh'].tolist()
            time_sat=pass_data_dict[num]['time'].tolist()
            #read flag data from calculated dictionary
            flags= pass_flag[num]['flag']
            lats_flag=pass_flag[num]['lat']
            
            flag_new = griddata(lats_flag, flags, lats, method='linear', fill_value=0)
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
            buffer_filtered_num, reef_filtered_num, difference_num = mean_of_flags_seastate_new(lons, lats, swhs, flag_new, ssl, ssh, time_sat)
            reef_filtered[num] = reef_filtered_num
            buffer_filtered[num] = buffer_filtered_num
            difference_swh[num] = difference_num
            
    return buffer_filtered, reef_filtered, difference_swh

def colorbars(fig_r, fig_b, fig_d, ax_r, ax_b, ax_d, hmax, hmin):
    norm = mpl.colors.Normalize(hmin, hmax) #vmax=max(swh))# normalizes all values existing in swh to fit the colorbar
    col_set=cm.ScalarMappable(norm=norm, cmap=plt.cm.jet) #gives colors and values to colorbar
    clb = fig_r.colorbar(col_set, ax=ax_r, orientation='horizontal', pad=0.05, aspect=35)
    clb.set_label('SWH [m]', fontsize=16)
    clb.ax.get_yaxis().label.set_rotation(-90)

    norm = mpl.colors.Normalize(hmin, hmax) #vmax=max(swh))# normalizes all values existing in swh to fit the colorbar
    col_set=cm.ScalarMappable(norm=norm, cmap=plt.cm.jet) #gives colors and values to colorbar
    clb_buf=fig_b.colorbar(col_set, ax=ax_b, orientation='horizontal', pad=0.05, aspect=35)
    clb_buf.set_label('SWH [m]', fontsize=16)
    clb_buf.ax.get_yaxis().label.set_rotation(-90)

    norm = mpl.colors.Normalize(vmin=-1, vmax=1) #vmax=max(swh))# normalizes all values existing in swh to fit the colorbar
    col_set=cm.ScalarMappable(norm=norm, cmap=plt.cm.bwr) #gives colors and values to colorbar
    clb_diff = fig_d.colorbar(col_set, ax=ax_d, orientation='horizontal', pad=0.05, aspect=35)
    clb_diff.set_label('difference [m]', fontsize=16)
    clb_diff.ax.get_yaxis().label.set_rotation(-90)

    return clb, clb_buf, clb_diff




def filter_by_flags_model(mean_swh_dict, pass_flags_sats, pass_num, satellite_number, ssl, ssh):
    reef_filtered_swh = {}
    buffer_filtered_swh = {}
    pass_flag=pass_flags_sats[satellite_number]
    #print(pass_num)
    if (pass_num in pass_flag) and (pass_num in mean_swh_dict):
        #print('check')
        lons = mean_swh_dict['lon'].tolist() 
        lats = mean_swh_dict['lat'].tolist() 
        swhs = mean_swh_dict['swh'].tolist()
        #print(pass_flags_sats)
        
        pass_flags_sats=pass_flags_sats[satellite_number]
        flags= pass_flags_sats[pass_num]['flag'].tolist()
        lats_flag=pass_flags_sats[pass_num]['lat'].tolist() 
        lons_flag=pass_flags_sats[pass_num]['lon'].tolist() 
        #define interpolation func
        #inter_func=interp2d(lons_flag, lats_flag, flags, bounds_error=False, fill_value=0)
        flag_new= griddata((lons_flag, lats_flag), flags, (lons, lats),method='linear',fill_value=0)                    
        flag_new=flag_new.round()
        if flag_new.any() == 1:            
            #filter out bad locations on buffer
            filtered_indexes_buf = []
            for i in range(2, len(flag_new)):#this works
                #if flag_new[i] == 2 and (all(value == 0 for value in flag_new[i-30:i]) or all(value == 0 for value in flag_new[i+1:i+30])): 
                if flag_new[i] == 2 and  (not any(value == 1 for value in flag_new[i-2:i]) or not any(value == 1 for value in flag_new[i+1:i+3])):              
                    filtered_indexes_buf.append(i)
                elif flag_new[i]==1:
                    check=2        
                else:
                    flag_new[i]=0
                        
                            
            lon_reef, lat_reef, swhs_reef, lon_buf, lat_buf, swhs_buf = mean_of_flags_seastate(lons, lats, swhs, flag_new,ssl, ssh)
            reef_filtered_swh = pd.DataFrame({'lon': lon_reef,
                                                'lat': lat_reef, 
                                                'swh': swhs_reef})
                            
            buffer_filtered_swh = pd.DataFrame({'lon': lon_buf,
                                                    'lat': lat_buf, 
                                                    'swh': swhs_buf})
    return(reef_filtered_swh, buffer_filtered_swh)

def haversine(lat1, lon1, lat2, lon2):
    # Radius of the Earth in kilometers
    earth_radius = 6371  # You can also use 3958.8 miles for distance in miles
    
    # Convert latitude and longitude from degrees to radians
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    
    # Haversine formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.asin(math.sqrt(a))
    
    # Calculate the distance
    distance = earth_radius * c
    
    return distance 
