import pandas as pd
import os
import geopandas as gpd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pickle
import time
from common_settings import lonmin, lonmax, latmin, latmax
from common_settings import satellites, sat_flags
from common_settings import months, years, information
from poly_funcs import sort_by_pass, filter_by_flags_seastate
from netCDF4 import Dataset
#%%
#variables
#calculate for south asia in year 2003(1,6) and 2015(0,2,4)
#australia 2008(1,2,6) and 2018(2,3,4,5)
#working file for sat state statistics derived from satellite data for 1 year (2017)
#month_folder='D:/thesis/ceda_data/2017/08'
#chek new differences closest point function
#East asia 2009: 40% 2017: 33.88%  
#moorea (should be hightly energetic) in the beginning of altimetry era and now
#%%
print(information)
#satellite_numbers=[3,5]
#months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
reef_folder="D:/thesis/reefs"  #folder contains reef polygons
reef_file=os.listdir(reef_folder)
reefs=[]
for file in reef_file:
    full_name=os.path.join(reef_folder, file)
    reefs.append(gpd.read_file(full_name))
polygons = reefs

sat_ss={}
std_ss=[]
mean_ss=[]
sat_ss_09 ={}
#%%
start_time=time.time()
#sea states 3-7 according to WMO sea state code
ssl=[0.5, 1.25, 2.5, 4, 6] #lower boundaries
ssh=[1.25, 2.5, 4, 6, 9] #upper boundaries

for ss in range(len(ssl)):
    differences_all=[]
    attenuation_perecent_all=[]
    for year in years:
        for month in months:
            print('Processing month = ', month)
            month_folder='D:/thesis/ceda_data/{}/{}'.format(year, month)
            for i in range(len(satellites)):
                pass_nums, pass_data_dict = sort_by_pass(month_folder, satellites[i], lonmin, lonmax, latmin, latmax)
                buffer_filtered, reef_filtered, difference_swh = filter_by_flags_seastate(pass_data_dict, sat_flags, pass_nums, i, ssl[ss],ssh[ss])        
                for num in difference_swh:
                    swhs_diff = difference_swh[num]['swh'].tolist()
                    att_percent = difference_swh[num]['attenuation_percent'].tolist()
                    differences_all.extend(swhs_diff)
                    attenuation_perecent_all.extend(att_percent)
     

    sat_ss[ss]= pd.DataFrame({'differences': differences_all,'attenuation': attenuation_perecent_all})

end_time=time.time()
print("Execution time = ", (end_time-start_time)/60 , 'minutes')
#%%
states=["Sea state 3 (0.5-1.25 m)", "Sea state 4 (1.25-2.5 m)", "Sea state 5 (2.5-4 m)", "Sea state 6 (4-6 m)", "Sea state 7 (6-9 m)"]
bins_nums=[100, 100, 50, 20, 10]
fig, ax= plt.subplots(2, 3, figsize=(15, 11))
a=0
b=0
pos_diff=[]
mean_att=[]
mean_att_meters=[]
mean_std=[]
for sea_state in range(len(ssl)):
    if sea_state ==3:
        a=1
        b=0
    bins_n=bins_nums[sea_state]
        
    diff_sat = sat_ss[sea_state]['differences'].tolist() #sat
    att_sat = sat_ss[sea_state]['attenuation'].tolist()
    if len(diff_sat)==0:
        diff_sat=[0]

    ax[a][b].hist(diff_sat, bins=bins_n)
    ax[a][b].title.set_text(states[sea_state])

    ax[a][b].grid(True)    
    b=b+1
    count=0
    for i in diff_sat:
        if i>=0:
            count+=1
    pos_diff.append((100/len(diff_sat))*count)
    
    sum_per=[]
    count=0
    for n in range(len(att_sat)):
        if not (np.isnan(att_sat[n])) and att_sat[n] > -100:
            sum_per.append(att_sat[n])
    mean_att.append(np.mean(sum_per))
    mean_std.append(np.nanstd(diff_sat))
    mean_att_meters.append(np.nanmean(diff_sat))

ax[1][2].axis('off')

#plots
sea_states=['3', '4', '5', '6', '7']
n=len(sea_states)
r=np.arange(n)
fig = plt.figure(figsize = (10, 5))
#pos_diff_alt[0]=50.012343

plt.bar(r -0.1, pos_diff, 0.2, alpha=0.8, facecolor='#773EA9')

plt.bar(r + 0.1, mean_att, 0.2, alpha=0.8, facecolor='#F9D405')

plt.xticks(r,['0.5-1.25 m', '1.25-2.5 m', '2.5-4 m', '4-6 m', '6-9 m'])
plt.xlabel("Sea states")
plt.ylabel("[%]")
plt.legend(['positive differences','mean attenuation'], loc='upper left')
plt.title("Effect of reefs on different sea states")
plt.show()

#%% save data to analyse later with model stats
with open('stats_alt.pickle', 'wb') as handle:
    pickle.dump(sat_ss , handle, protocol=pickle.HIGHEST_PROTOCOL)

#%% 