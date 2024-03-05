#%%
import pandas as pd
import os
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import pickle
import time
from netCDF4 import Dataset
from common_settings import lonmin, lonmax, latmin, latmax
from common_settings import satellites, sat_flags
from common_settings import months, years, information
from model_functions import model_filter_flags, model_filter_reef_buff
#%%
#upload ERA5 data
model='D:/thesis/model_data/2018model3hours.nc'

#upload reef polygons for plots
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

differences_all=[]
attenuation_perecent_all=[]
model_ss={}

#sea states 3-7 according to WMO sea state code
ssl=[0.5] #, 1.25, 2.5, 4, 6] #lower borders
ssh=[1.25]#, 2.5, 4, 6, 9] #upper borders 
fid = Dataset(model)
swh=fid.variables['swh'][:]
reef_indices, buf_indices, lon_flat,lat_flat, new_time, swh = model_filter_flags(model,reefs_meters)
for ss in range(len(ssl)):
    print(ss)
    differences_all=[]
    attenuation_perecent_all=[]
    reef_filtered_swh, buffer_filtered_swh, difference_swh = model_filter_reef_buff(reef_indices, buf_indices, lon_flat, lat_flat, ssl,ssh, swh)
    #for i in range(len(differences_all)):
    #    if differences_all[i]>100:
    #        differences_all[i]=float('nan')
    #        attenuation_perecent_all[i]=float('nan')

    for time_model in reef_filtered_swh:
        diff_model = difference_swh[time_model]['swh']
        att_model = difference_swh[time_model]['attenuation_percent']
        differences_all.extend(diff_model)
        attenuation_perecent_all.extend(att_model)
        
    model_ss[ss]= pd.DataFrame({'differences': differences_all,'attenuation': attenuation_perecent_all})


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
        
    diff_model = model_ss[sea_state]['differences'].tolist() #sat
    att_model = model_ss[sea_state]['attenuation'].tolist()
    if len(diff_model)==0:
        diff_model=[0]

    ax[a][b].hist(diff_model, bins=bins_n)
    ax[a][b].title.set_text(states[sea_state])

    ax[a][b].grid(True)    
    b=b+1
    count=0
    for i in diff_model:
        if i>=0:
            count+=1
    pos_diff.append((100/len(diff_model))*count)
    
    sum_per=[]
    count=0
    for n in range(len(att_model)):
        if not (np.isnan(att_model[n])) and att_model[n] > -100:
            sum_per.append(att_model[n])
    mean_att.append(np.mean(sum_per))
    mean_std.append(np.nanstd(diff_model))
    mean_att_meters.append(np.nanmean(diff_model))

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
#plt.title("Effect of reefs on different sea states")
plt.show()
#%%
#save data
with open('model_ss_final.pickle', 'wb') as handle:
    pickle.dump(model_ss, handle, protocol=pickle.HIGHEST_PROTOCOL)
#%%