import pandas as pd
import os
import geopandas as gpd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pickle
import time
from scipy.spatial.distance import cdist
from scipy.interpolate import griddata
from common_settings import lonmin, lonmax, latmin, latmax
from common_settings import satellites, sat_flags
from common_settings import months, years, information
import datetime
from poly_funcs import sort_by_pass, geo_plots, filter_by_flags_seastate, filter_by_flags_model
from model_functions import modeled_track_new, filter_by_flags_seastate_mt
#%%
#upload ERA5 data
model='D:/thesis/model_data/model_2018_03.nc'
#month_folder='D:/thesis/ceda_data/2017/08'

#upload reef polygons for plots
reef_folder="D:/thesis/reefs" 
reef_file=os.listdir(reef_folder)
reefs=[]
for file in reef_file:
    full_name=os.path.join(reef_folder, file)
    reefs.append(gpd.read_file(full_name))
polygons = reefs
ssl=[0.5, 1.25, 2.5, 4, 6]
ssh=[1.25, 2.5, 4, 6, 9]
density = 1025  # Water density in kg/m³
gravity = 9.81

#%%

#east asian seas : 90.7, 141.5, -10,  26.8
#pacific: 150, 180, -30,  30 and -180, -120, -30,  30
#australia: 108, 158, -30,  -10
#carribean: -89.2,-55.65, 9.96, 29.26
#west indian ocean : 34.9, 68.4, -25.8, 1.9
#red sea 30.1, 55.9, 0.7, 30
#south asia:66, 87.2, -7, 24.6
#et pacific: -119.9, -75.7, -2, 25
start_time=time.time()
energy_modeled_tr={}
swh_modeled_tr={}
corr=[]
ssl=[0.5, 1.25, 2.5, 4, 6]
ssh=[1.25, 2.5, 4, 6, 9]
month_folder='D:/thesis/ceda_data/2018/{}'.format(month)
for ss in range(len(ssl)):
    print("Processing sea state ", ss)
    differences_swh=[]
    attenuation_perecent_swh=[]
    differences_energy=[]
    attenuation_perecent_energy=[]
    corr_ss=[]
for month in months:
        print('Processing month ', month)
        month_folder='D:/thesis/ceda_data/2018/{}'.format(month)
        model='D:/thesis/model_data/model_2018_{}.nc'.format(month)
        reference_time = datetime.datetime(2018, 5, 1, 0, 0, 0)
        
        mt_dict={}
        for i in range(len(satellites)):      
                print(i)       
                pass_nums, pass_data_dict = sort_by_pass(month_folder, satellites[i], -180, 180, -30, 30) 
                for num in pass_data_dict:
                    modeled_track_df, rmse, corr = modeled_track_new(model, reference_time, pass_data_dict[num])  
                    mt_dict[num]=modeled_track_df
                    corr_ss.append(corr[0])
                buffer_filtered_mt, reef_filtered_mt, difference_swh_mt = filter_by_flags_seastate_mt(mt_dict, sat_flags, pass_nums, i, ssl=1, ssh=100)
               
                for num in difference_swh_mt:
                        swhs_diff = difference_swh_mt[num]['diff_swh'].tolist()
                        att_percent = difference_swh_mt[num]['att_swh'].tolist()
                        diff_energy =  difference_swh_mt[num]['diff_energy'].tolist()
                        att_energy =   difference_swh_mt[num]['att_energy'].tolist()
                        differences_swh.extend(swhs_diff)
                        attenuation_perecent_swh.extend(att_percent)
                        differences_energy.extend(diff_energy)
                        attenuation_perecent_energy.extend(att_energy)
                
        energy_modeled_tr[ss] = pd.DataFrame({'differences': differences_energy,
                                            'attenuation': attenuation_perecent_energy})     
        swh_modeled_tr[ss] = pd.DataFrame({'differences': differences_swh,
                                            'attenuation': attenuation_perecent_swh})
        #corr.append(corr_ss)

end_time=time.time()
print("Execution time = ", (end_time-start_time)/60 , 'minutes')
 #%%caclulate stats
#%%
states=["Sea state 3 (0.5-1.25 m)", "Sea state 4 (1.25-2.5 m)", "Sea state 5 (2.5-4 m)", "Sea state 6 (4-6 m)", "Sea state 7 (6-9 m)"]
bins_nums=[100, 100, 50, 20, 10]
fig, ax= plt.subplots(2, 3, figsize=(15, 11))
fig_energy, ax_energy= plt.subplots(2, 3, figsize=(15, 11))

a=0
b=0
pos_diff=[]
mean_att=[]
mean_att_meters=[]
mean_std=[]

pos_diff_energy=[]
mean_att_energy=[]
mean_att_meters_energy=[]
mean_std_energy=[]
for sea_state in range(len(ssl)):
    if sea_state ==3:
        a=1
        b=0
    bins_n=bins_nums[sea_state]
        
    diff_swh_mt = swh_modeled_tr[sea_state]['differences'].tolist() #sat
    att_swh_mt = swh_modeled_tr[sea_state]['attenuation'].tolist()

    diff_energy_mt = energy_modeled_tr[sea_state]['differences'].tolist()
    att_energy_mt = energy_modeled_tr[sea_state]['attenuation'].tolist()

    
    ax[a][b].hist(diff_swh_mt, bins=bins_n)
    ax[a][b].title.set_text(states[sea_state])
    ax[a][b].grid(True)    

    ax_energy[a][b].hist(diff_energy_mt, bins=bins_n)
    ax_energy[a][b].title.set_text(states[sea_state])
    ax_energy[a][b].grid(True)  

    b=b+1

    count=0
    for i in diff_swh_mt:
        if i>=0:
            count+=1
    pos_diff.append((100/len(diff_swh_mt))*count)
    
    sum_per=[]
    count=0
    for n in range(len(att_swh_mt)):
        if not (np.isnan(att_swh_mt[n])) and att_swh_mt[n] > -100:
            sum_per.append(att_swh_mt[n])
    mean_att.append(np.mean(sum_per))
    mean_std.append(np.nanstd(diff_swh_mt))
    mean_att_meters.append(np.nanmean(diff_swh_mt))

    count=0
    for i in diff_energy_mt:
        if i>=0:
            count+=1
    pos_diff_energy.append((100/len(diff_energy_mt))*count)
    
    sum_per=[]
    count=0
    for n in range(len(att_energy_mt)):
        if not (np.isnan(att_energy_mt[n])) and att_energy_mt[n] > -100:
            sum_per.append(att_energy_mt[n])
    mean_att_energy.append(np.mean(sum_per))
    mean_std_energy.append(np.nanstd(diff_energy_mt))
    mean_att_meters_energy.append(np.nanmean(diff_energy_mt))

ax[1][2].axis('off')
ax_energy[1][2].axis('off')
fig.suptitle('Histogram of SWH differences', y=0.96, fontsize=16)
fig_energy.suptitle('Histogram of energy differences', y=0.96, fontsize=16)

#plots swh
sea_states=['3', '4', '5', '6', '7']
n=len(sea_states)
r=np.arange(n)
fig = plt.figure(figsize = (10, 5))

plt.bar(r -0.1, pos_diff, 0.2, alpha=0.8, facecolor='#773EA9')
plt.bar(r + 0.1, mean_att, 0.2, alpha=0.8, facecolor='#F9D405')

plt.xticks(r,['0.5-1.25 m', '1.25-2.5 m', '2.5-4 m', '4-6 m', '6-9 m'])
plt.xlabel("Sea states")
plt.ylabel("[%]")
plt.legend(['positive differences','mean attenuation'], loc='upper left')
plt.title('Effect of coral reefs on SWH for different sea states')

plt.show()

#plot energy
fig_energy =plt.figure(figsize= (10,5))

plt.bar(r -0.1, pos_diff_energy, 0.2, alpha=0.8, facecolor='#773EA9')
plt.bar(r + 0.1, mean_att_energy, 0.2, alpha=0.8, facecolor='#F9D405')

plt.xticks(r,['0.5-1.25 m', '1.25-2.5 m', '2.5-4 m', '4-6 m', '6-9 m'])
plt.xlabel("Sea states")
plt.ylabel("[%]")
plt.legend(['positive differences','mean attenuation'], loc='upper left')
plt.title('Effect of coral reefs on wave energy for different sea states')
plt.show()

#%% save results
with open('swh_mt_stats.pickle', 'wb') as handle:
    pickle.dump(swh_modeled_tr, handle, protocol=pickle.HIGHEST_PROTOCOL) 
with open('energy_mt_stats.pickle', 'wb') as handle:
    pickle.dump(energy_modeled_tr, handle, protocol=pickle.HIGHEST_PROTOCOL) 
