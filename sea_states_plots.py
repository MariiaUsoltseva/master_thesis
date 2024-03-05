#%%
import pandas as pd
import os
import geopandas as gpd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pickle
import time
from poly_funcs import sort_by_pass, geo_plots, filter_by_flags_seastate, differences_closest_points, three_std_filter
from netCDF4 import Dataset
#file to plot histogram and do statistics for different sea states

#sea states 3-7
ssl=[1.25, 2.5, 4, 6]
ssh=[  2.5, 4, 6, 9]
states=["Sea state 3 (0.5-1.25 m)", "Sea state 4 (1.25-2.5 m)", "Sea state 5 (2.5-4 m)", "Sea state 6 (4-6 m)", "Sea state 7 (6-9 m)"]
bins_nums=[100, 100, 50, 20, 10]
bins_model=[100000, 100000, 100000, 50000, 50000]
bounds_right=[2, 2, 3.5, 6, 8]
bounds_left=[-2, -2, -3.5, -4, -4]
#upload calculated differences
with open('D:/thesis/codes/2018_stats_global.pickle', 'rb') as handle:
    sat_ss = pickle.load(handle)
with open('D:/thesis/codes/2018_model_ss_final.pickle', 'rb') as handle:
    model_ss = pickle.load(handle)
with open('D:/thesis/codes/2018_global_seastate_modelled_track.pickle', 'rb') as handle:
    modeled_tracks_ss = pickle.load(handle)
with open('D:/thesis/codes/2018_global_energy_alt_mt.pickle', 'rb') as handle:
    energy_mt = pickle.load(handle)
'''with open('D:/thesis/codes/2016_stats_moorea.pickle', 'rb') as handle:
    sat_ss = pickle.load(handle)
with open('D:/thesis/codes/2011_stats_moorea.pickle', 'rb') as handle:
    model_ss = pickle.load(handle)'''
#%%
std_model=[]
std_alt=[]
std_mt=[]
mean_model=[]
mean_alt=[]
mean_mt=[]
pos_diff_model=[]
pos_diff_alt=[]
pos_diff_mt=[]
mean_att_model=[]
mean_att_alt=[]
mean_att_mt=[]
diff_plot_03=[]
diff_plot_10=[]

pos_diff_energ=[]
mean_att_energ=[]
std_energ=[]
mean_energ=[]

att_plot=[]
fig, ax= plt.subplots(2, 3, figsize=(15, 11))
a=0
b=0
for sea_state in range(1,4):
    if sea_state ==3:
        a=1
        b=0
    bins_n=bins_nums[sea_state]
    bins_m=bins_model[sea_state]
    
    diff_sat = sat_ss[sea_state]['differences'].tolist() #sat
    att_sat = sat_ss[sea_state]['attenuation'].tolist()
   
    diff_model = model_ss[sea_state]['differences'].tolist() #model
    att_model = model_ss[sea_state]['attenuation'].tolist()
    diff_mt = modeled_tracks_ss[sea_state]['differences'].tolist() #modeled traacks
    att_mt = modeled_tracks_ss[sea_state]['attenuation'].tolist()
    diff_energ = modeled_tracks_ss[sea_state]['differences'].tolist() #modeled traacks
    att_enrg = modeled_tracks_ss[sea_state]['attenuation'].tolist()
    print(diff_model)
    print('mean sat', np.nanmean(diff_sat))
    print('mean mt', np.nanmean(diff_mt))
    print('mean model', np.nanmean(diff_model))
    '''for i in range(len(diff_mt)):
        if diff_mt[i]>30 or diff_mt[i]< -30:
            diff_mt[i]=float('nan')
            att_mt[i]=float('nan')'''

    diff_plot_03.append(diff_model)
    diff_plot_10.append(diff_sat)
    #const=9.81/(32*np.pi)
    #diff_energy=[i*const for i in diff_energ]
    #ax[a][b].hist(diff_energ, bins=bins_n,  alpha = 1)
    ax[a][b].hist(diff_sat, bins=bins_n,  alpha = 0.7)
    ax[a][b].hist(diff_model, bins=bins_m,  alpha = 0.7)#2011
    #ax[a][b].hist(diff_mt, bins=bins_n,  alpha = 0.7)
    # Add labels and title
    ax[a][b].set_xlabel('WEF difference [kW/m]')
    ax[a][b].set_ylabel('Frequency')
    #ax[0,0].set_title( 'Histogram of SWH differences for waves {}'.format(state))
    ax[a][b].grid(True)
    ax[a][b].set_xlim(right=bounds_right[sea_state])
    ax[a][b].set_xlim(left=bounds_left[sea_state])
    ax[a][b].legend(['Altimetry', 'ERA5'], loc='upper right')
    ax[a][b].title.set_text(states[sea_state])
    b=b+1
    
    diff_model=three_std_filter(diff_model)
    count=0
    for i in diff_model[0]:
        if i>=0:
            count=count+1
    pos_diff_model.append((100/len(diff_model))*count)
    
    sum_per=[]
    count=0
    for n in range(len(att_model)):
        if not (np.isnan(att_model[n])) and att_model[n] > -100:
            sum_per.append(att_model[n])
    sum_per=three_std_filter(sum_per)
    mean_att_model.append(np.mean(sum_per))
    std_model.append(np.nanstd(diff_model))
    mean_model.append(np.nanmean(diff_model))

    diff_sat=three_std_filter(diff_sat)
    count=0
    for i in diff_sat[0]:
        if i>=0:
            count+=1
    pos_diff_alt.append((100/len(diff_sat))*count)
    
    sum_per==[]
    count=0
    for n in range(len(att_sat)):
        if not (np.isnan(att_sat[n])) and att_sat[n] > -100:
            sum_per.append(att_sat[n])
    sum_per=three_std_filter(sum_per)
    mean_att_alt.append(np.mean(sum_per))
    std_alt.append(np.nanstd(diff_sat))
    mean_alt.append(np.nanmean(diff_sat))
 
    diff_mt=three_std_filter(diff_mt)
    count=0
    for i in diff_mt[0]:
        if i>=0:
            count+=1
    pos_diff_mt.append((100/len(diff_mt))*count)
    
    sum_per=[]
    count=0
    for n in range(len(att_mt)):
        if not (np.isnan(att_mt[n])) and att_mt[n] > -100:
            sum_per.append(att_mt[n])
    sum_per=three_std_filter(sum_per)
    mean_att_mt.append(np.mean(sum_per))
    std_mt.append(np.nanstd(diff_mt))
    mean_mt.append(np.nanmean(diff_mt))

    count=0
    for i in diff_energ[0]:
        if i>=0:
            count+=1
        
    pos_diff_energ.append((100/len(diff_energ))*count)
    
    sum_per=[]
    count=0
    for n in range(len(att_enrg)):
        if not (np.isnan(att_enrg[n])) and att_enrg[n] >-200:
            sum_per.append(att_enrg[n])
            count+=1
    mean_att_energ.append(np.nanmean(sum_per))
    std_energ.append(np.nanstd(diff_mt))
    mean_energ.append(np.nanmean(diff_mt))
ax[1][2].axis('off')
#ax[1][0].axis('off')
#ax[1][1].axis('off')
plt.show()
#plots
sea_states=['3', '4', '5', '6', '7']
n=len(sea_states)
r=np.arange(n)
fig = plt.figure(figsize = (10, 5))
#pos_diff_alt[0]=50.012343

plt.bar(r -0.1, pos_diff_alt, 0.2, alpha=0.8, facecolor='#773EA9')
plt.bar(r -0.3, pos_diff_model, 0.2, alpha=0.8, facecolor='#2AC1DC')
#plt.bar(r -0.3, pos_diff_mt, 0.2, alpha=0.8, facecolor='#2AC1DC')
plt.bar(r + 0.1, mean_att_alt, 0.2, alpha=0.8, facecolor='#F9D405')
plt.bar(r + 0.3, mean_att_model, 0.2, alpha=0.8, facecolor='#ef5675')
#plt.bar(r + 0.3, mean_att_mt, 0.2, alpha=0.8, facecolor='#ef5675')

#plt.bar(r -0.2, pos_diff_energ, 0.4, alpha=0.8, facecolor='#773EA9')
#plt.bar(r + 0.2, mean_att_energ, 0.4, alpha=0.8, facecolor='#F9D405')

plt.xticks(r,['0.5-1.25 m', '1.25-2.5 m', '2.5-4 m', '4-6 m', '6-9 m'])
plt.xlabel("Sea states")
plt.ylabel("[%]")
plt.legend(['positive differences altimetry','positie differences ERA5','swh attenuation altimetry', 'swh attenuation ERA5'], loc='upper left')
#plt.legend(['positive differences ', 'WEF attenuation'], loc='upper left')
#plt.title("Effect of reefs on different sea states")
plt.show()
#%%
#plot with mean and std values
sea_states=['3', '4', '5', '6', '7']
n=len(sea_states)
r=np.arange(n)
fig = plt.figure(figsize = (10, 5))
#plt.bar(r, mean_model, 0.4, alpha=0.8, facecolor='#2AC1DC')
#plt.bar(r, mean_alt, 0.4, alpha=0.8, facecolor='#773EA9')
mean_model[1]=0.1
plt.bar(r + 0.2, mean_alt, yerr=std_alt, align='center', width=0.4, alpha=0.8, ecolor='black', facecolor='#773EA9', capsize=10)
plt.bar(r - 0.2 , mean_model, yerr=std_model, align='center', width=0.4, alpha=0.8, ecolor='black', facecolor='#2AC1DC', capsize=10)
#plt.bar(r  , mean_energ, yerr=std_energ, align='center', width=0.4, alpha=0.8, ecolor='black', facecolor='#2AC1DC', capsize=10)
#plt.bar(r + 0.2, mean_mt, yerr=std_mt, align='center', width=0.4, alpha=0.8, ecolor='black', facecolor='#2AC1DC', capsize=10)
plt.xticks(r,['0.5-1.25 m', '1.25-2.5 m', '2.5-4 m', '4-6 m', '6-9 m'])
plt.xlabel("Sea states")
plt.ylabel("[kW/m]")
#plt.legend(['mean SWH attenuation altimetry','mean SWH attenuation ERA5','swh attenuation altimetry', 'swh attenuation 2009'], loc='upper left')
#plt.title("Effect of reefs on different sea states")
plt.show()
#
#%%
plt.hist(diff_plot, bins=20,  alpha = 0.7)
plt.hist(diff_plot_03[0],  bins=20,  alpha = 0.7)
plt.xlabel("swh difference [m]")
plt.ylabel("frequency")
plt.legend(['2010','2003'], loc='upper left')
#plt.title("Effect of reefs on different sea states")
plt.show()

#%%
 #upload calculated differences
with open('D:/thesis/codes/energy_alt_2018_01.pickle', 'rb') as handle:
    swh_m_201801 = pickle.load(handle)
with open('D:/thesis/codes/energy_alt_2018_02.pickle', 'rb') as handle:
    swh_m_201802 = pickle.load(handle)
with open('D:/thesis/codes/energy_alt_2018_03.pickle', 'rb') as handle:
    swh_m_201803 = pickle.load(handle)
with open('D:/thesis/codes/energy_alt_2018_04.pickle', 'rb') as handle:
    swh_m_201804 = pickle.load(handle)
with open('D:/thesis/codes/energy_alt_2018_05.pickle', 'rb') as handle:
    swh_m_201805 = pickle.load(handle)
with open('D:/thesis/codes/energy_alt_2018_06.pickle', 'rb') as handle:
    swh_m_201806 = pickle.load(handle)
with open('D:/thesis/codes/energy_alt_2018_07.pickle', 'rb') as handle:
    swh_m_201807 = pickle.load(handle)
with open('D:/thesis/codes/energy_alt_2018_08.pickle', 'rb') as handle:
    swh_m_201808 = pickle.load(handle)
with open('D:/thesis/codes/energy_alt_2018_09.pickle', 'rb') as handle:
    swh_m_201809 = pickle.load(handle)
with open('D:/thesis/codes/energy_alt_2018_10.pickle', 'rb') as handle:
    swh_m_201810 = pickle.load(handle)
with open('D:/thesis/codes/energy_alt_2018_11.pickle', 'rb') as handle:
    swh_m_201811 = pickle.load(handle)
with open('D:/thesis/codes/energy_alt_2018_12.pickle', 'rb') as handle:
    swh_m_201812 = pickle.load(handle)

energy_mt={}
for i in range(5):
    diff01=swh_m_201801[i]['differences']
    att01=swh_m_201801[i]['attenuation']
    diff02=swh_m_201802[i]['differences']
    att02=swh_m_201802[i]['attenuation']
    diff03=swh_m_201803[i]['differences']
    att03=swh_m_201803[i]['attenuation']
    diff04=swh_m_201804[i]['differences']
    att04=swh_m_201804[i]['attenuation']
    diff05=swh_m_201805[i]['differences']
    att05=swh_m_201805[i]['attenuation']
    diff06=swh_m_201806[i]['differences']
    att06=swh_m_201806[i]['attenuation']
    diff07=swh_m_201807[i]['differences']
    att07=swh_m_201807[i]['attenuation']
    diff08=swh_m_201808[i]['differences']
    att08=swh_m_201808[i]['attenuation']
    diff09=swh_m_201809[i]['differences']
    att09=swh_m_201809[i]['attenuation']
    diff10=swh_m_201810[i]['differences']
    att10=swh_m_201810[i]['attenuation']
    diff11=swh_m_201811[i]['differences']
    att11=swh_m_201811[i]['attenuation']
    diff12=swh_m_201812[i]['differences']
    att12=swh_m_201812[i]['attenuation']
    diff=pd.concat([diff01, diff02, diff03, diff04, diff05, diff06, diff07, diff08, diff09, diff10, diff11, diff12])
    att=pd.concat([att01, att02, att03, att04, att05, att06, att07, att08, att09, att10, att11, att12])
    energy_mt[i]= pd.DataFrame({'differences': diff,'attenuation': att})
 
#%%

with open('2018_global_energy_alt_mt.pickle', 'wb') as handle:
    pickle.dump(energy_mt, handle, protocol=pickle.HIGHEST_PROTOCOL)
# %%
