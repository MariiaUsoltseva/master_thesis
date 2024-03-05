#%%
import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import griddata
import time
from common_settings import lonmin, lonmax, latmin, latmax
from common_settings import satellites, sat_flags
from common_settings import months, years
from poly_funcs import sort_by_pass, filter_by_flags_seastate, mean_of_flags_seastate_new
#%%
#documnent to analyse SWH differences on regional scale for different coral cover
#variables
satellite_numbers_2018=[3,5]
satellite_numbers_2009=[2,6]
satellite_numbers_2003=[1,6] 
years=['2018', '2009', '2003']
ssl =0.5
ssh=100
sat_numbers=[satellite_numbers_2018, satellite_numbers_2009, satellite_numbers_2003]

with open('sat_1_flags_final.pickle', 'rb') as handle:
    sat_1_flags = pickle.load(handle)
with open('sat_2_flags_final.pickle', 'rb') as handle:
    sat_2_flags = pickle.load(handle)
with open('sat_3_flags_final.pickle', 'rb') as handle:
    sat_3_flags = pickle.load(handle)
with open('sat_5_flags_final.pickle', 'rb') as handle:
    sat_5_flags = pickle.load(handle)
with open('sat_6_flags_final.pickle', 'rb') as handle:
    sat_6_flags = pickle.load(handle)
pass_flags_sats_2018=[sat_3_flags, sat_5_flags]
pass_flags_sats_2009=[sat_2_flags, sat_6_flags]
pass_flags_sats_2003=[sat_1_flags, sat_6_flags]
pass_flags=[pass_flags_sats_2018, pass_flags_sats_2009, pass_flags_sats_2003]
#%%
start_time=time.time()
"""
Coordinates of reef areas
#east asian seas : 90.7, 141.5, -10,  26.8
#pacific: 150, 180, -30,  30 and -180, -120, -30,  30
#australia: 108, 158, -30,  -10
#carribean: -89.2,-55.65, 9.96, 29.26
#west indian ocean : 34.9, 68.4, -25.8, 1.9
#red sea 30.1, 55.9, 0.7, 30
#south asia: 66, 87.2, -7, 24.6
#east-tropical pacific: -119.9, -75.7, -2, 25
"""
months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

for i, year in enumerate(years):
    print('Processing year', year)
    satellite_numbers=sat_numbers[i]
    pass_flags_sats=pass_flags[i]
    differences_all=[]
    att_all=[]
    for month in months:
        month_folder='D:/thesis/ceda_data/{}/{}'.format(year, month)
        for i in range(len(satellite_numbers)):
            pass_nums, pass_data_dict = sort_by_pass(month_folder, satellite_numbers[i],90.7, 141.5, -10, 26.8)
            
            ''' buffer_filtered, reef_filtered, difference_swh = filter_by_flags_seastate(pass_data_dict, pass_flags_sats, pass_nums, i,ssl=1, ssh=100)
            for num in difference_swh:
               diff = difference_swh[num]['swh']
               att = difference_swh[num]['attenuation_percent']
               differences_all.extend(diff)
               att_all.extend(att)
            
    plt.hist(differences_all, bins=50)
    plt.xlabel('difference [m]')
    plt.ylabel('Frequency')
    plt.title('Histogram of SWH differences ')
    plt.grid(True)
    plt.xlim(right=2.5) 
    plt.xlim(left=-2.5)
    plt.show()

    count=0
    positive=[]
    for n,i in enumerate(differences_all):
        if i>=0:
            count=count+1
            positive.append(att_all[n])
    percents=(100/len(differences_all))*count #percentage of positive differences
    mean_alt_m=np.nanmean(differences_all) # mean difference in meters
    mean_att_pos=np.nanmean(positive) # mean positive attenuation in percents
    mean_att=np.nanmean(att_all) # mean  attenuation in percents


    #mean_positive=sum(positive)/count
    print("Positive differences ", percents, "%")
    print("mean of positive att", mean_att_pos, "%")
    print("Mean attenuation ", mean_att, "%")

end_time=time.time()
print("Execution time = ", (end_time-start_time)/60 , 'minutes')'''
#%%