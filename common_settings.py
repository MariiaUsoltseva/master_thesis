from datetime import datetime
import pickle

#Borders of significant reef regions (lonmin, lonmax, latmin, latmax)
#east asian seas : 90.7, 141.5, -10,  26.8
#pacific: 150, 180, -30,  30 and -180, -120, -30,  30
#australia: 108, 158, -30,  -10
#carribean: -89.2,-55.65, 9.96, 29.26
#west indian ocean : 34.9, 68.4, -25.8, 1.9
#red sea 30.1, 55.9, 0.7, 30
#south asia:66, 87.2, -7, 24.6
#et pacific: -119.9, -75.7, -2, 25

#define geographic area
lonmin = -180
lonmax = 180
latmin = -30
latmax = 30

#define timeframe 
start_date = datetime(2018, 4, 1) #year, month, day
end_date = datetime(2018, 4, 1)

#list of years and months for further processing
all_months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
#years included in sea state cci
all_years = ['2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010', 
             '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018', '2019', '2020']
start_month = int(start_date.month)
end_month = int(end_date.month)
start_year = int(start_date.year)-2002
end_year = int(end_date.year)-2001

months = all_months[(start_month-1):(end_month)]
years = all_years[(start_year):(end_year)]

#find stellites operating during the chosen timeframe
#satellite id as defined in SeaState CCI 0-cryosat, 1-jason1, 2-jason2, 
#3-jason3, 4-saral, 5- sentinel3a, 6-envisat
jason1_start = datetime(2002, 1, 15)
jason1_end = datetime(2009, 1, 26)

jason2_start = datetime(2008, 6, 4)
jason2_end = datetime(2016, 10, 17)

jason3_start = datetime(2016, 2, 17)
jason3_end = datetime(2020, 12, 31)
 
saral_start = datetime(2013, 3, 14)
saral_end = datetime(2016, 6, 4)

sentinel3a_start = datetime(2016, 6, 1)
sentinel3a_end = datetime(2020, 12, 31)

envisat_start = datetime(2002, 5, 14)
envisat_end = datetime(2010, 10, 22)

satellites = []
sat_flags = []
if (jason1_start <= start_date <= jason1_end) | (jason1_start <= end_date <= jason1_end):
    satellites.append(1)
    #append flags
    with open('sat_1_1cycle.pickle', 'rb') as handle:
        sat_1_flags = pickle.load(handle)
    sat_flags.append(sat_1_flags)

if (jason2_start <= start_date <= jason2_end) | (jason2_start <= end_date <= jason2_end):
    satellites.append(2)
    #append flags
    with open('sat_2_flags_final.pickle', 'rb') as handle:
        sat_2_flags = pickle.load(handle)
    sat_flags.append(sat_2_flags)

if (jason3_start <= start_date <= jason3_end) | (jason3_start <= end_date <= jason3_end):
    satellites.append(3)
    #append flags
    with open('sat_3_flags_final.pickle', 'rb') as handle:
        sat_3_flags = pickle.load(handle)
    sat_flags.append(sat_3_flags)

if (saral_start <= start_date <= saral_end) | (saral_start <= end_date <= saral_end):
    satellites.append(4)
    #append flags
    with open('sat_4_flags_final.pickle', 'rb') as handle:
        sat_4_flags = pickle.load(handle)
    sat_flags.append(sat_4_flags)

if (sentinel3a_start <= start_date <= sentinel3a_end) | (sentinel3a_start <= end_date <= sentinel3a_end):
    satellites.append(5)
    #append flags
    with open('sat_5_flags_final.pickle', 'rb') as handle:
        sat_5_flags = pickle.load(handle)
    sat_flags.append(sat_5_flags)

if (envisat_start <= start_date <= envisat_end) | (envisat_start <= end_date <= envisat_end):
    satellites.append(6)
    #append flags
    with open('sat_6_flags_final.pickle', 'rb') as handle:
        sat_6_flags = pickle.load(handle)
    sat_flags.append(sat_6_flags)

information = 'Processing for year(s) {}, month(s) {}, in area defined as lonmin={}, lonmax={}, latmin={}, latmax={}. To change the timeframe of geographical boundaries go to file common_settings.py'.format(years,months, lonmin, lonmax, latmin, latmax)
