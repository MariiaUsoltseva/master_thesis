## Description

Codes for the master's thesis titled "Effect of Coral Reefs on Wave Energy and Wave Height"

## Codes included

• common settings.py. - the file includes the year, region of interest and paths to
folders with satellite and model data and reef extends.

• satellite {IDnumber} flags.pickle- the file contains location flags calculated for each
ground track of satellite {IDnumber}. Satellite ID numbers are defined in the
Sea State CCI dataset: 0-Cryosat, 1-Jason1, 2-Jason2, 3-Jason3, 4-Saral, 5-Sentinel-3A, 6-Envisat.

• altimetry functions.py, model functions.py - the files include all functions used for
data processing from altimetry and model datasets.

• global maps.py - the file used to generate maps with mean SWH and wave energy flux
on coral reefs, right before coral reefs and with wave height and energy attenuation.

• statistics CCI.py, statistics ERA5.py, statistics ERA5 ontracks.py - calculates mean
SWH attenuation, attenuation percentage and percentage of positive differences for

each sea state ranging from 0.5 m to 9 m utilizing data from SeaState CCI dataset,
ERA5 on a grid and ERA5 interpolated on altimetry tracks respectively.

• coral cover.py - calculates mean SWH attenuation, attenuation percentage and per-
centage of positive differences for two years with different coral reefs structural
complexity for the region of interest.

## Installation

To create a python environment to run the codes use file environment.yml
