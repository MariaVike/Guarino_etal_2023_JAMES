import iris 
import iris.plot as iplt
import iris.quickplot as qplt

import iris.analysis.cartography 
import cartopy.crs as ccrs
import cartopy.feature as cfe

import os
from mpl_toolkits.basemap import Basemap, maskoceans, shiftgrid

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from iris.experimental.equalise_cubes import *
from iris.experimental.animate import*
import matplotlib.dates as mdates
import matplotlib.ticker as mtick
from matplotlib import ticker

import warnings
import glob
import os 

import pandas as pd
import scipy

from shiftedColorMap import*
from scipy import stats

from iris.cube import Cube

########################################### Function that compute monthly climatological means for each month in dataaset ####################################################################################

def Climatology (data_in):

 #compute long-term mean for each month
 month=data_in.shape[0]

 months=iris.cube.CubeList()
 stdev=iris.cube.CubeList()
 for i in range (0, 12):
  months.append(data_in[i:month:12].collapsed('time', iris.analysis.MEAN))  
  stdev.append(data_in[i:month:12].collapsed('time', iris.analysis.STD_DEV))  


 months=months.merge_cube()
 stdev=stdev.merge_cube()

 return(months, stdev)

###############################################  Extract selected latitudinal band or single lat (and compute zonal mean) #################################################################################################

def extract_lat_band(data_in, lat1, lat2):

 R=iris.Constraint(latitude=lambda lat:  lat1 <= lat <= lat2) 
 data_out=data_in.extract(R)

 #data_out=data_out.collapsed('time', iris.analysis.MEAN) 
 #data_zonal=data_out.collapsed('longitude', iris.analysis.MEAN)

 return(data_out) #data_zonal

def extract_lat(data_in, lat_in, SH):

 if SH:
  R=iris.Constraint(latitude=lambda lat:   lat_in-1 <= lat <= lat_in) 
 else:
  R=iris.Constraint(latitude=lambda lat:   lat_in <= lat <= lat_in+1) 

 data_out=data_in.extract(R)

 #data_out=data_out.collapsed('time', iris.analysis.MEAN) 
 #data_zonal=data_out.collapsed('longitude', iris.analysis.MEAN)

 return(data_out) #data_zonal

###############################################  Compute area weighed means  #################################################################################################

def area_weighted(data_in):

 data_in.coord('latitude').guess_bounds()
 data_in.coord('longitude').guess_bounds()
 weights_in = iris.analysis.cartography.area_weights(data_in)

 weighted_mean= data_in.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, 
					weights=weights_in)
 data_max= data_in.collapsed(['latitude', 'longitude'], iris.analysis.MAX)
 data_min= data_in.collapsed(['latitude', 'longitude'], iris.analysis.MIN)
 data_stdev= data_in.collapsed(['latitude', 'longitude'], iris.analysis.STD_DEV)

 return(weighted_mean, data_max, data_min, data_stdev)

############################## compute area means, non spatially weighed (for satellite obs no needed) ##############################
def area_mean(data_in):

 data_in.coord('latitude').guess_bounds()
 data_in.coord('longitude').guess_bounds()

 data_mean= data_in.collapsed(['latitude', 'longitude'], iris.analysis.MEAN)
 data_stdev= data_in.collapsed(['latitude', 'longitude'], iris.analysis.STD_DEV)

 return(data_mean, data_stdev)

########## Interpolate model data onto the observation vertical grid (pressure levels) ##############################

def to_pressure_obs(data_in, P_obs):

 #print data_in.coord('atmosphere_hybrid_sigma_pressure_coordinate')
 #print data_in.coord('vertical pressure')
 ref_p = [('vertical pressure', P_obs.data*100)] #convert obs from hPa to Pa
 interp_cube = data_in.interpolate(ref_p, iris.analysis.Linear(extrapolation_mode='mask'))

 return(interp_cube)

############### Plot vertical profiles ################################################################################################################

def plot_vertical_profile_obs_latband(data_in, data_ctrl, data_obs, P_obs, lat1, lat2, month, color1, color2, xlabel, label1, label2, title): 

 #extract latitudinal band for specified month
 data_ctrl=extract_lat_band(data_ctrl[month-1], lat1, lat2)
 data_in=extract_lat_band(data_in[month-1], lat1, lat2)


 #compute area weighted mean, stdev and max and min for each lat bands:
 data_ctrl, data_ctrl_max, data_ctrl_min, data_ctrl_stdev =area_weighted(data_ctrl)
 data_in, data_in_max, data_in_min, data_in_stdev=area_weighted(data_in)

 #do as above but for obs
 data_obs=extract_lat_band(data_obs[month-1,:,:,58:], lat1, lat2) #use only obs from about 70 km
 P_obs=extract_lat_band(P_obs[month-1,:,:,58:], lat1, lat2)
 data_obs, data_obs_stdev=area_mean(data_obs)
 P_obs, P_obs_stdev=area_mean(P_obs)

 #plot data
 data_in=data_in[:20] #use only levels from about 0 to 2 hPa (here hybrid coord is pure pressure)
 data_in_stdev=data_in_stdev[:20]
 data_ctrl=data_ctrl[:20]
 data_ctrl_stdev=data_ctrl_stdev[:20] 
 
 p_lev=(data_in[:20].coord('vertical pressure').points)/100

 fig=plt.figure(tight_layout=True, figsize=(5, 6))
 ax=plt.gca()
 plt.plot(data_in.data, p_lev , c=color1, linewidth=4, linestyle='-', alpha=0.8, label=label1)
 plt.fill_betweenx(p_lev, data_in.data+(data_in_stdev.data)*2, data_in.data-(data_in_stdev.data)*2, alpha=0.3, facecolor=color1)

 p_lev=(data_ctrl[:20].coord('vertical pressure').points)/100
 plt.plot(data_ctrl.data, p_lev , c=color2, linewidth=4, linestyle='-', alpha=0.8, label=label2)
 plt.fill_betweenx( p_lev, data_ctrl.data+(data_ctrl_stdev.data)*2, data_ctrl.data-(data_ctrl_stdev.data)*2, alpha=0.3, facecolor=color2)

 plt.plot(data_obs.data, P_obs.data, '-o', c='black', linewidth=1, alpha=0.8, label='MIPAS')
 plt.fill_betweenx(P_obs.data, data_obs.data+(data_obs_stdev.data)*2, data_obs.data-(data_obs_stdev.data)*2, alpha=0.5, facecolor='black')

 ax.set_yscale('log')
 ax.invert_yaxis()
 plt.ylim(10E-3, 10E-5)
 plt.xlabel(xlabel, fontsize=14)
 plt.ylabel('Pressure (hPa)', fontsize=14)
 plt.title(title, fontsize=16)
 plt.legend()
 plt.grid()

 return()

############################################################### #Load datasets and call functions ######################################################

####### Plot CO2 within lat bands ############################################################
file_path_hist_ctrl='/path_to_data/CTRL_HIST_monthlyclim_20082011_CO2TZ3_monthly.nc'
file_path_hist_Kdyn='/path_to_data/Kdyn_HIST_monthlyclim_20082011_CO2TZ3_monthly.nc'

file_path_MIPAS_CO2='/path_to_data/MIPAS_CO2_20082011_661_monthly_p.nc'


CO2_ctrl=iris.load_cube(file_path_hist_ctrl, 'CO2')
CO2=iris.load_cube(file_path_hist_Kdyn, 'CO2')

CO2_MIPAS=iris.load_cube(file_path_MIPAS_CO2, 'CO2')
P_MIPAS=iris.load_cube(file_path_MIPAS_CO2, 'pressure')

#compute mean (and stdev) and return climatology
CO2, CO2_stdev=Climatology(CO2)
CO2_ctrl, CO2_ctrl_stdev=Climatology(CO2_ctrl)
CO2_MIPAS, CO2_MIPAS_stdev=Climatology(CO2_MIPAS)
P_MIPAS, P_MIPAS_stdev=Climatology(P_MIPAS)


save_fig=False

#June NH
month=6 
lat1=-30
lat2=30
plot_vertical_profile_obs_latband(CO2*1E6, CO2_ctrl*1E6, CO2_MIPAS, P_MIPAS, lat1, lat2, month, 'blue', 'red', 'ppm', 'K_DYN', 'CTRL', title='CO2 month:{month}, Lat: [{lat1}, {lat2}]'.format(month=month, lat1=lat1, lat2=lat2))
if save_fig:
 plt.savefig('CO2_latband_{lat1}S{lat2}N_month_{month}_p.jpg'.format(lat1=abs(lat1), lat2=lat2, month=month))

lat1=30
lat2=60
plot_vertical_profile_obs_latband(CO2*1E6, CO2_ctrl*1E6, CO2_MIPAS, P_MIPAS, lat1, lat2, month, 'blue', 'red', 'ppm', 'K_DYN', 'CTRL', title='CO2 month:{month}, Lat: [{lat1}, {lat2}]'.format(month=month, lat1=lat1, lat2=lat2))
if save_fig:
 plt.savefig('CO2_latband_{lat1}N{lat2}N_month_{month}_p.jpg'.format(lat1=lat1, lat2=lat2, month=month))

lat1=60
lat2=90
plot_vertical_profile_obs_latband(CO2*1E6, CO2_ctrl*1E6, CO2_MIPAS, P_MIPAS, lat1, lat2, month, 'blue', 'red', 'ppm', 'K_DYN', 'CTRL', title='CO2 month:{month}, Lat: [{lat1}, {lat2}]'.format(month=month, lat1=lat1, lat2=lat2))
if save_fig:
 plt.savefig('CO2_latband_{lat1}N{lat2}N_month_{month}_p.jpg'.format(lat1=lat1, lat2=lat2, month=month))

plt.show()



