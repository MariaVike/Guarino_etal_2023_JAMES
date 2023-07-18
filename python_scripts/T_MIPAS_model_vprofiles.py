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


############################################## Function that compute monthly climatological means for each month in dataaset ####################################################################################

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

##### Compute zonal mean of input data at specific latitude ####

def zmean(data_in, lat_in):
  
 #compute zonal mean
 data_in=data_in.collapsed('longitude', iris.analysis.MEAN) 

 if lat_in > 0: #NH
  latitude=iris.Constraint(latitude=lambda lat:   lat_in <= lat <= lat_in+1) 
 else: #SH
  latitude=iris.Constraint(latitude=lambda lat:   lat_in-1 <= lat <= lat_in) 

 data_out=data_in.extract(latitude)

 return(data_out)


#############################################################

def plot_vertical_profile_12_panels(data_in, data_in_stdev, data_ctrl, data_ctrl_stdev, data_obs, data_obs_stdev, p_obs, lat_in, color1, color2, xlabel, label1, label2): 


 fig=plt.figure(tight_layout=True, figsize=(13.5, 6))
 n_r=2
 n_c=6

 month_label=['Jan', 'Feb', 'Mar','Apr', 'May','Jun', 'July', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

 month=np.arange(1,13,1)
 for im in month:

  ax = fig.add_subplot(str(n_r), str(n_c), im) #n_of_row, n_of_columns and plot_id 
  ax=plt.gca()

  data_in_im=data_in[im-1, :20]      #use only levels from about 0 to 2 hPa (here hybrid coord is pure pressure)
  data_in_stdev_im=data_in_stdev[im-1, :20]
  data_ctrl_im=data_ctrl[im-1, :20]
  data_ctrl_stdev_im=data_ctrl_stdev[im-1, :20]
  data_obs_im=data_obs[im-1, 58:]    #use only obs from about 70 km
  data_obs_stdev_im=data_obs_stdev[im-1, 58:]
  p_obs_im=p_obs[im-1, 58:]

  p_lev=(data_in_im.coord('vertical pressure').points)/100
  plt.plot(data_in_im.data, p_lev , c=color1, linewidth=4, linestyle='-', alpha=0.8, label=label1)
  plt.fill_betweenx(p_lev, data_in_im.data+(data_in_stdev_im.data)*2, data_in_im.data-(data_in_stdev_im.data)*2, alpha=0.3, facecolor=color1)

  plt.plot(data_ctrl_im.data, p_lev , c=color2, linewidth=4, linestyle='-', alpha=0.8, label=label2)
  plt.fill_betweenx( p_lev, data_ctrl_im.data+(data_ctrl_stdev_im.data)*2, data_ctrl_im.data-(data_ctrl_stdev_im.data)*2, alpha=0.3, facecolor=color2)
 

  plt.plot(data_obs_im.data, p_obs_im.data, '-o', c='black', linewidth=1, alpha=0.8, label='MIPAS')
  plt.fill_betweenx(p_obs_im.data, data_obs_im.data+(data_obs_stdev_im.data)*2, data_obs_im.data-(data_obs_stdev_im.data)*2, alpha=0.5, facecolor='black')

  ax.set_yscale('log')
  ax.invert_yaxis()
  plt.ylim(10E-3, 10E-5)
  plt.title('{month_label}'.format(month_label=month_label[im-1]), fontsize=12)

  if abs(lat_in) == 80:
   x_lim=(100,300) 
  if abs(lat_in) == 59:
   x_lim=(100,300)
  if abs(lat_in) == 40:
   x_lim=(135,250)
  if abs(lat_in) == 19:
   x_lim=(150,250)

  plt.xlim(x_lim)

  if im == 1 or im == 7: 
   plt.ylabel('Pressure (hPa)', fontsize=12)
  if 7 <= im <= 12:
   plt.xlabel(xlabel, fontsize=12)

 return()


###################################################################################################################
#'''
#Na_Fe runs (CTRL & Kdyn) 
file_path_hist_ctrl='/path_to_data/CTRL_HIST_monthlyclim_20082011_CO2TZ3_monthly.nc'
file_path_hist_Kdyn='/path_to_data/Kdyn_HIST_monthlyclim_20082011_CO2TZ3_monthly.nc'
file_path_MIPAS_T_2='/path_to_data/MIPAS_T_20082011_661_monthly_p.nc'

T_ctrl=iris.load_cube(file_path_hist_ctrl, 'T')
T=iris.load_cube(file_path_hist_Kdyn, 'T')
T_MIPAS=iris.load_cube(file_path_MIPAS_T_2, 'Temperature')
P_MIPAS=iris.load_cube(file_path_MIPAS_T_2, 'pressure')


#lat_in= -80 / 80
#lat_in= -59 / 59
#lat_in= -40 / 40
#lat_in= -19 / 19

#compute zonal means at specified lat for each month
save_fig=True

lat_in=-40
T=zmean(T, lat_in)
T_ctrl=zmean(T_ctrl, lat_in)
T_MIPAS=zmean(T_MIPAS, lat_in)
P_MIPAS=zmean(P_MIPAS, lat_in)

#for high-res sat data there might be multiple datapoints within 1deg of lat
if T_MIPAS.shape[1] >1: 
  T_MIPAS=T_MIPAS[:,0,:]
  P_MIPAS=P_MIPAS[:,0,:]

#print T_MIPAS.data
#compute mean and stdev and return climatology
T, T_stdev=Climatology(T)
T_ctrl, T_ctrl_stdev=Climatology(T_ctrl)
T_MIPAS, T_MIPAS_stdev=Climatology(T_MIPAS)
P_MIPAS, P_MIPAS_stdev=Climatology(P_MIPAS)


#Plot 12 panel figures
plot_vertical_profile_12_panels(T, T_stdev, T_ctrl, T_ctrl_stdev, T_MIPAS, T_MIPAS_stdev, P_MIPAS, lat_in , 'blue', 'red', 'T(K)', 'K_DYN', 'CTRL')
if save_fig:
 if lat_in >0:
  plt.savefig('/path_to_directory/Temp_{lat}N_12months_p.pdf'.format(lat=lat_in))
 else:
  plt.savefig('/path_to_directory/Temp_{lat}S_12months_p.pdf'.format(lat=abs(lat_in)))
plt.show()

plt.show()


