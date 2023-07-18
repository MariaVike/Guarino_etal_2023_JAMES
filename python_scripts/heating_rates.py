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

from collections import OrderedDict


#### Different linestyles ######

linestyles_dict = OrderedDict(
    [('solid',               (0, ())),
     ('loosely dotted',      (0, (1, 10))),
     ('dotted',              (0, (1, 5))),
     ('densely dotted',      (0, (1, 1))),

     ('loosely dashed',      (0, (5, 10))),
     ('dashed',              (0, (5, 5))),
     ('densely dashed',      (0, (5, 1))),

     ('loosely dashdotted',  (0, (3, 10, 1, 10))),
     ('dashdotted',          (0, (3, 5, 1, 5))),
     ('densely dashdotted',  (0, (3, 1, 1, 1))),

     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))])

###############################################  Extract selected latitudinal band #################################################################################################

def extract_lat_band(data_in, lat1, lat2):

 R=iris.Constraint(latitude=lambda lat:  lat1 <= lat <= lat2) 
 data_out=data_in.extract(R)

 return(data_out) 

###############################################  Compute area weighted means  #################################################################################################

def area_weighted(data_in):

 data_in.coord('latitude').guess_bounds()
 data_in.coord('longitude').guess_bounds()
 weights_in = iris.analysis.cartography.area_weights(data_in)

 weighted_mean= data_in.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, 
								weights=weights_in)

 return(weighted_mean)

############## Plot vertical profile within lat band ################################################################################################################

def vertical_profile_mltp_lat_band(D1, GPH, month, xmin, xmax, color, style, label, xlabel, title, single, lat1, lat2): 


 if single:
  #Extract data within lat band 
  D1=extract_lat_band(D1, lat1, lat2)
  GPH=extract_lat_band(GPH, lat1, lat2)
  D1=area_weighted(D1)
  GPH=area_weighted(GPH)
  D1_data=D1[month-1,:70].data
  GPH_data=GPH[month-1,:70].data/1000

  fig=plt.figure(tight_layout=True, figsize=(5, 6))
  ax=plt.gca()
  plt.plot(D1_data, GPH_data, c=color, linewidth=4, linestyle='-', alpha=0.8, label=label)
  
 else:
  #Extract data within lat band 
  D1_data=iris.cube.CubeList()
  for ic in range(0, len(D1)):
   D1_lats=extract_lat_band(D1[ic], lat1, lat2)
   D1_mean=area_weighted(D1_lats)
   D1_data.append(D1_mean[month-1,:70].data)
  GPH=extract_lat_band(GPH, lat1, lat2)
  GPH=area_weighted(GPH)
  GPH_data=GPH[month-1,:70].data/1000
  #sigma_l=D1[0][0,:70].coord('atmosphere_hybrid_sigma_pressure_coordinate').points

  fig=plt.figure(tight_layout=True, figsize=(5, 6))
  ax=plt.gca()
  plt.axvline(x=0, c='gray', linestyle='dashdot')  
  for i in range(0, len(D1_data)):
   plt.plot(D1_data[i], GPH_data, c=color[i], linewidth=4, linestyle=style[i], alpha=0.8, label=label[i])
 

 #fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
 #fmt.set_powerlimits((0, 0))

 plt.ylim(70,100)
 plt.xlim(xmin, xmax)
 plt.xlabel(xlabel, fontsize=14)
 plt.ylabel('GPH (Km)', fontsize=14)
 #ax.xaxis.set_major_formatter(fmt)
 plt.title(title, fontsize=16)
 plt.legend()
 plt.grid()



############################################################### #Load datasets and call functions ######################################################
#Na_Fe runs (MIF/2 CTRL & Kdyn) 
file_path_hist_ctrl='/path_to_data/CTRL_HIST_monthlyclim_20082011.nc'
file_path_hist_Kdyn='/path_to_data/Kdyn_HIST_monthlyclim_20082011.nc'
file_path_hist_ctrl_heating='/path_to_data/CTRL_HIST_monthlyclim_20082011_heating.nc'
file_path_hist_Kdyn_heating='/path_to_data/Kdyn_HIST_monthlyclim_20082011_heating.nc'
file_path_hist_ctrl_GWdrag='/path_to_data/CTRL_HIST_monthlyclim_20082011_GWdrag.nc'
file_path_hist_Kdyn_GWdrag='/path_to_data/Kdyn_HIST_monthlyclim_20082011_GWdrag.nc'
file_path_hist_ctrl_GWdrag_2='/path_to_data/CTRL_HIST_monthlyclim_20082011_GWdrag_2.nc'
file_path_hist_Kdyn_GWdrag_2='/path_to_data/Kdyn_HIST_monthlyclim_20082011_GWdrag_2.nc'

QRS_TOT=iris.load_cube(file_path_hist_Kdyn_heating, 'QRS_TOT') 
QCP=iris.load_cube(file_path_hist_Kdyn_heating, 'QCP')
QRS_EUV=iris.load_cube(file_path_hist_Kdyn_heating, 'QRS_EUV')
QRS_CO2NIR=iris.load_cube(file_path_hist_Kdyn_heating, 'QRS_CO2NIR')
QRS_AUR=iris.load_cube(file_path_hist_Kdyn_heating, 'QRS_AUR')
QTHERMAL=iris.load_cube(file_path_hist_Kdyn_heating, 'QTHERMAL'))
QRL_TOT=iris.load_cube(file_path_hist_Kdyn_heating, 'QRL_TOT') #Total LW heating
DTCOND=iris.load_cube(file_path_hist_Kdyn_heating, 'DTCOND')   #T tendency - moist processes (i.e. condensational heating)
DTV=iris.load_cube(file_path_hist_Kdyn_heating, 'DTV')	       #T vertical diffusion - molecular diffusion 
TTGW=iris.load_cube(file_path_hist_Kdyn_GWdrag, 'TTGW')	       #T tendency - gravity wave drag (heating due to parameterized GW)
TTGWSDF=iris.load_cube(file_path_hist_Kdyn_GWdrag_2, 'TTGWSDF')	 #component of TTGW due to vertical diffusion of Temp
TTGWSKE=iris.load_cube(file_path_hist_Kdyn_GWdrag_2, 'TTGWSKE')	 #component of TTGW due to conversion of kin energ into heat
GPH=iris.load_cube(file_path_hist_Kdyn, 'Z3')
T=iris.load_cube(file_path_hist_Kdyn, 'T')
P=iris.load_cube(file_path_hist_Kdyn, 'PMID')
Kdyn=iris.load_cube(file_path_hist_Kdyn, 'k_tot_tot')
Kzz_Kdyn=iris.load_cube(file_path_hist_Kdyn, 'EKGW')
Kwave_Kdyn=iris.load_cube(file_path_hist_Kdyn, 'k_wave_tot')


QRS_TOT_ctrl=iris.load_cube(file_path_hist_ctrl_heating, 'QRS_TOT')
QCP_ctrl=iris.load_cube(file_path_hist_ctrl_heating, 'QCP')
QRS_EUV_ctrl=iris.load_cube(file_path_hist_ctrl_heating, 'QRS_EUV')
QRS_CO2NIR_ctrl=iris.load_cube(file_path_hist_ctrl_heating, 'QRS_CO2NIR')
QRS_AUR_ctrl=iris.load_cube(file_path_hist_ctrl_heating, 'QRS_AUR')
QTHERMAL_ctrl=iris.load_cube(file_path_hist_ctrl_heating, 'QTHERMAL')
QRL_TOT_ctrl=iris.load_cube(file_path_hist_ctrl_heating, 'QRL_TOT')
DTCOND_ctrl=iris.load_cube(file_path_hist_ctrl_heating, 'DTCOND')
DTV_ctrl=iris.load_cube(file_path_hist_ctrl_heating, 'DTV')
TTGW_ctrl=iris.load_cube(file_path_hist_ctrl_GWdrag, 'TTGW')
TTGWSDF_ctrl=iris.load_cube(file_path_hist_ctrl_GWdrag_2, 'TTGWSDF')	 
TTGWSKE_ctrl=iris.load_cube(file_path_hist_ctrl_GWdrag_2, 'TTGWSKE')	
GPH_ctrl=iris.load_cube(file_path_hist_ctrl, 'Z3')
T_ctrl=iris.load_cube(file_path_hist_ctrl, 'T')
P_ctrl=iris.load_cube(file_path_hist_ctrl, 'PMID')
Kzz=iris.load_cube(file_path_hist_ctrl, 'EKGW')


#Plot QRS_TOT, total diabatic heating, GW heating and diffusion at specified latitude
im=6
#lat1=60
#lat2=90
lat1=-90
lat2=90
sec_to_day=86400


save_fig=True

##1. Plot QRS_TOT and budget terms
cube_list=iris.cube.CubeList([QRS_TOT*sec_to_day, QCP*sec_to_day, QTHERMAL*sec_to_day, QRS_EUV*sec_to_day, QRS_CO2NIR*sec_to_day, QRS_AUR*sec_to_day,   QRS_TOT_ctrl*sec_to_day, QCP_ctrl*sec_to_day, QTHERMAL_ctrl*sec_to_day, QRS_EUV_ctrl*sec_to_day, QRS_CO2NIR_ctrl*sec_to_day, QRS_AUR_ctrl*sec_to_day])
vertical_profile_mltp_lat_band(cube_list, GPH, month=im,  xmin=0, xmax=40, color=['black','black','black','black','black','black',  'gray','gray','gray','gray','gray','gray'], style=['-','dashdot','dotted','dashed',linestyles_dict['densely dashdotdotted'],linestyles_dict['loosely dashed'], '-','dashdot','dotted','dashed',linestyles_dict['densely dashdotdotted'], linestyles_dict['loosely dashed'] ], label=['Q_SW_TOT', 'QCP', 'QTHERMAL', 'Q_EUV', 'Q_CO2NIR', 'Q_AUR',  '', '', '', '', '', ''], xlabel='Heating rate (K/day)', title=' SW budget month: {im}, [{lat1},{lat2}]'.format(im=im, lat1=lat1, lat2=lat2), single=False, lat1=lat1, lat2=lat2)
if save_fig:
 plt.savefig('SW_heating_budget_month_{month}_lat_{lat1}{lat2}N.jpg'.format(month=im, lat1=abs(lat1), lat2=abs(lat2) ) )


##2. Plot Total Diabatic Heating Budget
Q_NET=QRS_TOT+QRL_TOT
Diab_TOT=Q_NET+DTCOND #total diabatic heating is combinatioon of total radiation and condensation terms
Q_NET_ctrl=QRS_TOT_ctrl+QRL_TOT_ctrl
Diab_TOT_ctrl=Q_NET_ctrl+DTCOND_ctrl
#in meosphere DTCOND=0, thus Diab_TOT=Q_NET, total diabatic heating  = total (net) radiative heating
cube_list=iris.cube.CubeList([ Diab_TOT*sec_to_day, QRS_TOT*sec_to_day, QRL_TOT*sec_to_day, Diab_TOT_ctrl*sec_to_day, QRS_TOT_ctrl*sec_to_day, QRL_TOT_ctrl*sec_to_day] )
vertical_profile_mltp_lat_band(cube_list, GPH, month=im,  xmin=-40, xmax=40, color=['black','black', 'black', 'gray', 'gray', 'gray'], style=['-','dashdot','dotted', '-','dashdot','dotted'], label=['Radiative_TOT', 'Q_SW', 'Q_LW', '', '',''], xlabel='Heating rate (K/day)', title='Radiative budget month: {im}, [{lat1},{lat2}]'.format(im=im, lat1=lat1, lat2=lat2), single=False, lat1=lat1, lat2=lat2)
if save_fig:
 plt.savefig('TotRadiative_heating_month_{month}_lat_{lat1}{lat2}N.jpg'.format(month=im, lat1=abs(lat1), lat2=abs(lat2) ) )

#3. Plot gravity wave heating rate (TTGW) and its components: TTGWSDF, TTGWSKE
cube_list=iris.cube.CubeList([TTGW*sec_to_day, TTGWSDF*sec_to_day, TTGWSKE*sec_to_day, TTGW_ctrl*sec_to_day, TTGWSDF_ctrl*sec_to_day, TTGWSKE_ctrl*sec_to_day] )
vertical_profile_mltp_lat_band(cube_list, GPH, month=im,  xmin=-40, xmax=40, color=['black','black', 'black', 'gray', 'gray', 'gray'], style=['-','dotted','dashdot', '-','dotted','dashdot'], label=['Total GW heating', 'Temp_diff', 'KinEn_conv', '', '', ''], xlabel='Heating rate (K/day)', title='GW heating budget month: {im}, [{lat1},{lat2}]'.format(im=im, lat1=lat1, lat2=lat2), single=False, lat1=lat1, lat2=lat2)
if save_fig:
 plt.savefig('TotalGWheating_month_{month}_lat_{lat1}{lat2}N.jpg'.format(month=im, lat1=abs(lat1), lat2=abs(lat2) ) )


#4. Plot K_zz and K_Dyn
cube_list=iris.cube.CubeList([Kzz, Kdyn] )
vertical_profile_mltp_lat_band(cube_list, GPH, month=im,  xmin=0, xmax=300, color=['gray','black'], style=['-','-'], label=['Kzz_CTRL', 'Kdyn'], xlabel='(m2/s)', title='Kzz and Kdyn: {im} at Lat: [{lat1},{lat2}]'.format(im=im, lat1=lat1, lat2=lat2), single=False, lat1=lat1, lat2=lat2)
if save_fig:
 plt.savefig('KzzKdyn_month_{month}_lat_{lat1}{lat2}N.jpg'.format(month=im, lat1=abs(lat1), lat2=abs(lat2) ) )




plt.show()

