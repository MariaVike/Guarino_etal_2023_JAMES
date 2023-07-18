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


########### Interpolate data on a Geopotential height vertical coordinate #################################################
#This function compute the zonal mean of the geopotential height (assuming GPH does not change much over the latitude) and for each grid-point in the dataset (data_in) interpolates the (GPH,latitude) profiles to a common, pre-defined GPH vertical coordinate. Returns a cube with dimensions (GPH,latitude) to be used to plot the zonal mean of data_in.

def to_gph(data_in,GPH):

 GPH=GPH.collapsed('longitude', iris.analysis.MEAN) 
 ##create the new GPH coord:
 cube_list=iris.cube.CubeList()
 lats=GPH.shape[1]
 for i in range(0, lats): 
  ##create a new cube for each latitude and add a new gph_coord to each
  gph_coord=iris.coords.DimCoord(GPH[:,i].data, long_name='gph')
  new_cube=iris.cube.Cube(data_in[:,i].data,  dim_coords_and_dims=[(gph_coord,0)])
  lat_coord= iris.coords.AuxCoord(data_in[:,i].coord('latitude').points, long_name='latitude')
  new_cube.add_aux_coord(lat_coord)
  ##interpolate each cube to a same reference vertical coordinate
  ref_gph = [('gph', np.linspace(100000, 40000, 100))]#meters (40 - 100 km)
  interp_cube = new_cube.interpolate(ref_gph, iris.analysis.Linear(extrapolation_mode='mask'))
  ##add all cubes to a list
  cube_list.append(interp_cube)

 ##merge list into one single cube with (latitude, gph) dimensions
 cube_gph=cube_list.merge_cube()

 return(cube_gph)
 
############################# Function to plot vertical- cross sections #########################################################################################################

def Plot_cross(data_in, GPH, title, color_map, mid_point, levels, levels2, c_label, mole,  Log_scale, model_levels, gph_coord, geom_alt):
    
 fig=plt.figure(figsize=(5,5))
 ax=plt.gca()

 #compute zonal mean of input data
 data_in=data_in.collapsed('longitude', iris.analysis.MEAN) 

 if model_levels:
  hgt=data_in.coord('atmosphere_hybrid_sigma_pressure_coordinate').points
  lat=data_in.coord('latitude').points
  lat,hgt=np.meshgrid(lat,hgt)
  data_in_data=data_in.data
  ylabel='Hybrid sigma (hPa)'
  ax.invert_yaxis()
  ax.set_ylim(0.0005,0.0000001)

 if gph_coord:
  data_in=to_gph(data_in, GPH)
  hgt=data_in.coord('gph').points
  lat=data_in.coord('latitude').points
  lat,hgt=np.meshgrid(lat,hgt/1000.)
  data_in_data=(data_in.data).transpose([1,0])
  ylabel='GPH (Km)'
  ax.set_ylim(50,100) #(75,120) #(50,100) (80,120)

 if geom_alt: #convert gph into geometric height (h) for comparison with obs. h=r*gph/r-gph
  data_in=to_gph(data_in, GPH)
  r=6371008.8  #earth's radius in m
  gph=data_in.coord('gph').points
  hgt=(r*gph)/(r-gph)
  lat=data_in.coord('latitude').points
  lat,hgt=np.meshgrid(lat,hgt/1000.)
  data_in_data=(data_in.data).transpose([1,0])
  ylabel='Altitude (Km)'
  ax.set_ylim(50,100) #(75,120) #(50,100) (80,120)
  
 if Log_scale:
  shifted_cmap=shiftedColorMap(color_map, midpoint=mid_point, name='shifted')
  contours=plt.contourf(lat, hgt, data_in_data, levels, cmap='shifted', locator=ticker.LogLocator())
 else:
  shifted_cmap=shiftedColorMap(color_map, midpoint=mid_point, name='shifted')
  contours=plt.contourf(lat, hgt, data_in_data, levels, cmap='shifted')                     
  contours_2=plt.contour(lat, hgt, data_in_data, levels2, colors='black', linewidths=1)
  plt.clabel(contours_2)

 if mole:
  fmt = matplotlib.ticker.ScalarFormatter(useMathText=True) #use this for mol/mol 
  fmt.set_powerlimits((0, 0))
 else:
  fmt=None

 cbar = plt.colorbar(contours, orientation='horizontal', format=fmt)
 cbar.set_label(c_label)
 plt.xlabel('Latitude', fontsize=14)
 plt.ylabel(ylabel, fontsize=14)
 plt.title(title, fontsize=16)


 return(data_in)

##### Plot MIPAS zonal means for T from 2 different products ####
def Plot_cross_MIPAS_T(data_in, data_in_2, title, color_map, mid_point, levels, levels2, c_label):
    
 fig=plt.figure(figsize=(5,5))
 ax=plt.gca()

 #compute zonal mean of input data
 data_in=data_in.collapsed('longitude', iris.analysis.MEAN) 
 data_in_2=data_in_2.collapsed('longitude', iris.analysis.MEAN) 
 hgt=data_in.coord('altitude').points
 lat=data_in.coord('latitude').points
 hgt_2=data_in_2.coord('altitude').points
 lat_2=data_in_2.coord('latitude').points
 lat,hgt=np.meshgrid(lat,hgt)
 lat_2,hgt_2=np.meshgrid(lat_2,hgt_2)
 data_in_data=(data_in.data).transpose([1,0])
 data_in_data_2=(data_in_2.data).transpose([1,0])
  

 shifted_cmap=shiftedColorMap(color_map, midpoint=mid_point, name='shifted')
 contours=plt.contourf(lat, hgt, data_in_data, levels, cmap='shifted')   
 contours_2=plt.contourf(lat_2, hgt_2, data_in_data_2, levels, cmap='shifted')                     
 contours_lev=plt.contour(lat, hgt, data_in_data, levels2, colors='black', linewidths=1)
 contours_lev_2=plt.contour(lat_2, hgt_2, data_in_data_2, levels2, colors='black', linewidths=1)
 plt.clabel(contours_lev)
 plt.clabel(contours_lev_2)


 cbar = plt.colorbar(contours, orientation='horizontal')
 #cbar = plt.colorbar(contours_2, orientation='vertical')
 cbar.set_label(c_label)
 ax.set_ylim(50,100) #(80,110) 
 plt.xlabel('Latitude', fontsize=14)
 plt.ylabel('Altitude (Km)', fontsize=14)
 plt.title(title, fontsize=16)


 return(data_in, data_in_2)

############################ Function to plot vertical- cross sections of anomalies #########################################################################################################

def Plot_cross_anom(data_in, title, color_map, mid_point, levels, levels2, c_label):
    
 fig=plt.figure(figsize=(5,5))
 ax=plt.gca()

 hgt=data_in.coord('gph').points
 lat=data_in.coord('latitude').points
 lat,hgt=np.meshgrid(lat,hgt/1000)
 data_in_data=(data_in.data).transpose([1,0])
 ylabel='GPH (Km)'
 ax.set_ylim(50,100) #(80,110)
  
 shifted_cmap=shiftedColorMap(color_map, midpoint=mid_point, name='shifted')
 contours=plt.contourf(lat, hgt, data_in_data, levels, cmap='shifted')                     
 #contours_2=plt.contour(lat, hgt, data_in_data, levels2, colors='black', linewidths=0.5)
 #plt.clabel(contours_2)

 cbar = plt.colorbar(contours, orientation='horizontal')
 cbar.set_label(c_label)
 plt.xlabel('Latitude', fontsize=14)
 plt.ylabel(ylabel, fontsize=14)
 plt.title(title, fontsize=16)


 return()

###### Find minimum in vertical column ####### 
def find_min(T, GPH, lat_in, MIPAS):

 if MIPAS:
  T=T.collapsed('longitude', iris.analysis.MEAN)
  where=[('latitude', lat_in)]
  T=T.interpolate(where, iris.analysis.Linear())
  T_min=T.collapsed('altitude', iris.analysis.MIN) 
  for il in range(0, T.shape[0]):
   if T[il].data==T_min.data:
    level_min=T[il].coord('altitude').points
  level_gh=level_min
 else:
  #compute zonal mean
  T=T.collapsed('longitude', iris.analysis.MEAN) 
  GPH=GPH.collapsed('longitude', iris.analysis.MEAN) 
  where=[('latitude', lat_in)]
  T=T.interpolate(where, iris.analysis.Linear())
  T_min=T.collapsed('atmosphere_hybrid_sigma_pressure_coordinate', iris.analysis.MIN)
  ##get model level where T is T_min
  for il in range(0, T.shape[0]):
   if T[il].data==T_min.data:
    level_min=T[il].coord('atmosphere_hybrid_sigma_pressure_coordinate').points

  where_gph=[('latitude', lat_in),('atmosphere_hybrid_sigma_pressure_coordinate', level_min)]
  GPH_min=GPH.interpolate(where_gph, iris.analysis.Linear())
  ##compute corrisponding geom_height
  r=6371008.8 
  level_gh=(r*GPH_min.data)/(r-GPH_min.data)

 return(T_min.data, level_gh)

##########################################################################################################################
#Load datasets and call functions ######################################################

#'''
#Na_Fe runs (CTRL & Kdyn) 
file_path_hist_ctrl='/path_to_data/CTRL_HIST_MIF_div2/CTRL_HIST_monthlyclim_20082011.nc'
file_path_hist_Kdyn='/path_to_data/Kdyn_HIST_monthlyclim_20082011.nc'

file_path_MIPAS_T= '/path_to_data/MIPAS_T_monthlyclim_20082011_V8R_662.nc'
file_path_MIPAS_T_2='/path_to_data/MIPAS_T_monthlyclim_20082011_661.nc'
file_path_MIPAS_CO2='/path_to_data/MIPAS_CO2_monthlyclim_20082011_V8R.nc'

CO2_ctrl=iris.load_cube(file_path_hist_ctrl, 'CO2')
T_ctrl=iris.load_cube(file_path_hist_ctrl, 'T')
GPH_ctrl=iris.load_cube(file_path_hist_ctrl, 'Z3')

CO2=iris.load_cube(file_path_hist_Kdyn, 'CO2')
T=iris.load_cube(file_path_hist_Kdyn, 'T')
GPH=iris.load_cube(file_path_hist_Kdyn, 'Z3')

CO2_MIPAS=iris.load_cube(file_path_MIPAS_CO2, 'CO2')
T_MIPAS=iris.load_cube(file_path_MIPAS_T, 'Temperature')
T_MIPAS_2=iris.load_cube(file_path_MIPAS_T_2, 'Temperature')

#Plot zonal means
im=6 #6,12

Kdyn=Plot_cross(T[im-1,:70], GPH[im-1,:70], title='month: {im} T zonal mean K_DYN'.format(im=im), color_map=matplotlib.cm.RdYlBu_r, mid_point=0.5, levels=np.arange(100,350,5), levels2=[120,140,160,180,200], c_label='K', mole=False, Log_scale=False, model_levels=False, gph_coord=False, geom_alt=True)

CTRL=Plot_cross(T_ctrl[im-1,:70], GPH_ctrl[im-1,:70], title='month: {im} T zonal mean CTRL'.format(im=im), color_map=matplotlib.cm.RdYlBu_r, mid_point=0.5, levels=np.arange(100,350,5), levels2=[120,140,160,180,200], c_label='K', mole=False, Log_scale=False, model_levels=False, gph_coord=False, geom_alt=True)

T_obs, T_obs_2=Plot_cross_MIPAS_T(T_MIPAS[im-1], T_MIPAS_2[im-1], title='month: {im} T zonal mean MIPAS'.format(im=im), color_map=matplotlib.cm.RdYlBu_r, mid_point=0.5, levels=np.arange(100,350,5), levels2=[120,140,160,180,200], c_label='K')


#find coldest point 
lat=85
T_min, GPH_min=find_min(T[im-1,:70], GPH[im-1,:70], lat, MIPAS=False)
print 'T min KDYN at', lat, '=', T_min, 'K at km', GPH_min
T_min_ctrl, GPH_min_ctrl=find_min(T_ctrl[im-1,:70], GPH_ctrl[im-1,:70], lat, MIPAS=False)
print 'T min CTRL at', lat, '=', T_min_ctrl, 'K at km', GPH_min_ctrl
T_min_MIPAS, hgt_min_MIPAS=find_min(T_MIPAS_2[im-1], None, lat, MIPAS=True)
print 'T min MIPAS at', lat, '=', T_min_MIPAS, 'K at km', hgt_min_MIPAS

plt.show()
