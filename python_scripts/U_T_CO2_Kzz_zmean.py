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
#This function compute the zonal mean of the geopotential height (assuming GPH does not change much for a fixed latitude) and for each grid-point in the dataset (data_in) interpolates the (GPH,latitude) profiles to a common, pre-defined GPH vertical coordinate. Returns a cube with dimensions (GPH,latitude) to be used to plot the zonal mean of data_in.

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
 
############################# Function to plot vertical cross-sections #########################################################################################################

def Plot_cross(data_in, GPH, title, color_map, mid_point, levels, levels2, c_label, Log_scale, model_levels, gph_coord):
    
 fig=plt.figure(figsize=(5,5))
 ax=plt.gca()

 #compute zonal mean of input data
 data_in=data_in.collapsed('longitude', iris.analysis.MEAN) 

 if model_levels:
  hgt=data_in.coord('atmosphere_hybrid_sigma_pressure_coordinate').points
  lat=data_in.coord('latitude').points
  lat,hgt=np.meshgrid(lat,hgt)
  data_in_data=data_in.data
  ylabel='(hPa)'
  ax.invert_yaxis()
  ax.set_ylim(10E-3, 10E-6)
  ax.set_yscale('log')

 if gph_coord:
  data_in=to_gph(data_in, GPH)
  hgt=data_in.coord('gph').points
  lat=data_in.coord('latitude').points
  lat,hgt=np.meshgrid(lat,hgt/1000.)
  data_in_data=(data_in.data).transpose([1,0])
  ylabel='GPH (Km)'
  ax.set_ylim(50,100)
  
 if Log_scale:
  shifted_cmap=shiftedColorMap(color_map, midpoint=mid_point, name='shifted')
  contours=plt.contourf(lat, hgt, data_in_data, levels, cmap='shifted', locator=ticker.LogLocator())
 else:
  shifted_cmap=shiftedColorMap(color_map, midpoint=mid_point, name='shifted')
  contours=plt.contourf(lat, hgt, data_in_data, levels, cmap='shifted')                     
  #contours_2=plt.contour(lat, hgt, data_in_data, levels2, colors='black', linewidths=1)
  #plt.clabel(contours_2)


 cbar = plt.colorbar(contours, orientation='horizontal')
 cbar.set_label(c_label)
 plt.xlabel('Latitude', fontsize=14)
 plt.ylabel(ylabel, fontsize=14)
 plt.title(title, fontsize=16)


 return(data_in)

################### Plot anomalies #############################################################
def Plot_cross_anom(data_in, title, color_map, mid_point, levels, levels2, c_label, contour2):
    
 fig=plt.figure(figsize=(5,5))
 ax=plt.gca()

 hgt=data_in.coord('gph').points
 lat=data_in.coord('latitude').points
 lat,hgt=np.meshgrid(lat,hgt/1000)
 data_in_data=(data_in.data).transpose([1,0])
 ylabel='GPH (Km)'
 ax.set_ylim(50,100) 
  
 shifted_cmap=shiftedColorMap(color_map, midpoint=mid_point, name='shifted')
 contours=plt.contourf(lat, hgt, data_in_data, levels, cmap='shifted')
 if contour2:                     
  contours_2=plt.contour(lat, hgt, data_in_data, levels2, colors='black', linewidths=0.5)
  plt.clabel(contours_2)

 cbar = plt.colorbar(contours, orientation='horizontal')
 cbar.set_label(c_label)
 plt.xlabel('Latitude', fontsize=14)
 plt.ylabel(ylabel, fontsize=14)
 plt.title(title, fontsize=16)


 return()

##### Plot Lat vs months ##################################################################

def plot_lat_months(data_in, title, color_map, mid_point, levels, levels2, c_label, Log_scale):

 fig=plt.figure(figsize=(6,4))
 ax=plt.gca()

 #mask zero values so they don't mess up zonal mean 
 data_in.data=np.ma.masked_where(data_in.data == 0, data_in.data)
 #compute zonal mean
 data_in=data_in.collapsed('longitude', iris.analysis.MEAN) 

 #Create meshgrid
 lat=data_in.coord('latitude').points
 mnth=data_in.coord('time').points
 mnth,lat=np.meshgrid(mnth,lat)
 data_in_data=(data_in.data).transpose([1,0])

 if Log_scale:
  shifted_cmap=shiftedColorMap(color_map, midpoint=mid_point, name='shifted')
  contours=plt.contourf(mnth, lat, data_in_data, levels, cmap='shifted', locator=ticker.LogLocator()) #cmap=color_map
  contours_2=plt.contour(mnth, lat, data_in_data, levels2, colors='black', locator=ticker.LogLocator())
  plt.clabel(contours_2)
 else:
  shifted_cmap=shiftedColorMap(color_map, midpoint=mid_point, name='shifted')
  contours=plt.contourf(mnth,lat, data_in_data, levels, cmap='shifted' , extend='max')  #extend='both'                   
  contours_2=plt.contour(mnth, lat, data_in_data, levels2, colors='black', linewidths=2)
  plt.clabel(contours_2)

 cbar = plt.colorbar(contours, orientation='vertical')
 cbar.set_label(c_label)
 plt.xticks(data_in.coord('time').points,['J', 'F', 'M', 'A', 'M', 'J','J','A','S','O','N','D']) 
 plt.xlabel('Months', fontsize=14)
 plt.ylabel('Latitude (deg)', fontsize=14)
 plt.title(title, fontsize=16)

 return(data_in)

############################################################### #Load datasets and call functions ######################################################
#Na_Fe runs (CTRL & Kdyn) 
file_path_hist_ctrl='/path_to_data/CTRL_HIST_monthlyclim_20082011.nc'
file_path_hist_Kdyn='/path_to_data/Kdyn_HIST_monthlyclim_20082011.nc'
file_path_hist_ctrl_2='/path_to_data/CTRL_HIST_monthlyclim_20082011_2.nc'
file_path_hist_Kdyn_2='/path_to_data/Kdyn_HIST_monthlyclim_20082011_2.nc'

K_zz=iris.load_cube(file_path_hist_Kdyn, 'EKGW')
K_zz_ctrl=iris.load_cube(file_path_hist_ctrl, 'EKGW')
K_dyn=iris.load_cube(file_path_hist_Kdyn, 'k_tot_tot')
K_wave=iris.load_cube(file_path_hist_Kdyn, 'k_wave_tot')

GPH=iris.load_cube(file_path_hist_Kdyn, 'Z3')
CO2=iris.load_cube(file_path_hist_Kdyn, 'CO2')
T=iris.load_cube(file_path_hist_Kdyn, 'T')
U=iris.load_cube(file_path_hist_Kdyn_2, 'U')
P=iris.load_cube(file_path_hist_Kdyn, 'PMID')

GPH_ctrl=iris.load_cube(file_path_hist_ctrl, 'Z3')
CO2_ctrl=iris.load_cube(file_path_hist_ctrl, 'CO2')
T_ctrl=iris.load_cube(file_path_hist_ctrl, 'T')
U_ctrl=iris.load_cube(file_path_hist_ctrl_2, 'U')
P_ctrl=iris.load_cube(file_path_hist_ctrl, 'PMID')


############################################# Months vs Lat plots ######################################################
#find index of vertical coord closest to specified pressure in hPa
pressure=0.005  #this is ~82km for Kdyn and ~83km for CTRL (results do not change using one or the other)
p_lev=GPH.coord('atmosphere_hybrid_sigma_pressure_coordinate').nearest_neighbour_index(pressure)
altitude=np.round(GPH[0,p_lev,:,0].data/1000)
print p_lev, altitude

save_fig=False

#set type of scale for colorbar and levels 
Log_scale=False
levels=np.arange(0,300,1)
levels2=5
#Log_scale=True
#levels=np.logspace(-1,3,13)
#levels2=[1,11,51, 151, 201]

K_wave_out=plot_lat_months(K_wave[:,p_lev], title='Kwave for Kdyn @ P: {pressure} hPa ~{altitude} km @'.format(pressure=pressure, altitude=altitude), color_map=matplotlib.cm.RdYlGn_r, mid_point=0.5, levels=levels, levels2=levels2, c_label='m2/s' , Log_scale=Log_scale)
if save_fig:
 plt.savefig('Kwave_Kdyn_LatvsMonth_at_{altitude}Km.jpg'.format(altitude=altitude))

K_zz_out=plot_lat_months(K_zz[:,p_lev], title='Kzz for Kdyn @ P: {pressure} hPa ~{altitude} km @'.format(pressure=pressure, altitude=altitude), color_map=matplotlib.cm.RdYlGn_r, mid_point=0.5, levels=levels, levels2=levels2, c_label='m2/s' , Log_scale=Log_scale)
if save_fig:
 plt.savefig('Kzz_Kdyn_LatvsMonth_at_{altitude}Km.jpg'.format(altitude=altitude))

K_dyn_out=plot_lat_months(K_dyn[:,p_lev], title='Kdyn for Kdyn @ P: {pressure} hPa ~{altitude} km @'.format(pressure=pressure, altitude=altitude), color_map=matplotlib.cm.RdYlGn_r, mid_point=0.5, levels=levels, levels2=levels2, c_label='m2/s' , Log_scale=Log_scale)
if save_fig:
 plt.savefig('Kdyn_Kdyn_LatvsMonth_at_{altitude}Km.jpg'.format(altitude=altitude))

K_zz_ctrl_out=plot_lat_months(K_zz_ctrl[:,p_lev], title='Kzz for CTRL @ P: {pressure} hPa ~{altitude} km @'.format(pressure=pressure, altitude=altitude), color_map=matplotlib.cm.RdYlGn_r, mid_point=0.5, levels=levels, levels2=levels2, c_label='m2/s' , Log_scale=Log_scale)
if save_fig:
 plt.savefig('Kzz_CTRL_LatvsMonth_at_{altitude}Km.jpg'.format(altitude=altitude))

#compute and plot kdyn/kzz ratio
diffsv_ratio=(K_dyn_out/K_zz_ctrl_out)
plot_lat_months_anom(diffsv_ratio, title='Kdyn/Kzz_ctrl ratio @ P: {pressure} hPa ~{altitude} km @'.format(pressure=pressure, altitude=altitude), color_map=matplotlib.cm.cubehelix, mid_point=0.5, levels=np.arange(0,7,0.2), levels2=[1,2,4,6], c_label='times' , Log_scale=False)
if save_fig:
 plt.savefig('KdynKzz_ctrl_ratio_LatvsMonth_at_{altitude}Km.jpg'.format(altitude=altitude))
#compute and plot kzz/kzz_ctrl ratio
diffsv_ratio_Kzz=(K_zz_out/K_zz_ctrl_out)
plot_lat_months_anom(diffsv_ratio_Kzz, title='Kzz/Kzz_ctrl ratio @ P: {pressure} hPa ~{altitude} km @'.format(pressure=pressure, altitude=altitude), color_map=matplotlib.cm.cubehelix, mid_point=0.5, levels=np.arange(0,3,0.1), levels2=[0.5,1,2], c_label='times' , Log_scale=False)
if save_fig:
 plt.savefig('Kzz_ratio_LatvsMonth_at_{altitude}Km.jpg'.format(altitude=altitude))

############################################# Zonal means  ############################################################

im=6
save_fig=False

#'''
#CO2
Kdyn=Plot_cross(CO2[im-1,:70]*1E6, GPH[im-1,:70], title='month: {im} CO2 zonal mean Kdyn'.format(im=im), color_map=matplotlib.cm.gist_rainbow_r, mid_point=0.5, levels=np.arange(90,400,2), levels2=[120,160,200,240,280,320,360], c_label='ppm', fmt=False, Log_scale=False, model_levels=False, gph_coord=True)
if save_fig:
 plt.savefig('CO2_zmean_Kdyn_month_{month}.jpg'.format(month=im))
CTRL=Plot_cross(CO2_ctrl[im-1,:70]*1E6, GPH_ctrl[im-1,:70], title='month: {im} CO2 zonal mean CTRL'.format(im=im), color_map=matplotlib.cm.gist_rainbow_r, mid_point=0.5, levels=np.arange(90,400,2), levels2=[120,160,200,240,280,320,360], c_label='ppm', fmt=False, Log_scale=False, model_levels=False, gph_coord=True)
if save_fig:
 plt.savefig('CO2_zmean_CTRL_month_{month}.jpg'.format(month=im))
#plot anomalies
#levels=np.arange(-8,38,1) #June
#mid_point=0.2
levels=50 #Dec
mid_point=0.1
Plot_cross_anom(Kdyn-CTRL, title='month: {im} Kdyn-CTRL CO2 anomalies'.format(im=im), color_map=matplotlib.cm.RdBu_r, mid_point=mid_point, levels=levels, levels2=[0], c_label='ppm', fmt=False, contour2=False)
if save_fig:
 plt.savefig('CO2_zmean_KdynCTRL_anom_month_{month}.jpg'.format(month=im))
#'''

#T
'''
Kdyn=Plot_cross(T[im-1,:70], GPH[im-1,:70], title='month: {im} T zonal mean Kdyn'.format(im=im), color_map=matplotlib.cm.RdYlBu_r, mid_point=0.5, levels=np.arange(100,350,5), levels2=[120,140,160,180,200], c_label='K', fmt=False, Log_scale=False, model_levels=False, gph_coord=True)
if save_fig:
 plt.savefig('T_zmean_Kdyn_month_{month}.jpg'.format(month=im))
CTRL=Plot_cross(T_ctrl[im-1,:70], GPH_ctrl[im-1,:70], title='month: {im} T zonal mean CTRL'.format(im=im), color_map=matplotlib.cm.RdYlBu_r, mid_point=0.5, levels=np.arange(100,350,5), levels2=[120,140,160,180,200], c_label='K', fmt=False, Log_scale=False, model_levels=False, gph_coord=True)
if save_fig:
 plt.savefig('T_zmean_CTRL_month_{month}.jpg'.format(month=im))
plot anamoalies
levels=np.arange(-45, 48, 1) #June
#mid_point=0.5
#levels=np.arange(-40, 30, 1) #Dec
#mid_point=0.6
Plot_cross_anom(Kdyn-CTRL, title='month: {im} Kdyn-CTRL T anomalies'.format(im=im), color_map=matplotlib.cm.RdBu_r, mid_point=mid_point, levels=levels, levels2=[0], c_label='K', fmt=False, contour2=False)
if save_fig:
 plt.savefig('T_zmean_KdynCTRL_anom_month_{month}.jpg'.format(month=im))
'''

#U
'''
Kdyn=Plot_cross(U[im-1,:70], GPH[im-1,:70], title='month: {im} U zonal mean Kdyn'.format(im=im), color_map=matplotlib.cm.RdBu_r, mid_point=0.4, levels=np.arange(-82,122,2), levels2=5, c_label='m/s', fmt=False, Log_scale=False, model_levels=False, gph_coord=True)
if save_fig:
 plt.savefig('U_zmean_Kdyn_month_{month}.jpg'.format(month=im))
CTRL=Plot_cross(U_ctrl[im-1,:70], GPH_ctrl[im-1,:70], title='month: {im} U zonal mean CTRL'.format(im=im), color_map=matplotlib.cm.RdBu_r, mid_point=0.4, levels=np.arange(-82,122,2), levels2=5, c_label='m/s', fmt=False, Log_scale=False, model_levels=False, gph_coord=True)
if save_fig:
 plt.savefig('U_zmean_CTRL_month_{month}.jpg'.format(month=im))
#plot anamoalies
Plot_cross_anom(Kdyn-CTRL, title='month: {im} Kdyn-CTRL U anomalies'.format(im=im), color_map=matplotlib.cm.RdBu_r, mid_point=0.5, levels=np.arange(-30,30,2), levels2=5, c_label='m/s', fmt=False, contour2=True)
if save_fig:
 plt.savefig('U_zmean_KdynCTRL_anom_month_{month}.jpg'.format(month=im))
'''

plt.show()



