import iris 
import iris.plot as iplt
import iris.quickplot as qplt
from iris.cube import Cube

import os
from mpl_toolkits.basemap import Basemap, maskoceans, shiftgrid

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from iris.experimental.equalise_cubes import *
import matplotlib.dates as mdates
import matplotlib.ticker as mtick
from matplotlib import ticker

import warnings
import glob
import os 

import pandas as pd
import scipy
from scipy import stats
from shiftedColorMap import*

############################################################

def plot_months(data_in, hgt, title, color_map, mid_point, levels, levels2, c_label, Log_scale, PMC):

 fig=plt.figure(figsize=(6,4))
 ax=plt.gca()

 #Create mesh_grid
 mnth=np.arange(0,len(data_in),1)
 mnth,hgt=np.meshgrid(mnth,hgt)
 data=(data_in).transpose([1,0])

 if Log_scale:
  shifted_cmap=shiftedColorMap(color_map, midpoint=mid_point, name='shifted')
  contours=plt.contourf(mnth, hgt, data, levels, cmap='shifted', locator=ticker.LogLocator()) #cmap=color_map
 else:
  shifted_cmap=shiftedColorMap(color_map, midpoint=mid_point, name='shifted')
  contours=plt.contourf(mnth, hgt, data/1000, levels, cmap='shifted' , extend='max')  #extend='both'                   
  contours_2=plt.contour(mnth, hgt, data/1000, levels2, colors='black', linewidths=2)
  plt.clabel(contours_2)

 if PMC: 
  PMC_season=data_in[:]   #take a copy of the input array
  PMC_season[2:10]=np.nan #fill in months non in PMC seasons with NaNs
  PMC_season=(PMC_season).transpose([1,0])
  contours_h=plt.contourf(mnth, hgt, PMC_season/1000, colors='grey', hatches=['///', '///']) #colors='none'


 cbar = plt.colorbar(contours, orientation='vertical')
 cbar.set_label(c_label)
 ax.set_ylim(70,110) 
 plt.xticks([0,1,2,3,4,5,6,7,8,9,10,11],['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A','S','O','N','D'])
 plt.xlabel('Months', fontsize=14)
 plt.ylabel('Height (Km)', fontsize=14)
 plt.title(title, fontsize=16)

 return()

##### Plot Lat vs months ##################################################################

def plot_lat_months_obs(data_in, lat, title, color_map, mid_point, levels, levels2, c_label , Log_scale):

 fig=plt.figure(figsize=(6,4))
 ax=plt.gca()

 #Create meshgrid
 mnth=np.arange(0,len(data_in),1)
 mnth,lat=np.meshgrid(mnth,lat)
 data_in=(data_in).transpose([1,0])

 if Log_scale:
  shifted_cmap=shiftedColorMap(color_map, midpoint=mid_point, name='shifted')
  contours=plt.contourf(mnth, lat, data_in, levels, cmap='shifted', locator=ticker.LogLocator()) #cmap=color_map
 else:
  shifted_cmap=shiftedColorMap(color_map, midpoint=mid_point, name='shifted')
  contours=plt.contourf(mnth,lat, data_in, levels, cmap='shifted' , extend='max')  #extend='both'                   
  contours_2=plt.contour(mnth, lat, data_in, levels2, colors='black', linewidths=2)
  plt.clabel(contours_2)

 cbar = plt.colorbar(contours, orientation='vertical')
 cbar.set_label(c_label)
 plt.xticks([0,1,2,3,4,5,6,7,8,9,10,11],['J', 'F', 'M', 'A', 'M', 'J','J','A','S','O','N','D']) 
 plt.xlabel('Months', fontsize=14)
 plt.ylabel('Latitude (deg)', fontsize=14)
 plt.title(title, fontsize=16)

 return()


##################################################### LOADING DATA #####################################################################################################
#Fe at Rothera 
Fe_Rothera=iris.load_cube('Fe_Rothera_2003_2004_2005.nc')

#Fe at South Pole 
Fe_SP=iris.load_cube('Fe_Southpole_MM.nc')

#Fe Urbana
Fe_Urbana=np.loadtxt('Fe_Urbana_19951996')
Fe_Urbana_hgt=Fe_Urbana[:,0]
#create multi-dim array with months  
Fe_Urbana_mnth=[]
for im in range (1,13,1): 
 month=Fe_Urbana[:,im]
 Fe_Urbana_mnth.append(month)
Fe_Urbana_mnth=np.array(Fe_Urbana_mnth)

#Na at South Pole
Na_SP=iris.load_cube('/Na_Southpole_MM.nc')

#Na at Urbana
Na_Urbana=iris.load_cube('Na_Urbana_19961998.nc', 'Na')

#Na Lat vs month data
Na_lat_month=np.loadtxt('Na_Ref_atm')
Na_lat=Na_lat_month[:,0]
Na_month=[]
for im in range (1,13,1): 
 month=Na_lat_month[:,im]
 Na_month.append(month)
Na_month=np.array(Na_month)


save_fig=False

#Rothera
levels_Fe_1=np.arange(0,15,0.5) #against obs
levels_Fe_2=np.arange(0,30,0.5) #against CTRL
plot_months(Fe_Rothera[:,:,0,0].data, Fe_Rothera.coord('height').points, title='Fe at Rothera', color_map=matplotlib.cm.rainbow, mid_point=0.5, levels=levels_Fe_2, levels2=7, c_label='10^3 cm-3', Log_scale=False, PMC=True)
if save_fig:
 plt.savefig('Fe_Rothera_levels2_obs.jpg')

#Urbana
levels_Fe_1=np.arange(0,15,0.5) #against obs
levels_Fe_2=np.arange(0,24,0.5) #against CTRL
levels_Na_1=np.arange(0,6,0.5)
levels_Na_2=np.arange(0,10,0.5)
plot_months(Fe_Urbana_mnth, Fe_Urbana_hgt, title='Fe at Urbana', color_map=matplotlib.cm.rainbow, mid_point=0.5, levels=levels_Fe_1, levels2=7, c_label='10^3 cm-3', Log_scale=False, PMC=False)
if save_fig:
 plt.savefig('Fe_Urbana_levels1_obs.jpg')
plot_months(Na_Urbana.data, Na_Urbana.coord('altitude').points, title='Na at Urbana', color_map=matplotlib.cm.rainbow, mid_point=0.5, levels=levels_Na_1, levels2=7, c_label='10^3 cm-3', Log_scale=False, PMC=False)
if save_fig:
 plt.savefig('Na_Urbana_levels1_obs.jpg')

#SouthPole
levels_Fe_1=np.arange(0,12,0.5) #against obs
levels_Fe_2=np.arange(0,36,0.5) #against CTRL
levels_Na_1=np.arange(0,8,0.5)
levels_Na_2=np.arange(0,14,0.5)
plot_months(Fe_SP[:,:,0,0].data, Fe_SP.coord('height').points, title='Fe at South Pole', color_map=matplotlib.cm.rainbow, mid_point=0.5, levels=levels_Fe_2, levels2=7, c_label='10^3 cm-3', Log_scale=False, PMC=True)
if save_fig:
 plt.savefig('Fe_SouthPole_levels2_obs.jpg')
plot_months(Na_SP[:,:,0,0].data, Na_SP.coord('height').points, title='Na at South Pole', color_map=matplotlib.cm.rainbow, mid_point=0.5, levels=levels_Na_2, levels2=7, c_label='10^3 cm-3', Log_scale=False, PMC=True)
if save_fig:
 plt.savefig('Na_SouthPole_levels2_obs.jpg')

#Na column abundance  Lat vs month plot 
levels_Na_1=np.arange(0,9,0.5) #against obs
levels_Na_2=np.arange(0,20,0.5) #against CTRL
plot_lat_months_obs(Na_month, Na_lat, title=' Na Reference Atmosphere Coulmn Abundance', color_map=matplotlib.cm.rainbow, mid_point=0.5, levels=levels_Na_1, levels2=7, c_label='10^9 cm-2' , Log_scale=False)
if save_fig:
 plt.savefig('Na_LatvsMonth_levels1_obs.jpg')

plt.show()
