import iris 
import iris.plot as iplt
import iris.quickplot as qplt
from iris.cube import Cube

import os
os.environ['PROJ_LIB'] = './.conda/pkgs/pyproj-1.9.4-py27_0/lib/python2.7/site-packages/pyproj/data/'
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

########## Interpolate data on a Geopotential height vertical coordinate #################################################
def to_gph_months(data_in,GPH):

 ##create the new GPH coord:
 cube_list=iris.cube.CubeList()
 months=GPH.shape[0]
 for i in range(0, months): 
  ##create a new cube for each latitude and add a new gph_coord to each
  gph_coord=iris.coords.DimCoord(GPH[i,:].data, long_name='gph')
  new_cube=iris.cube.Cube(data_in[i,:].data,  dim_coords_and_dims=[(gph_coord,0)])
  time_coord= iris.coords.AuxCoord(data_in[i,:].coord('time').points, long_name='time')
  new_cube.add_aux_coord(time_coord)
  ##interpolate each cube to a same reference vertical coordinate
  ref_gph = [('gph', np.linspace(140000, 0, 100))]#meters (0 - 140 km)
  interp_cube = new_cube.interpolate(ref_gph, iris.analysis.Linear(extrapolation_mode='mask'))
  ##add all cubes to a list
  cube_list.append(interp_cube)

 ##merge list into one single cube with (latitude, gph) dimensions
 cube_gph=cube_list.merge_cube()

 return(cube_gph)

#############################################################

def plot_months(data_in, GPH, lat, lon, title, color_map, mid_point, levels, levels2, c_label, ymin, ymax, extend, color_c, Log_scale, PMC):

 fig=plt.figure(figsize=(6,4))
 ax=plt.gca()

 #Extract value at specific location 
 site=[('latitude', lat), ('longitude', lon)]
 data_in=data_in.interpolate(site, iris.analysis.Nearest())
 GPH=GPH.interpolate(site, iris.analysis.Nearest())

 #Interpolate vertical profile for each month onto a gph vertical coordinate
 data_in=to_gph_months(data_in,GPH)
 #convert gph into geometric height (h) for comparison with obs. h=r*gph/r-gph
 r=6371008.8  #earth's radius in m
 gph=data_in.coord('gph').points
 hgt=(r*gph)/(r-gph)

 mnth=data_in.coord('time').points
 mnth,hgt=np.meshgrid(mnth,hgt/1000.)
 data_in_data=(data_in.data).transpose([1,0])

 if Log_scale:
  shifted_cmap=shiftedColorMap(color_map, midpoint=mid_point, name='shifted')
  contours=plt.contourf(mnth, hgt, data_in_data, levels, cmap='shifted', locator=ticker.LogLocator()) #cmap=color_map
 else:
  shifted_cmap=shiftedColorMap(color_map, midpoint=mid_point, name='shifted')
  contours=plt.contourf(mnth, hgt, data_in_data, levels, cmap='shifted' , extend=extend)  #extend='both' 
  #contours.cmap.set_over('yellow')  #use only for Var_T plots                 
  contours_2=plt.contour(mnth, hgt, data_in_data, levels2, colors=color_c, linewidths=2)
  plt.clabel(contours_2)

 if PMC: 
  PMC_season=data_in.data[:].copy()  #take a copy of the input array
  PMC_season[2:10]=np.nan            #fill in months non in PMC seasons with NaNs
  PMC_season=(PMC_season).transpose([1,0])
  contours_h=plt.contourf(mnth, hgt, PMC_season/1000, colors='grey', hatches=['///', '///']) #colors='none'

 cbar = plt.colorbar(contours, orientation='vertical')
 cbar.set_label(c_label)
 ax.set_ylim(ymin, ymax)
 #plt.xticks([31,59,90,120,151,181,212,243,273,304,334,365],['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A','S','O','N','D']) #'J'
 plt.xticks(data_in.coord('time').points,['J', 'F', 'M', 'A', 'M', 'J','J','A','S','O','N','D']) 
 plt.xlabel('Months', fontsize=14)
 plt.ylabel('Altitude (Km)', fontsize=14)
 plt.title(title, fontsize=16)

 return(data_in)


############# Load data ###########
#Na_Fe runs (CTRL & Kdyn) 
file_path_hist_ctrl='/path_to_data/CTRL_HIST_monthlyclim_20082011.nc'
file_path_hist_Kdyn='/path_to_data/Kdyn_HIST_monthlyclim_20082011.nc'

#Kdyn
T=iris.load_cube(file_path_hist_Kdyn, 'T')
GPH=iris.load_cube(file_path_hist_Kdyn, 'Z3')
Na=iris.load_cube(file_path_hist_Kdyn, 'Na')
Fe=iris.load_cube(file_path_hist_Kdyn, 'Fe')
P=iris.load_cube(file_path_hist_Kdyn, 'PMID')
#CTRL
T_ctrl=iris.load_cube(file_path_hist_ctrl, 'T')
GPH_ctrl=iris.load_cube(file_path_hist_ctrl, 'Z3')
Na_ctrl=iris.load_cube(file_path_hist_ctrl, 'Na')
Fe_ctrl=iris.load_cube(file_path_hist_ctrl, 'Fe')
P_ctrl=iris.load_cube(file_path_hist_ctrl, 'PMID')



############# Plot Na and Fe ###########
##at Rothera
#'''
latit=-68
longit=292
levels_Fe_1=np.arange(0,15,0.5) #against obs
levels_Fe_2=np.arange(0,30,0.5) #against CTRL
levels_Fe=levels_Fe_2
levels_Fe_anom=np.arange(-22,4,0.5)
mid_point_Fe=0.85 
levels_Na_1=np.arange(0,8,0.5)
levels_Na_2=np.arange(0,12,0.5)
levels_Na=levels_Na_2
levels_Na_anom=np.arange(-10,2,0.25)
mid_point_Na=0.85
PMC=True
site='Rothera'
#'''

##at Urbana
'''
latit=40
longit=-88
levels_Fe_1=np.arange(0,15,0.5) #against obs
levels_Fe_2=np.arange(0,24,0.5) #against CTRL
levels_Fe=levels_Fe_2
levels_Fe_anom=np.arange(-16,1,0.1)
mid_point_Fe=0.95 
levels_Na_1=np.arange(0,6,0.5)
levels_Na_2=np.arange(0,10,0.5)
levels_Na=levels_Na_2
levels_Na_anom=np.arange(-7,1,0.1)
mid_point_Na=0.95
PMC=False
site='Urbana'
'''

## at South Pole
'''
latit=-90
longit=59
levels_Fe_1=np.arange(0,12,0.5) #against obs
levels_Fe_2=np.arange(0,36,0.5) #against CTRL
levels_Fe=levels_Fe_2
levels_Fe_anom=np.arange(-28,8,0.5)
mid_point_Fe=0.75 
levels_Na_1=np.arange(0,8,0.5)
levels_Na_2=np.arange(0,14,0.5)
levels_Na=levels_Na_2
levels_Na_anom=np.arange(-12,4,0.25)
mid_point_Na=0.75
PMC=True
site='SouthPole'
'''

ymin=70 
ymax=110

save_fig=False

## !!! Fe !!! ###
#'''
##plot concentration in cm-3 computing n=P*Av/R*T
R=8.31446261815324
Av=6.023*1E23 
to_cm3=1E-6
n=(P*Av/(R*T))*to_cm3
n_ctrl=(P_ctrl*Av/(R*T_ctrl))*to_cm3

CTRL = plot_months((Fe_ctrl[:,:70,:,:]*n_ctrl)/1000., GPH_ctrl[:,:70,:,:], lat=latit, lon=longit, title=' Fe CTRL, Lat: {latit} , Lon: {longit}'.format(latit=latit, longit=longit), color_map=matplotlib.cm.rainbow, mid_point=0.5, levels=levels_Fe, levels2=7, c_label='10^3 cm-3' , ymin=ymin, ymax=ymax, extend='max', color_c='black', Log_scale=False, PMC=PMC)
if save_fig:
 plt.savefig('Fe_{site}_levels2_model_CTRL_geom_alt.jpg'.format(site=site))

Kdyn = plot_months((Fe[:,:70,:,:]*n)/1000., GPH[:,:70,:,:], lat=latit, lon=longit, title=' Fe K_DYN, Lat: {latit} , Lon: {longit}'.format(latit=latit, longit=longit), color_map=matplotlib.cm.rainbow, mid_point=0.5, levels=levels_Fe, levels2=7, c_label='10^3 cm-3' , ymin=ymin, ymax=ymax, extend='max', color_c='black', Log_scale=False, PMC=PMC)
if save_fig:
 plt.savefig('Fe_{site}_levels2_model_Kdyn_geom_alt.jpg'.format(site=site))


## !!! Na !!! ###
##plot concentration in cm-3 computing n=P*Av/R*T
R=8.31446261815324
Av=6.023*1E23 
to_cm3=1E-6
n=(P*Av/(R*T))*to_cm3
n_ctrl=(P_ctrl*Av/(R*T_ctrl))*to_cm3

CTRL = plot_months((Na_ctrl[:,:70,:,:]*n_ctrl)/1000., GPH_ctrl[:,:70,:,:], lat=latit, lon=longit, title=' Na CTRL, Lat: {latit} , Lon: {longit} '.format(latit=latit, longit=longit), color_map=matplotlib.cm.rainbow, mid_point=0.5, levels=levels_Na, levels2=7, c_label='10^3 cm-3' , ymin=ymin, ymax=ymax, extend='max', color_c='black', Log_scale=False, PMC=PMC)
if save_fig:
 plt.savefig('Na_{site}_levels2_model_CTRL_geom_alt.jpg'.format(site=site))

Kdyn = plot_months((Na[:,:70,:,:]*n)/1000., GPH[:,:70,:,:], lat=latit, lon=longit, title=' Na K_DYN, Lat: {latit} , Lon: {longit} '.format(latit=latit, longit=longit), color_map=matplotlib.cm.rainbow, mid_point=0.5, levels=levels_Na, levels2=7, c_label='10^3 cm-3' , ymin=ymin, ymax=ymax, extend='max', color_c='black', Log_scale=False, PMC=PMC)
if save_fig:
 plt.savefig('Na_{site}_levels2_model_Kdyn_geom_alt.jpg'.format(site=site))


plt.show()
