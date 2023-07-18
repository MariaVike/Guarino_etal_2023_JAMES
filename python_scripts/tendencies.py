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

import warnings
import glob
import os 

import pandas as pd
import scipy
from scipy import stats


############## Compute N^2 at latitude #########################################################################
def compute_Nsq(TH, GPH, lat, lon, zmean):

 g=9.81
 if zmean:
  TH=TH.collapsed('longitude', iris.analysis.MEAN) 
  GPH=GPH.collapsed('longitude', iris.analysis.MEAN) 
  site=[('latitude', lat)]
  TH=TH.interpolate(site, iris.analysis.Linear())
  GPH=GPH.interpolate(site, iris.analysis.Linear())
 else:
  site=[('latitude', lat), ('longitude', lon)]
  TH=TH.interpolate(site, iris.analysis.Linear())
  GPH=GPH.interpolate(site, iris.analysis.Linear())


 Nsq_time=[]
 Ri_time=[]

 for it in range(0,GPH.shape[0]):
   Nsq=[]
   Ri=[]
   ilev=GPH.shape[1]-1 #python indexes: 0-69
   for iz in range(ilev-1, 0, -1):
     dth_dz=(TH[it,iz].data-TH[it,iz+1].data)/(GPH[it,iz].data-GPH[it,iz+1].data)
     Nsq.append( (g/TH.data[it,iz+1])*dth_dz )


   Nsq_time.append(Nsq[::-1]) #reverse array so that first element in array is model top


 hgt=GPH[:,1:69] #using centred differences

 return(Nsq_time, hgt) 
 
############### Plot vertical profile #######################################################################################################

def vertical_profile(D1, D2, change, hgt, color1, color2, color3, label1, label2, label3, title): 

  fig=plt.figure(tight_layout=True, figsize=(7, 6))
  fig.add_subplot(121)
  ax=plt.gca()
  plt.plot(D1, hgt/1000., c=color1, linewidth=4, linestyle='-', alpha=0.8, label=label1)
  plt.plot(D2, hgt/1000., c=color2, linewidth=4, linestyle='-', alpha=0.8, label=label2)
  plt.xlabel('1/s2', fontsize=14)
  plt.ylabel('GPH (Km)', fontsize=14)
  plt.title(title, fontsize=16)
  plt.ylim(50,100)
  plt.xlim(0,0.0016)
  fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
  fmt.set_powerlimits((0, 0))
  ax.xaxis.set_major_formatter(fmt)
  plt.legend()

  fig.add_subplot(122)
  plt.plot(change, hgt/1000., c=color3, linewidth=4, linestyle='-', alpha=0.8, label=label3)
  plt.xlabel('%', fontsize=14)
  plt.ylabel('GPH (Km)', fontsize=14)
  plt.ylim(50,100)
  plt.title('% change', fontsize=16)
  #plt.xlim()
  plt.legend()

  return()

############## Plot vertical profile at specified lat ################################################################################################################

def vertical_profile_mltp(D1, GPH, month, lat, color, style, label, xlabel, title, single): 

 if single:
  #compute zonal mean
  D1=D1.collapsed('longitude', iris.analysis.MEAN) 
  GPH=GPH.collapsed('longitude', iris.analysis.MEAN)
  #Extract value at specific lat
  site=[('latitude', lat)]
  D1=D1.interpolate(site, iris.analysis.Nearest())
  GPH=GPH.interpolate(site, iris.analysis.Nearest())
  D1_data=D1[month-1,:70].data
  GPH_data=GPH[month-1,:70].data/1000

  fig=plt.figure(tight_layout=True, figsize=(4, 6))
  ax=plt.gca()
  plt.plot(D1_data, GPH_data, c=color, linewidth=4, linestyle='-', alpha=0.8, label=label)
  
 else:
  #Compute zonal mean and extract value at specific lat 
  site=[('latitude', lat)]
  D1_data=iris.cube.CubeList()
  for ic in range(0, len(D1)):
   D1_zmean=D1[ic].collapsed('longitude', iris.analysis.MEAN) 
   D1_site=D1_zmean.interpolate(site, iris.analysis.Nearest())
   D1_data.append(D1_site[month-1,:70].data)
  GPH=GPH.collapsed('longitude', iris.analysis.MEAN)
  GPH=GPH.interpolate(site, iris.analysis.Nearest())
  GPH_data=GPH[month-1,:70].data/1000


  fig=plt.figure(tight_layout=True, figsize=(4, 6))
  ax=plt.gca()
  for i in range(0, len(D1_data)):
   plt.plot(D1_data[i], GPH_data, c=color[i], linewidth=4, linestyle=style[i], alpha=0.8, label=label[i])
   

 fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
 fmt.set_powerlimits((0, 0))

 plt.ylim(50,100)
 plt.xlabel(xlabel, fontsize=14)
 plt.ylabel('GPH (Km)', fontsize=14)
 ax.xaxis.set_major_formatter(fmt)
 plt.title(title, fontsize=16)
 plt.legend()
 plt.grid()

 return()

######################################################
#Na_Fe runs (CTRL & Kdyn) 
file_path_hist_ctrl='/path_to_data/CTRL_HIST_monthlyclim_20082011.nc'
file_path_hist_Kdyn='/path_to_data/Kdyn_HIST_monthlyclim_20082011.nc'
file_path_hist_ctrl_GWdrag='/path_to_data/CTRL_HIST_monthlyclim_20082011_GWdrag.nc'
file_path_hist_Kdyn_GWdrag='/path_to_data/Kdyn_HIST_monthlyclim_20082011_GWdrag.nc'


GPH=iris.load_cube(file_path_hist_Kdyn, 'Z3')
T=iris.load_cube(file_path_hist_Kdyn, 'T')
GPH_ctrl=iris.load_cube(file_path_hist_ctrl, 'Z3')
T_ctrl=iris.load_cube(file_path_hist_ctrl, 'T')

UTGWSPEC=iris.load_cube(file_path_hist_Kdyn_GWdrag, 'UTGWSPEC')
VTGWSPEC=iris.load_cube(file_path_hist_Kdyn_GWdrag, 'VTGWSPEC')
TTGWSPEC=iris.load_cube(file_path_hist_Kdyn_GWdrag, 'TTGWSPEC')   
UTEND1=iris.load_cube(file_path_hist_Kdyn_GWdrag, 'UTEND1')      #C&M U tendency   c < -40
UTEND2=iris.load_cube(file_path_hist_Kdyn_GWdrag, 'UTEND2')      #C&M U tendency  -40 < c < -15
UTEND3=iris.load_cube(file_path_hist_Kdyn_GWdrag, 'UTEND3')      #C&M U tendency  -15 < c <  15
UTEND4=iris.load_cube(file_path_hist_Kdyn_GWdrag, 'UTEND4')      #C&M U tendency   15 < c <  40
UTEND5=iris.load_cube(file_path_hist_Kdyn_GWdrag, 'UTEND5') 
BUTGWSPEC=iris.load_cube(file_path_hist_Kdyn_GWdrag, 'BUTGWSPEC')
BVTGWSPEC=iris.load_cube(file_path_hist_Kdyn_GWdrag, 'BVTGWSPEC')  
BTTGWSPEC=iris.load_cube(file_path_hist_Kdyn_GWdrag, 'BTTGWSPEC')   
BUTEND1=iris.load_cube(file_path_hist_Kdyn_GWdrag, 'BUTEND1')    
BUTEND2=iris.load_cube(file_path_hist_Kdyn_GWdrag, 'BUTEND2')   
BUTEND3=iris.load_cube(file_path_hist_Kdyn_GWdrag,'BUTEND3')     
BUTEND4=iris.load_cube(file_path_hist_Kdyn_GWdrag, 'BUTEND4')   
BUTEND5=iris.load_cube(file_path_hist_Kdyn_GWdrag, 'BUTEND5') 

UTGWSPEC_ctrl=iris.load_cube(file_path_hist_ctrl_GWdrag, 'UTGWSPEC')
VTGWSPEC_ctrl=iris.load_cube(file_path_hist_ctrl_GWdrag, 'VTGWSPEC')
TTGWSPEC_ctrl=iris.load_cube(file_path_hist_ctrl_GWdrag, 'TTGWSPEC')   
UTEND1_ctrl=iris.load_cube(file_path_hist_ctrl_GWdrag, 'UTEND1')      #C&M U tendency   c < -40
UTEND2_ctrl=iris.load_cube(file_path_hist_ctrl_GWdrag, 'UTEND2')      #C&M U tendency  -40 < c < -15
UTEND3_ctrl=iris.load_cube(file_path_hist_ctrl_GWdrag, 'UTEND3')      #C&M U tendency  -15 < c <  15
UTEND4_ctrl=iris.load_cube(file_path_hist_ctrl_GWdrag, 'UTEND4')      #C&M U tendency   15 < c <  40
UTEND5_ctrl=iris.load_cube(file_path_hist_ctrl_GWdrag, 'UTEND5') 

BUTGWSPEC_ctrl=iris.load_cube(file_path_hist_ctrl_GWdrag, 'BUTGWSPEC')
BVTGWSPEC_ctrl=iris.load_cube(file_path_hist_ctrl_GWdrag, 'BVTGWSPEC')  
BTTGWSPEC_ctrl=iris.load_cube(file_path_hist_ctrl_GWdrag, 'BTTGWSPEC')   
BUTEND1_ctrl=iris.load_cube(file_path_hist_ctrl_GWdrag, 'BUTEND1')    
BUTEND2_ctrl=iris.load_cube(file_path_hist_ctrl_GWdrag, 'BUTEND2')   
BUTEND3_ctrl=iris.load_cube(file_path_hist_ctrl_GWdrag,'BUTEND3')     
BUTEND4_ctrl=iris.load_cube(file_path_hist_ctrl_GWdrag, 'BUTEND4')   
BUTEND5_ctrl=iris.load_cube(file_path_hist_ctrl_GWdrag, 'BUTEND5') 

im=6
latit=40
longit=None

#compute and plot Nsq for Kdyn and CTRL
Nsq, hgt= compute_Nsq(TH, GPH,lat=latit ,lon=longit, zmean=True)
Nsq_ctrl, hgt_ctrl= compute_Nsq(TH_ctrl, GPH_ctrl, lat=latit, lon=longit, zmean=True)
#compute %change 
change=((np.array(Nsq[im-1])-np.array(Nsq_ctrl[im-1]))/np.array(Nsq_ctrl[im-1]))*100
#plot Nsq and %change
vertical_profile(Nsq[im-1], Nsq_ctrl[im-1], change, hgt[im-1].data, color1='red', color2='orange', color3='blue', label1='Nsq', label2='Nsq_ctrl', label3='%change', title='month: {im} , Lat: {latit}'.format(im=im,latit=latit))

#Plot tendencies
#Kdyn and CTRl - WAVES GENERATED BY FRONTS
cube_list=iris.cube.CubeList([UTGWSPEC, UTGWSPEC_ctrl, UTEND1, UTEND1_ctrl, UTEND5, UTEND5_ctrl])
vertical_profile_mltp(cube_list, GPH, month=im, lat=latit, color=['black','black','blue','blue','red','red'], style=['-','--', '-','--','-','--'], label=['UTEND_TOT', 'UTEND_TOT_ctrl', 'c < -40','c < -40 ctrl', 'c > 40', 'c > 40 ctrl'], xlabel='m/s2', title='Kdyn month: {im} at Lat: [{latit}]'.format(im=im, latit=latit), single=False)

plt.show()


