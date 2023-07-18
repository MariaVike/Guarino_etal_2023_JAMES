import iris
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
from scipy import optimize
from scipy import stats
from shiftedColorMap import*

############ Function to apply harmonic fit to input data: a1 annual amplitude, b1 annual phase constant, a2 semiannual amplitude, b2 semiannual phase constant, t time coordinate  ########################

def harmonic_function(t, a1, b1, a2, b2):

    A = a1*np.cos(2*np.pi*(t-b1)/365.)+ a2*np.cos(4*np.pi*(t-b2)/365.) #for year-around data
    #A = a1*np.cos(2*np.pi*(t-b1)/120.)+ a2*np.cos(4*np.pi*(t-b2)/120.)
    return(A) 

#############  Fit data_in using harmonic function defined above ################################################ 
def harmonic_fit(data_in): 

 #select data in specified range of latitudes 

 data_in_list=iris.cube.CubeList()
 lats=np.arange(-90,90,5)
 for il in lats:
  site=[('latitude', il)]
  data_in_list.append(data_in.interpolate(site, iris.analysis.Nearest()) )

 data_in=data_in_list.merge_cube()
 #print data_in.coords('latitude')

 harmonic_fit=[]
 harmonic_fit_ann=[]
 harmonic_fit_semiann=[]
 anomalies=[]

 annual_amp=[]
 annual_phase=[]
 semiann_amp=[]
 semiann_phase=[]

 #loop over latitudes
 for il in range(0, data_in.shape[0]):
  mean=np.mean(data_in[il,:].data)
  y=data_in[il,:].data-mean #compute anomalies

  #define x-axis, this is 12 months 
  x=np.arange(1,13,1)

  #Interpolate the data using a cubic spline to "new_length" samples
  new_length = 365 #(use 365 to go from months to days)  
  new_x = np.linspace(x.min(), x.max(), new_length)
  new_y = scipy.interpolate.interp1d(x, y, kind='cubic')(new_x)

  #x=new_x
  x=np.arange(0, len(new_y), 1)
  y=new_y


  popt,cov = scipy.optimize.curve_fit(harmonic_function, x,y)

  #take abs value of amplitude and compute phase
  if popt[0] >= 0:
   a1 = popt[0]
   b1 = np.mod(popt[1], 365) 
  else:
   a1 = abs(popt[0])
   b1 = np.mod(popt[1] + 182.5, 365)        

  if popt[2] >= 0:
   a2 = popt[2]
   b2 = np.mod(popt[3], 182.5) 
  else:
   a2 = abs(popt[2])
   b2 = np.mod(popt[3] + 92.5, 185)  


  #use coefficents just found to plot the corresponding harmonic fit (this is for plotting purposes)
  x_fit = np.arange(min(x), max(x)+1, 1)
  y_fit = harmonic_function(x_fit, a1, b1, a2, b2)
  #as above for annual oscillations only
  y_fit_ann=harmonic_function(x_fit,  a1, b1, 0., 0.)
  #as above for semi-annual oscillations only
  y_fit_semiann=harmonic_function(x_fit, 0., 0., a2, b2)

  harmonic_fit.append(y_fit)
  harmonic_fit_ann.append(y_fit_ann)
  harmonic_fit_semiann.append(y_fit_semiann)
  anomalies.append(y)

  annual_amp.append(a1)
  annual_phase.append(b1)
  semiann_amp.append(a2)
  semiann_phase.append(b2)

 return(lats, x_fit, anomalies, harmonic_fit, harmonic_fit_ann, harmonic_fit_semiann, annual_amp, annual_phase, semiann_amp, semiann_phase)

############### Plot amplitude of CTRL and Kdyn and compute Pearson corr coeff  #############################################################

def plot_amplitude(data_in, data_ctrl, data_MIF, lats, color, color2, label, style, ylabel, ylabel2, title, phase):

 if phase==False: 
  fig=plt.figure(tight_layout=True, figsize=(9, 4)) 
  ax=plt.gca()
  ax.plot(lats, data_in, c=color[0], linewidth=4, linestyle=style[0], alpha=0.8, label=label[0])
  ax.plot(lats, data_ctrl, c=color[1], linewidth=4, linestyle=style[1], alpha=0.8, label=label[1])
 
  plt.xlim(-80,80)
  ax.set_xlabel('Latitude (degrees)', fontsize=14)
  ax.set_ylabel(ylabel, fontsize=14, color=color[0])
  plt.title(title, fontsize=16)
  plt.legend()
  plt.grid()


 if phase==True: 
  fig=plt.figure(tight_layout=True, figsize=(9, 4)) 
  ax=plt.gca()
  ax.plot(lats, data_in, c=color[0], linewidth=4, linestyle=style[0], alpha=0.8, label=label[0])
  ax.plot(lats, data_ctrl, c=color[1], linewidth=4, linestyle=style[1], alpha=0.8, label=label[1])
 
  plt.xlim(-80,80)
  ax.set_xlabel('Latitude (degrees)', fontsize=14)
  ax.set_ylabel(ylabel, fontsize=14)
  plt.title(title, fontsize=16)
  plt.legend()
  plt.grid()


 return()

############# Load data ############Kdyn
#Na_Fe runs (CTRL & Kdyn) 
file_path_hist_ctrl='/path_to_data/CTRL_HIST_monthlyclim_20082011.nc'
file_path_hist_Kdyn='/path_to_data/Kdyn_HIST_monthlyclim_20082011.nc'

#KDYN
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


##Fit harmonic function to find annual and semi-annual oscillations of column abundance

#Load already computed column abundance for Na and Fe at all lats
Na_abn=iris.load_cube('Na_cAbn_Kdyn.nc') 
Na_ctrl_abn=iris.load_cube('Na_cAbn_ctrl.nc') 
Fe_abn=iris.load_cube('Fe_cAbn_Kdyn.nc') 
Fe_ctrl_abn=iris.load_cube('Fe_cAbn_ctrl.nc') 

##Na
'''
#Fit harmonic fucntion MIF 
lats, time, MIF_anomalies, MIF_harmonic_fit, MIF_harmonic_fit_ann, MIF_harmonic_fit_semiann, MIF_annual_amp, MIF_annual_phase, MIF_semiann_amp, MIF_semiann_phase=harmonic_fit(Na_MIF)
#Fit harmonic fucntion Column abnd Kdyn
lats, time, CAb_anomalies, CAb_harmonic_fit, CAb_harmonic_fit_ann, CAb_harmonic_fit_semiann, CAb_annual_amp, CAb_annual_phase, CAb_semiann_amp, CAb_semiann_phase=harmonic_fit(Na_abn/1E9)
#Fit harmonic fucntion Column abnd CTRL
lats, time, CAb_anomalies_ctrl, CAb_harmonic_fit_ctrl, CAb_harmonic_fit_ann_ctrl, CAb_harmonic_fit_semiann_ctrl, CAb_annual_amp_ctrl, CAb_annual_phase_ctrl, CAb_semiann_amp_ctrl, CAb_semiann_phase_ctrl=harmonic_fit(Na_ctrl_abn/1E9)
'''

##Fe
#'''
#Fit harmonic fucntion MIF 
lats, time, MIF_anomalies, MIF_harmonic_fit, MIF_harmonic_fit_ann, MIF_harmonic_fit_semiann, MIF_annual_amp, MIF_annual_phase, MIF_semiann_amp, MIF_semiann_phase=harmonic_fit(Fe_MIF)
#Fit harmonic function Column abnd Kdyn
lats, time, CAb_anomalies, CAb_harmonic_fit, CAb_harmonic_fit_ann, CAb_harmonic_fit_semiann, CAb_annual_amp, CAb_annual_phase, CAb_semiann_amp, CAb_semiann_phase=harmonic_fit(Fe_abn/1E9)
#Fit harmonic function Column abnd CTRL
lats, time, CAb_anomalies_ctrl, CAb_harmonic_fit_ctrl, CAb_harmonic_fit_ann_ctrl, CAb_harmonic_fit_semiann_ctrl, CAb_annual_amp_ctrl, CAb_annual_phase_ctrl, CAb_semiann_amp_ctrl, CAb_semiann_phase_ctrl=harmonic_fit(Fe_ctrl_abn/1E9)
#'''


#Plot annual amplitude and phase against lats for Kdyn and CTRL
metal='Fe'
plot_amplitude(CAb_annual_amp, CAb_annual_amp_ctrl, MIF_annual_phase, lats, color=['black','black'], color2='blue', label=['Kdyn', 'CTRL'], style=['-', '--'], ylabel='Column abnd (10^9 cm-2)', ylabel2='Injection rate (atoms cm-2 s-1)', title='Fe annual amplitude of anomalies', phase=False)
plot_amplitude(CAb_annual_phase, CAb_annual_phase_ctrl, MIF_annual_phase, lats, color=['black','black'], color2='blue', label=['Kdyn', 'CTRL'], style=['-', '--'], ylabel='angle(deg)', ylabel2='', title='Fe annual phase of anomalies', phase=True)


plt.show()
