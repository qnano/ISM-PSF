# -*- coding: utf-8 -*-
"""
Created on Sun Sep 25 18:54:32 2022

@author: localadmin
"""
import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt



import matplotlib.gridspec as gridspec
fig = plt.figure(figsize=(9,3),constrained_layout=True)
gs = gridspec.GridSpec(ncols=2, nrows=1, wspace=0.0, hspace=0.0, figure=fig)


f3_ax1 = fig.add_subplot(gs[0, 0])
f3_ax1.set_xlabel("NA")
f3_ax1.set_ylabel("FWHM of PSF (nm)")

NA_arr=np.linspace(0.1,1.3,31)
color1="r"
color2="navy"
color3="teal"

fn="C:/research/spin disk/figure/figure3/ISM_PSFs_width_over_NA_vect.mat"
mat=loadmat(fn)
ISM_v=mat['PSF_im_fitting_arr'][0,1:]*2.355

fn="C:/research/spin disk/figure/figure3/ISM_PSFs_width_over_NA_scalar.mat"
mat=loadmat(fn)
ISM_s=mat['PSF_im_fitting_arr'][0,1:]*2.355

f3_ax1.plot(NA_arr[1:],ISM_v,label="vectorial",linewidth=3)
f3_ax1.plot(NA_arr[1:],ISM_s,label="scalar",linewidth=3)
f3_ax1.legend()


f3_ax1 = fig.add_subplot(gs[0, 1])
f3_ax1.set_xlabel("NA")



#f3_ax1.set_ylabel("FWHM of PSF (nm)")

#NA_arr=np.linspace(0.1,1.3,15)
color1="r"
color2="navy"
color3="teal"

fn="C:/research/spin disk/figure/figure3/ISM_PSFs_width_over_NA_vect.mat"
mat=loadmat(fn)
ISM_v=mat['PSF_im_fitting_arr'][0,1:]*2.355

fn="C:/research/spin disk/figure/figure3/ISM_PSFs_width_over_NA_scalar.mat"
mat=loadmat(fn)
ISM_s=mat['PSF_im_fitting_arr'][0,1:]*2.355

f3_ax1.plot(NA_arr[1:],(ISM_v-ISM_s)/ISM_s,linewidth=3,color="r")
f3_ax1.tick_params(axis ='y', labelcolor = 'red') 

ax2 = f3_ax1.twinx() 
ax2.set_ylabel('diff', color = 'blue') 
ax2.plot(NA_arr[1:], ISM_v-ISM_s,linewidth=3, color = 'blue') 
ax2.tick_params(axis ='y', labelcolor = 'blue') 

fn="C:/research/spin disk/figure/figure3/all.svg"
plt.savefig(fn)
