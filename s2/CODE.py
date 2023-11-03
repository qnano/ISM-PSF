# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 15:22:38 2023

@author: localadmin
"""

import numpy as np
from scipy.io import loadmat

import matplotlib.pyplot as plt
import scipy.ndimage

fn="C:/research/spin disk/figure new/supply-napsf/scalarPSF.mat"
mat=loadmat(fn)
PSF_arr_scalar=mat['PSF_im_arr']

fn="C:/research/spin disk/figure new/supply-napsf/vecPSF.mat"
mat=loadmat(fn)
PSF_arr_vec=mat['PSF_im_arr']

import matplotlib.gridspec as gridspec
fig = plt.figure(figsize=(10,10),constrained_layout=True)
gs = gridspec.GridSpec(ncols=5, nrows=6, wspace=0.0, hspace=0.0, figure=fig)

if True:
    for i in range(5):
        f3_ax1 = fig.add_subplot(gs[0, i])
        f3_ax1.imshow(abs(PSF_arr_scalar[:,:,i])/np.amax(abs(PSF_arr_scalar[:,:,i])),cmap="afmhot")
        f3_ax1.axis("off")
        
        f3_ax1 = fig.add_subplot(gs[1, i])
        f3_ax1.imshow(abs(PSF_arr_vec[:,:,i])/np.amax(abs(PSF_arr_vec[:,:,i])),cmap="afmhot")
        f3_ax1.axis("off")
        

        f3_ax1 = fig.add_subplot(gs[2, i])
        
        f3_ax1.set_xlim([-700,700])
        f3_ax1.set_ylim([0,1.05])
        coor=np.linspace(-256*5,256*5,512)
        f3_ax1.set_yticks([0,0.5,1.0])
        xvec=abs(PSF_arr_vec[:,:,i])/np.amax(abs(PSF_arr_vec[:,:,i]))
        f3_ax1.plot(coor,xvec[256,:],color="k",label="vectorial")
        
        xscalar=abs(PSF_arr_scalar[:,:,i])/np.amax(abs(PSF_arr_scalar[:,:,i]))
        f3_ax1.plot(coor,xscalar[256,:],color="r",label="scalar")
        if i==0:
            plt.legend()


for i in range(5):
    f3_ax1 = fig.add_subplot(gs[3, i])
    
    f3_ax1.imshow(abs(PSF_arr_scalar[:,:,i+5])/np.amax(abs(PSF_arr_scalar[:,:,i+5])),cmap="afmhot")
    f3_ax1.axis("off")
    
    f3_ax1 = fig.add_subplot(gs[4, i])
    f3_ax1.imshow(abs(PSF_arr_vec[:,:,i+5])/np.amax(abs(PSF_arr_vec[:,:,i+5])),cmap="afmhot")
    f3_ax1.axis("off")
    
    
    f3_ax1 = fig.add_subplot(gs[5, i])
    f3_ax1.set_ylim([0,1.05])
    f3_ax1.set_yticks([0,0.5,1.0])
    f3_ax1.set_xlim([-300,300])
    coor=np.linspace(-256*5,256*5,512)
    xvec=abs(PSF_arr_vec[:,:,i+5])/np.amax(abs(PSF_arr_vec[:,:,i+5]))
    f3_ax1.plot(coor,xvec[256,:],color="k",label="vectorial")
    
    xscalar=abs(PSF_arr_scalar[:,:,i+5])/np.amax(abs(PSF_arr_scalar[:,:,i+5]))
    f3_ax1.plot(coor,xscalar[256,:],color="r",label="scalar")

fn="fig.svg"
plt.savefig(fn)