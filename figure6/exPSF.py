# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 12:20:22 2022

@author: localadmin
"""
import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt



import matplotlib.gridspec as gridspec
fig = plt.figure(figsize=(8,3),constrained_layout=True)
gs = gridspec.GridSpec(ncols=6, nrows=3, wspace=0.0, hspace=0.0, figure=fig)

fn="C:/research/spin disk/figure/figure6/multi_photon_exPSF.mat"
mat=loadmat(fn)
exPSF=mat['PSF_ex_arr']

for k in range(3):
    f3_ax1 = fig.add_subplot(gs[0, k])
    f3_ax1.axis("off")
    f3_ax1.imshow(exPSF[k],cmap="afmhot")

for k in range(3):
    f3_ax1 = fig.add_subplot(gs[1, k])
    f3_ax1.axis("off")
    f3_ax1.imshow(exPSF[k]**(k+1),cmap="afmhot")

fn="C:/research/spin disk/figure/figure6/multi_photon_ISM_PSF_img.mat"
mat=loadmat(fn)
ismPSF=mat['ism_PSF_arr']

for k in range(3):
    f3_ax1 = fig.add_subplot(gs[2, k])
    f3_ax1.axis("off")
    f3_ax1.imshow(ismPSF[k],cmap="afmhot")

x=np.linspace(-256*5,256*5,512)
f3_ax1 = fig.add_subplot(gs[:, 3:])

color1="r"
color2="navy"
color3="teal"
c=[color1,color2,color3]

for k in range(3):
    f3_ax1.plot(x,exPSF[k,256,:]**(k+1)/np.amax(exPSF[k,256,:]**(k+1)),color=c[k],label=str(k+1)+"-photon ex.")
plt.legend()
fn="C:/research/spin disk/figure/figure6/ex_ism_PSF.svg"
plt.savefig(fn)