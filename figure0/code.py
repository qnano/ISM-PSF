# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 15:46:36 2022

@author: localadmin
"""
import numpy as np
from scipy.io import loadmat

import matplotlib.pyplot as plt
import scipy.ndimage

fn="C:/research/spin disk/figure new/figure0/xlinear_expsf.mat"
mat=loadmat(fn)
PSFx=mat['PSF_ex']
PSFx=abs(PSFx)/np.amax(abs(PSFx))
fn="C:/research/spin disk/figure new/figure0/ylinear_expsf.mat"
mat=loadmat(fn)
PSFy=mat['PSF_ex']
PSFy=abs(PSFy)/np.amax(abs(PSFy))
fn="C:/research/spin disk/figure new/figure0/circular_expsf.mat"
mat=loadmat(fn)
PSFc=mat['PSF_ex']
PSFc=abs(PSFc)/np.amax(abs(PSFc))
fn="C:/research/spin disk/figure new/figure0/scalar_expsf.mat"
mat=loadmat(fn)
PSFs=mat['PSF_ex']
PSFs=abs(PSFs)/np.amax(abs(PSFs))
import matplotlib.gridspec as gridspec
fig = plt.figure(figsize=(10,4),constrained_layout=True)
gs = gridspec.GridSpec(ncols=4, nrows=2, wspace=0.0, hspace=0.0, figure=fig)

f3_ax1 = fig.add_subplot(gs[0, 0])
f3_ax1.imshow(PSFs,cmap="afmhot")
f3_ax1.axis("off")

f3_ax1 = fig.add_subplot(gs[0, 1])
f3_ax1.imshow(PSFx,cmap="afmhot")
f3_ax1.axis("off")

f3_ax1 = fig.add_subplot(gs[0, 2])
f3_ax1.imshow(PSFy,cmap="afmhot")
f3_ax1.axis("off")

f3_ax1 = fig.add_subplot(gs[0, 3])
f3_ax1.imshow(PSFc,cmap="afmhot")
f3_ax1.axis("off")


fn="C:/research/spin disk/figure new/figure0/scalar_empsf.mat"
mat=loadmat(fn)
PSFsem=mat['PSF_em']
PSFsem=abs(PSFsem)/np.amax(abs(PSFsem))

fn="C:/research/spin disk/figure new/figure0/vectorial_empsf.mat"
mat=loadmat(fn)
PSFvem=mat['PSF_em']
PSFvem=abs(PSFvem)/np.amax(abs(PSFvem))
f3_ax1 = fig.add_subplot(gs[1, 0])
f3_ax1.imshow(PSFsem,cmap="afmhot")
f3_ax1.axis("off")

f3_ax1 = fig.add_subplot(gs[1, 1])
f3_ax1.imshow(PSFvem,cmap="afmhot")
f3_ax1.axis("off")

unit=640/(1.3)
pos=np.linspace(-5*100,5*100,200)/(unit)
f3_ax1 = fig.add_subplot(gs[1, 2])
f3_ax1.plot(pos,PSFs[256-100:256+100,256],color="r",label="scalar")
f3_ax1.plot(pos,PSFx[256-100:256+100,256],color="g",label="x-linear")
f3_ax1.plot(pos,PSFy[256-100:256+100,256],color="b",label="y-linear")
f3_ax1.plot(pos,PSFc[256-100:256+100,256],color="k",label="circular")
f3_ax1.set_title("excitation PSF")
f3_ax1.set_ylim([0,1.1])
plt.legend()

c=-6
unit=640/(1.3)
pos=np.linspace(-5*100,5*100,200)/(unit)
f3_ax1 = fig.add_subplot(gs[1, 3])
f3_ax1.plot(pos,PSFsem[256-100-c:256+100-c,256]/np.amax(PSFsem[256-100-c:256+100-c,256]),color="r",label="scalar")
f3_ax1.plot(pos,PSFvem[256-100-c:256+100-c,256]/np.amax(PSFvem[256-100-c:256+100-c,256]),color="k",label="vectorial")
f3_ax1.set_ylim([0,1.1])
f3_ax1.set_title("emission PSF")
plt.legend()




fn="C:/research/spin disk/figure new/figure0/PSF_show.svg"
plt.savefig(fn)