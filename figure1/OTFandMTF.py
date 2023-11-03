# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 21:35:50 2022

@author: localadmin
"""
import numpy as np
from scipy.io import loadmat

import matplotlib.pyplot as plt
import scipy.ndimage

fn="C:/research/spin disk/simulation/PSF/OTFandPSF.mat"
mat=loadmat(fn)
OTF_ex=mat['OTF_ex']
OTF_em=mat['OTF_em']
OTF_im=mat['OTF_im']
PSF_ex=mat['PSF_ex']
PSF_em=mat['PSF_em']
PSF_im=mat['PSF_im']


import matplotlib.gridspec as gridspec
fig = plt.figure(figsize=(10,4),constrained_layout=True)
gs = gridspec.GridSpec(ncols=6, nrows=3, wspace=0.0, hspace=0.0, figure=fig)

f3_ax1 = fig.add_subplot(gs[0, 0])
f3_ax1.imshow(abs(OTF_ex)/np.amax(abs(OTF_ex)),cmap="afmhot")
f3_ax1.axis("off")


f3_ax1 = fig.add_subplot(gs[0, 1])
f3_ax1.imshow(abs(OTF_em)/np.amax(abs(OTF_em)),cmap="afmhot")
f3_ax1.axis("off")

f3_ax1 = fig.add_subplot(gs[0, 2])
f3_ax1.imshow(abs(OTF_im)/np.amax(abs(OTF_im)),cmap="afmhot")
f3_ax1.axis("off")

f3_ax1 = fig.add_subplot(gs[1, 0])
f3_ax1.imshow(abs(PSF_ex[256-50:256+50,256-50:256+50])/np.amax(abs(PSF_ex[256-50:256+50,256-50:256+50])),cmap="afmhot")
f3_ax1.axis("off")


f3_ax1 = fig.add_subplot(gs[1, 1])
f3_ax1.imshow(abs(PSF_em[256-50:256+50,256-50:256+50])/np.amax(abs(PSF_em[256-50:256+50,256-50:256+50])),cmap="afmhot")
f3_ax1.axis("off")

f3_ax1 = fig.add_subplot(gs[1, 2])
pc=f3_ax1.imshow(abs(PSF_im[256-50:256+50,256-50:256+50])/np.amax(abs(PSF_im[256-50:256+50,256-50:256+50])),cmap="afmhot")
f3_ax1.axis("off")


fn="C:/research/spin disk/simulation/PSF/OTFandPSF_scalar.mat"
mat=loadmat(fn)
OTF_exs=mat['OTF_em']
OTF_ems=mat['OTF_em']
OTF_ims=mat['OTF_im']
PSF_exs=mat['PSF_em']
PSF_ems=mat['PSF_em']
PSF_ims=mat['PSF_im']

f3_ax1 = fig.add_subplot(gs[0, 3])
f3_ax1.imshow(abs(OTF_exs)/np.amax(abs(OTF_exs)),cmap="afmhot")
f3_ax1.axis("off")


f3_ax1 = fig.add_subplot(gs[0, 4])
f3_ax1.imshow(abs(OTF_ems)/np.amax(abs(OTF_ems)),cmap="afmhot")
f3_ax1.axis("off")

f3_ax1 = fig.add_subplot(gs[0, 5])
f3_ax1.imshow(abs(OTF_ims)/np.amax(abs(OTF_ims)),cmap="afmhot")
f3_ax1.axis("off")

f3_ax1 = fig.add_subplot(gs[1, 3])
f3_ax1.imshow(abs(PSF_exs[256-50:256+50,256-50:256+50])/np.amax(abs(PSF_exs[256-50:256+50,256-50:256+50])),cmap="afmhot")
f3_ax1.axis("off")


f3_ax1 = fig.add_subplot(gs[1, 4])
f3_ax1.imshow(abs(PSF_ems[256-50:256+50,256-50:256+50])/np.amax(abs(PSF_ems[256-50:256+50,256-50:256+50])),cmap="afmhot")
f3_ax1.axis("off")

f3_ax1 = fig.add_subplot(gs[1, 5])
pc=f3_ax1.imshow(abs(PSF_ims[256-50:256+50,256-50:256+50])/np.amax(abs(PSF_ims[256-50:256+50,256-50:256+50])),cmap="afmhot")
f3_ax1.axis("off")


f3_ax1 = fig.add_subplot(gs[2, 0:2])
f3_ax1.plot(np.linspace(-50*5,50*5,100),scipy.ndimage.shift(abs(PSF_ex[256,256-50:256+50])/np.amax(abs(PSF_ex[256,256-50:256+50])),(1)),color="k",label="vectorial")
f3_ax1.plot(np.linspace(-50*5,50*5,100),abs(PSF_exs[256,256-50:256+50])/np.amax(abs(PSF_exs[256,256-50:256+50])),color="red",label="Scalar")

f3_ax1 = fig.add_subplot(gs[2, 2:4])
f3_ax1.plot(np.linspace(-70*5,70*5,140),scipy.ndimage.shift(abs(PSF_em[256,256-70:256+70])/np.amax(abs(PSF_em[256,256-70:256+70])),(-0.5)),color="k",label="vectorial")
f3_ax1.plot(np.linspace(-70*5,70*5,140),abs(PSF_ems[256,256-70:256+70])/np.amax(abs(PSF_ems[256,256-70:256+70])),color="red",label="Scalar")

f3_ax1 = fig.add_subplot(gs[2, 4:])
f3_ax1.plot(np.linspace(-50*5,50*5,100),scipy.ndimage.shift(abs(PSF_im[256,256-50:256+50])/np.amax(abs(PSF_im[256,256-50:256+50])),(-0.5)),color="k",label="Vectorial")
f3_ax1.plot(np.linspace(-50*5,50*5,100),abs(PSF_ims[256,256-50:256+50])/np.amax(abs(PSF_ims[256,256-50:256+50])),color="red",label="Scalar")
plt.legend()
fn="C:/research/spin disk/figure/OTFandPSF/PSFandOTF.svg"
#plt.savefig(fn)