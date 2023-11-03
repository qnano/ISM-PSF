# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 13:45:10 2022

@author: localadmin
"""
import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt



import matplotlib.gridspec as gridspec
fig = plt.figure(figsize=(12,8),constrained_layout=True)
gs = gridspec.GridSpec(ncols=2, nrows=2, wspace=0.0, hspace=0.0, figure=fig)


f3_ax1 = fig.add_subplot(gs[1, 0])
f3_ax1.set_xlabel("pinhole size (A.U.)")
f3_ax1.set_ylabel("PSF width (nm)")

pinhole_scale=np.linspace(0.3,2,11);
color1="r"
color2="navy"
color3="teal"
color4="yellowgreen"

fn="C:/research/spin disk/figure new/figure6/0ml_imaging_scanning.mat"
mat=loadmat(fn)
PSF_im_fitting_arr=mat['PSF_im_fitting_arr']*2.355
PSF_em_fitting_arr=mat['PSF_em_fitting_arr']*2.355
f3_ax1.plot(pinhole_scale,np.mean(PSF_im_fitting_arr,axis=0),"-",linewidth=3,color=color1,label=r"one photon")
f3_ax1.fill_between(pinhole_scale,np.mean(PSF_im_fitting_arr,axis=0)-np.std(PSF_im_fitting_arr,axis=0),np.mean(PSF_im_fitting_arr,axis=0)+np.std(PSF_im_fitting_arr,axis=0),color=color1,alpha=0.3)


fn="C:/research/spin disk/figure new/figure6/two_photon_PSF_width_imaging.mat"
mat=loadmat(fn)
PSF_im_fitting_arr=mat['PSF_im_fitting_arr']*2.355
PSF_em_fitting_arr=mat['PSF_em_fitting_arr']*2.355
f3_ax1.plot(pinhole_scale,np.mean(PSF_im_fitting_arr,axis=0),"-",linewidth=3,color=color2,label=r"two photon")
f3_ax1.fill_between(pinhole_scale,np.mean(PSF_im_fitting_arr,axis=0)-np.std(PSF_im_fitting_arr,axis=0),np.mean(PSF_im_fitting_arr,axis=0)+np.std(PSF_im_fitting_arr,axis=0),color=color2,alpha=0.3)

fn="C:/research/spin disk/figure new/figure6/two_photon_PSF_width_imaging_optimal.mat"
mat=loadmat(fn)
PSF_im_fitting_arr=mat['PSF_im_fitting_arr']*2.355
PSF_em_fitting_arr=mat['PSF_em_fitting_arr']*2.355
f3_ax1.plot(pinhole_scale,np.mean(PSF_im_fitting_arr,axis=0),"*-",linewidth=3,color=color2,label=r"optimal two photon")
f3_ax1.fill_between(pinhole_scale,np.mean(PSF_im_fitting_arr,axis=0)-np.std(PSF_im_fitting_arr,axis=0),np.mean(PSF_im_fitting_arr,axis=0)+np.std(PSF_im_fitting_arr,axis=0),color=color2,alpha=0.3)

print("FWHM of optimal 2 P ISM:"+str(np.mean(PSF_im_fitting_arr,axis=0)[4]))
fn="C:/research/spin disk/figure new/figure6/three_photon_PSF_width_imaging.mat"
mat=loadmat(fn)
PSF_im_fitting_arr=mat['PSF_im_fitting_arr']*2.355
PSF_em_fitting_arr=mat['PSF_em_fitting_arr']*2.355
f3_ax1.plot(pinhole_scale,np.mean(PSF_im_fitting_arr,axis=0),"-",linewidth=3,color=color3,label=r"three photon")
f3_ax1.fill_between(pinhole_scale,np.mean(PSF_im_fitting_arr,axis=0)-np.std(PSF_im_fitting_arr,axis=0),np.mean(PSF_im_fitting_arr,axis=0)+np.std(PSF_im_fitting_arr,axis=0),color=color3,alpha=0.3)

fn="C:/research/spin disk/figure new/figure6/three_photon_PSF_width_imaging_optimal.mat"
mat=loadmat(fn)
PSF_im_fitting_arr=mat['PSF_im_fitting_arr']*2.355
PSF_em_fitting_arr=mat['PSF_em_fitting_arr']*2.355
f3_ax1.plot(pinhole_scale,np.mean(PSF_im_fitting_arr,axis=0),"*-",linewidth=3,color=color3,label=r"optimal three photon")
f3_ax1.fill_between(pinhole_scale,np.mean(PSF_im_fitting_arr,axis=0)-np.std(PSF_im_fitting_arr,axis=0),np.mean(PSF_im_fitting_arr,axis=0)+np.std(PSF_im_fitting_arr,axis=0),color=color2,alpha=0.3)

print("FWHM of optimal 3 P ISM:"+str(np.mean(PSF_im_fitting_arr,axis=0)[4]))



fn="C:/research/spin disk/figure new/figure6/0ml_imaging_scanning.mat"
mat=loadmat(fn)
PSF_im_fitting_arr=mat['PSF_im_fitting_arr']*2.355
PSF_em_fitting_arr=mat['PSF_em_fitting_arr']*2.355
f3_ax1.plot(pinhole_scale,np.mean(PSF_em_fitting_arr,axis=0),"-",linewidth=3,color=color4,label=r"WF PSF")
f3_ax1.fill_between(pinhole_scale,np.mean(PSF_em_fitting_arr,axis=0)-np.std(PSF_em_fitting_arr,axis=0),np.mean(PSF_em_fitting_arr,axis=0)+np.std(PSF_em_fitting_arr,axis=0),color=color4,alpha=0.3)
plt.legend()


f3_ax1 = fig.add_subplot(gs[1, 1])
f3_ax1.set_xlabel("pinhole size (A.U.)")
f3_ax1.set_ylabel("improvement factor")

pinhole_scale=np.linspace(0.3,2,11);
color1="r"
color2="navy"
color3="teal"
color4="cyan"

fn="C:/research/spin disk/figure new/figure6/0ml_imaging_scanning.mat"
mat=loadmat(fn)
PSF_im_fitting_arr=mat['PSF_im_fitting_arr']*2.355
PSF_em_fitting_arr=mat['PSF_em_fitting_arr']*2.355
f3_ax1.plot(pinhole_scale,np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0),"-",linewidth=3,color=color1,label=r"1-photon ex. ISM")
f3_ax1.fill_between(pinhole_scale,np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0)-np.std(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0),np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0)+np.std(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0),color=color1,alpha=0.3)


fn="C:/research/spin disk/figure new/figure6/two_photon_PSF_width_imaging.mat"
mat=loadmat(fn)
PSF_im_fitting_arr=mat['PSF_im_fitting_arr']*2.355

PSF_em_fitting_arr=mat['PSF_em_fitting_arr']*2.355
f3_ax1.plot(pinhole_scale,np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0),"-",linewidth=3,color=color2,label=r"2-photon ex. ISM")
f3_ax1.fill_between(pinhole_scale,np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0)-np.std(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0),np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0)+np.std(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0),color=color2,alpha=0.3)

print("best improvement factor of 2P ISM")
print(np.amax(np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0)))

fn="C:/research/spin disk/figure new/figure6/two_photon_PSF_width_imaging_optimal.mat"
mat=loadmat(fn)
PSF_im_fitting_arr=mat['PSF_im_fitting_arr']*2.355

PSF_em_fitting_arr=mat['PSF_em_fitting_arr']*2.355
f3_ax1.plot(pinhole_scale,np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0),"*-",linewidth=3,color=color2,label=r"optimal 2-photon ex. ISM")
f3_ax1.fill_between(pinhole_scale,np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0)-np.std(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0),np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0)+np.std(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0),color=color2,alpha=0.3)

print("best improvement factor of optimal 2P ISM")
print(np.amax(np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0)))

fn="C:/research/spin disk/figure new/figure6/three_photon_PSF_width_imaging.mat"
mat=loadmat(fn)
PSF_im_fitting_arr=mat['PSF_im_fitting_arr']*2.355

PSF_em_fitting_arr=mat['PSF_em_fitting_arr']*2.355
f3_ax1.plot(pinhole_scale,np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0),"-",linewidth=3,color=color3,label=r"3-photon ex. ISM")
f3_ax1.fill_between(pinhole_scale,np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0)-np.std(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0),np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0)+np.std(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0),color=color3,alpha=0.3)

print("best improvement factor of 3P ISM")
print(np.amax(np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0)))

fn="C:/research/spin disk/figure new/figure6/three_photon_PSF_width_imaging_optimal.mat"
mat=loadmat(fn)
PSF_im_fitting_arr=mat['PSF_im_fitting_arr']*2.355

PSF_em_fitting_arr=mat['PSF_em_fitting_arr']*2.355
f3_ax1.plot(pinhole_scale,np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0),"*-",linewidth=3,color=color3,label=r"optimal 3-photon ex. ISM")
f3_ax1.fill_between(pinhole_scale,np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0)-np.std(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0),np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0)+np.std(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0),color=color2,alpha=0.3)

print("best improvement factor of optimal 3P ISM")
print(np.amax(np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0)))

f3_ax1.plot(pinhole_scale,np.ones([11])*np.sqrt(2),"k--",label=r"therotical expectation $\sqrt{2}$")
plt.legend()


f3_ax1 = fig.add_subplot(gs[0, 0])
f3_ax1.set_xlabel("swipe factor")
f3_ax1.set_ylabel("PSF width (nm)")

pinhole_scale=np.linspace(1.0,3.0,31);
color1="r"
color2="navy"
color3="teal"
color4="yellowgreen"

fn="C:/research/spin disk/figure new/figure6/one_photon_PSF_width_imaging_swipfactor.mat"
mat=loadmat(fn)
PSF_im_fitting_arr=abs(mat['PSF_im_fitting_arr'])*2.355
PSF_em_fitting_arr=abs(mat['PSF_em_fitting_arr'])*2.355
f3_ax1.plot(pinhole_scale,np.mean(PSF_im_fitting_arr,axis=0),"-",linewidth=3,color=color1,label=r"one photon")
f3_ax1.fill_between(pinhole_scale,np.mean(PSF_im_fitting_arr,axis=0)-np.std(PSF_im_fitting_arr,axis=0),np.mean(PSF_im_fitting_arr,axis=0)+np.std(PSF_im_fitting_arr,axis=0),color=color1,alpha=0.3)


fn="C:/research/spin disk/figure new/figure6/two_photon_PSF_width_imaging_swipfactor.mat"
mat=loadmat(fn)
PSF_im_fitting_arr=abs(mat['PSF_im_fitting_arr'])*2.355
PSF_em_fitting_arr=abs(mat['PSF_em_fitting_arr'])*2.355
f3_ax1.plot(pinhole_scale,np.mean(PSF_im_fitting_arr,axis=0),"-",linewidth=3,color=color2,label=r"two photon")
f3_ax1.fill_between(pinhole_scale,np.mean(PSF_im_fitting_arr,axis=0)-np.std(PSF_im_fitting_arr,axis=0),np.mean(PSF_im_fitting_arr,axis=0)+np.std(PSF_im_fitting_arr,axis=0),color=color2,alpha=0.3)

ind_two=np.argmax(np.mean(PSF_im_fitting_arr,axis=0))
print("optimal sweep of 2 P ISM"+str(pinhole_scale[ind_two]))

fn="C:/research/spin disk/figure new/figure6/three_photon_PSF_width_imaging_swipfactor.mat"
mat=loadmat(fn)
PSF_im_fitting_arr=abs(mat['PSF_im_fitting_arr'])*2.355
PSF_em_fitting_arr=abs(mat['PSF_em_fitting_arr'])*2.355
f3_ax1.plot(pinhole_scale,np.mean(PSF_im_fitting_arr,axis=0),"-",linewidth=3,color=color3,label=r"three photon")
f3_ax1.fill_between(pinhole_scale,np.mean(PSF_im_fitting_arr,axis=0)-np.std(PSF_im_fitting_arr,axis=0),np.mean(PSF_im_fitting_arr,axis=0)+np.std(PSF_im_fitting_arr,axis=0),color=color3,alpha=0.3)

ind_three=np.argmax(np.mean(PSF_im_fitting_arr,axis=0))
print("optimal sweep of 3 P ISM"+str(pinhole_scale[ind_three]))


fn="C:/research/spin disk/figure new/figure6/one_photon_PSF_width_imaging_swipfactor.mat"
mat=loadmat(fn)
PSF_im_fitting_arr=abs(mat['PSF_im_fitting_arr'])*2.355
PSF_em_fitting_arr=abs(mat['PSF_em_fitting_arr'])*2.355
f3_ax1.plot(pinhole_scale,np.mean(PSF_em_fitting_arr,axis=0),"-",linewidth=3,color=color4,label=r"WF PSF")
f3_ax1.fill_between(pinhole_scale,np.mean(PSF_em_fitting_arr,axis=0)-np.std(PSF_em_fitting_arr,axis=0),np.mean(PSF_em_fitting_arr,axis=0)+np.std(PSF_em_fitting_arr,axis=0),color=color4,alpha=0.3)
plt.legend()


f3_ax1 = fig.add_subplot(gs[0, 1])
f3_ax1.set_xlabel("swipe factor")
f3_ax1.set_ylabel("improvement factor")

color1="r"
color2="navy"
color3="teal"
color4="cyan"

fn="C:/research/spin disk/figure new/figure6/one_photon_PSF_width_imaging_swipfactor.mat"
mat=loadmat(fn)
PSF_im_fitting_arr=abs(mat['PSF_im_fitting_arr'])*2.355
PSF_em_fitting_arr=abs(mat['PSF_em_fitting_arr'])*2.355
f3_ax1.plot(pinhole_scale,np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0),"-",linewidth=3,color=color1,label=r"1-photon ex. ISM")
f3_ax1.fill_between(pinhole_scale,np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0)-np.std(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0),np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0)+np.std(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0),color=color1,alpha=0.3)


fn="C:/research/spin disk/figure new/figure6/two_photon_PSF_width_imaging_swipfactor.mat"
mat=loadmat(fn)
PSF_im_fitting_arr=abs(mat['PSF_im_fitting_arr'])*2.355
PSF_em_fitting_arr=abs(mat['PSF_em_fitting_arr'])*2.355
f3_ax1.plot(pinhole_scale,np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0),"-",linewidth=3,color=color2,label=r"2-photon ex. ISM")
f3_ax1.fill_between(pinhole_scale,np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0)-np.std(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0),np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0)+np.std(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0),color=color2,alpha=0.3)
print("bets swipt factor: "+str(pinhole_scale[np.argmax((PSF_em_fitting_arr/PSF_im_fitting_arr)[0])]))


fn="C:/research/spin disk/figure new/figure6/three_photon_PSF_width_imaging_swipfactor.mat"
mat=loadmat(fn)
PSF_im_fitting_arr=abs(mat['PSF_im_fitting_arr'])*2.355
PSF_em_fitting_arr=abs(mat['PSF_em_fitting_arr'])*2.355
f3_ax1.plot(pinhole_scale,np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0),"-",linewidth=3,color=color3,label=r"3-photon ex. ISM")
f3_ax1.fill_between(pinhole_scale,np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0)-np.std(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0),np.mean(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0)+np.std(PSF_em_fitting_arr/PSF_im_fitting_arr,axis=0),color=color3,alpha=0.3)
f3_ax1.plot(pinhole_scale,np.ones([31])*np.sqrt(2),"k--",label=r"therotical expectation $\sqrt{2}$")
print("bets swipt factor: "+str(pinhole_scale[np.argmax((PSF_em_fitting_arr/PSF_im_fitting_arr)[0])]))

plt.legend()


fn="C:/research/spin disk/figure new/figure6/multi_photon.svg"
#plt.savefig(fn)
