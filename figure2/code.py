# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 18:04:17 2022

@author: localadmin
"""
import numpy as np
from scipy.io import loadmat

from skimage import io
import matplotlib.pyplot as plt
import scipy.ndimage

fn="C:/research/spin disk/figure new/figure2/new/PSFs_width_oversizeof pinhole.mat"
mat=loadmat(fn)
PSF_im_fitting_arr=abs(mat['PSF_im_fitting_arr'][0])
PSF_em_fitting_arr=abs(mat['PSF_em_fitting_arr'][0])
pinhole_scale=mat['pinhole_scale'][0]

fn="C:/research/spin disk/figure new/figure2/new/PSFs_width_oversizeof pinhole_vect.mat"
mat=loadmat(fn)
PSF_im_fitting_arr_vect=abs(mat['PSF_im_fitting_arr'][0])
PSF_em_fitting_arr_vect=abs(mat['PSF_em_fitting_arr'][0])
#pinhole_scale=mat['pinhole_scale'][0]


fn="C:/research/spin disk/figure new/figure2/new/confocal_PSFs_width_oversizeof pinhole_vect.mat";
mat=loadmat(fn)
PSF_conf_fitting_arr_vect=abs(mat['PSF_conf_fitting_arr'][0])
fn="C:/research/spin disk/figure new/figure2/new/confocal_PSFs_width_oversizeof pinhole_scalar.mat";
mat=loadmat(fn)
PSF_conf_fitting_arr_scalar=abs(mat['PSF_conf_fitting_arr'][0])

import matplotlib.gridspec as gridspec
markersize=4

fig = plt.figure(figsize=(10,10.0),constrained_layout=True)
gs = gridspec.GridSpec(ncols=3, nrows=3, wspace=0.1, hspace=0.0, figure=fig)

f3_ax1 = fig.add_subplot(gs[0, 0])

#f3_ax1.plot(pinhole_scale,PSF_im_fitting_arr*2.355,"-",linewidth=3,color="k",label="vectorial")
#f3_ax1.plot(pinhole_scale,PSF_im_fitting_arr*2.355,"-",linewidth=3,color="r",label="scalar")

f3_ax1.plot(pinhole_scale,PSF_im_fitting_arr*2.355,"o",linewidth=3,markersize=markersize,color="r")
f3_ax1.plot(pinhole_scale,PSF_em_fitting_arr*2.355,"--",linewidth=3,color="r")
f3_ax1.plot(pinhole_scale,PSF_conf_fitting_arr_scalar*2.355,"s",markersize=markersize,linewidth=3,color="r")

f3_ax1.plot(pinhole_scale,PSF_im_fitting_arr_vect*2.355,"o",linewidth=3,markersize=markersize,color="k",label="ISM PSF")
f3_ax1.plot(pinhole_scale,PSF_em_fitting_arr_vect*2.355,"--",linewidth=3,color="k",label="WF PSF")

f3_ax1.plot(pinhole_scale,PSF_conf_fitting_arr_vect*2.355,"s",linewidth=3,markersize=markersize,color="k",label="confocal PSF")

f3_ax1.legend()

f3_ax1.set_xlabel("pinhole size (A.U.)")
f3_ax1.set_ylabel("FWHM of PSF (nm)")


f3_ax1 = fig.add_subplot(gs[0, 1])
f3_ax1.plot(pinhole_scale,PSF_em_fitting_arr/PSF_im_fitting_arr,"-",linewidth=3,markersize=markersize,color="r",label="scalar")
f3_ax1.plot(pinhole_scale,PSF_em_fitting_arr_vect/PSF_im_fitting_arr_vect,"-",linewidth=3,markersize=markersize,color="k",label="vectorial")
f3_ax1.plot(pinhole_scale,np.ones([31])*np.sqrt(2),"b-",linewidth=2,label=r"therotical expectation $\sqrt{2}$")
f3_ax1.set_ylim([1.,1.7])
f3_ax1.set_ylabel("Improvement factor (ISM-WF)")
f3_ax1.set_xlabel("pinhole size (A.U.)")
f3_ax1.legend()


f3_ax1 = fig.add_subplot(gs[0, 2])
f3_ax1.plot(pinhole_scale,PSF_conf_fitting_arr_scalar/PSF_im_fitting_arr,"-",markersize=markersize,linewidth=3,color="r",label="scalar")
f3_ax1.plot(pinhole_scale,PSF_conf_fitting_arr_vect/PSF_im_fitting_arr_vect,"-",markersize=markersize,linewidth=3,color="k",label="vectorial")
f3_ax1.set_ylim([1.,1.7])

print(np.amax(PSF_conf_fitting_arr_scalar/PSF_im_fitting_arr))
print(np.amax(PSF_conf_fitting_arr_vect/PSF_im_fitting_arr_vect))

f3_ax1.plot(pinhole_scale,np.ones([31])*np.sqrt(2),"b-",linewidth=2,label=r"therotical expectation $\sqrt{2}$")
f3_ax1.legend()
#plt.legend()
f3_ax1.set_xlabel("pinhole size (A.U.)")
f3_ax1.set_ylabel("Improvement factor (ISM-confocal)")

if False:
    fn="C:/research/spin disk/figure new/figure2/new/overdepth/PSFs_width_oversizeof pinhole.mat"
    mat=loadmat(fn)
    PSF_im_fitting_arr=abs(mat['PSF_im_fitting_arr'][0])
    PSF_em_fitting_arr=abs(mat['PSF_em_fitting_arr'][0])
    pinhole_scale=np.linspace(-200,200,31)

    fn="C:/research/spin disk/figure new/figure2/new/overdepth/PSFs_width_oversizeof pinhole_vect.mat"
    mat=loadmat(fn)
    PSF_im_fitting_arr_vect=abs(mat['PSF_im_fitting_arr'][0])
    PSF_em_fitting_arr_vect=abs(mat['PSF_em_fitting_arr'][0])
    #pinhole_scale=mat['pinhole_scale'][0]


    fn="C:/research/spin disk/figure new/figure2/new/overdepth/confocal_PSFs_width_oversizeof pinhole_vect.mat";
    mat=loadmat(fn)
    PSF_conf_fitting_arr_vect=abs(mat['PSF_conf_fitting_arr'][0])
    fn="C:/research/spin disk/figure new/figure2/new/overdepth/confocal_PSFs_width_oversizeof pinhole_scalar.mat";
    mat=loadmat(fn)
    PSF_conf_fitting_arr_scalar=abs(mat['PSF_conf_fitting_arr'][0])


    f3_ax1 = fig.add_subplot(gs[1, 0])

    f3_ax1.plot(pinhole_scale,PSF_im_fitting_arr*2.355,"-",linewidth=3,color="k",label="vectorial")
    f3_ax1.plot(pinhole_scale,PSF_im_fitting_arr*2.355,"-",linewidth=3,color="r",label="scalar")

    f3_ax1.plot(pinhole_scale,PSF_im_fitting_arr*2.355,"o-",linewidth=3,markersize=10,color="r")
    f3_ax1.plot(pinhole_scale,PSF_em_fitting_arr*2.355,"--",linewidth=3,color="r")
    f3_ax1.plot(pinhole_scale,PSF_conf_fitting_arr_scalar*2.355,"s-",markersize=10,linewidth=3,color="r")

    f3_ax1.plot(pinhole_scale,PSF_im_fitting_arr_vect*2.355,"o-",linewidth=3,markersize=10,color="k",label="ISM PSF")
    f3_ax1.plot(pinhole_scale,PSF_em_fitting_arr_vect*2.355,"--",linewidth=3,color="k",label="WF PSF")

    f3_ax1.plot(pinhole_scale,PSF_conf_fitting_arr_vect*2.355,"s-",linewidth=3,markersize=10,color="k",label="confocal PSF")

    #f3_ax1.legend()

    f3_ax1.set_xlabel("pinhole size (A.U.)")
    f3_ax1.set_ylabel("FWHM of PSF (nm)")


    f3_ax1 = fig.add_subplot(gs[1, 1])
    f3_ax1.plot(pinhole_scale,PSF_em_fitting_arr/PSF_im_fitting_arr,"-",linewidth=3,markersize=10,color="r",label="scalar")
    f3_ax1.plot(pinhole_scale,PSF_em_fitting_arr_vect/PSF_im_fitting_arr_vect,"-",linewidth=3,markersize=10,color="k",label="vectorial")
    f3_ax1.plot(pinhole_scale,np.ones([31])*np.sqrt(2),"b-",linewidth=2,label=r"therotical expectation $\sqrt{2}$")
    #f3_ax1.set_ylim([1.,1.7])
    f3_ax1.set_ylabel("Improvement factor (ISM-WF)")
    f3_ax1.set_xlabel("pinhole size (A.U.)")
    #f3_ax1.legend()

    f3_ax1 = fig.add_subplot(gs[1, 2])
    f3_ax1.plot(pinhole_scale,PSF_conf_fitting_arr_scalar/PSF_im_fitting_arr,"-",markersize=10,linewidth=3,color="r",label="scalar")
    f3_ax1.plot(pinhole_scale,PSF_conf_fitting_arr_vect/PSF_im_fitting_arr_vect,"-",markersize=10,linewidth=3,color="k",label="vectorial")
    #f3_ax1.set_ylim([1.,1.7])

    f3_ax1.plot(pinhole_scale,np.ones([31])*np.sqrt(2),"b-",linewidth=2,label=r"therotical expectation $\sqrt{2}$")
    #f3_ax1.legend()
    #plt.legend()
    f3_ax1.set_xlabel("pinhole size (A.U.)")
    f3_ax1.set_ylabel("Improvement factor (ISM-confocal)")


if True:
    fn="C:/research/spin disk/figure new/figure2/new/overswipfactor/PSFs_width_oversizeof pinhole.mat"
    mat=loadmat(fn)
    PSF_im_fitting_arr=abs(mat['PSF_im_fitting_arr'][0])
    PSF_em_fitting_arr=abs(mat['PSF_em_fitting_arr'][0])
    pinhole_scale=np.linspace(1.0,3.0,31)

    fn="C:/research/spin disk/figure new/figure2/new/overswipfactor/PSFs_width_oversizeof pinhole_vect.mat"
    mat=loadmat(fn)
    PSF_im_fitting_arr_vect=abs(mat['PSF_im_fitting_arr'][0])
    PSF_em_fitting_arr_vect=abs(mat['PSF_em_fitting_arr'][0])
    #pinhole_scale=mat['pinhole_scale'][0]


    f3_ax1 = fig.add_subplot(gs[1, 0])

    #f3_ax1.plot(pinhole_scale,PSF_im_fitting_arr*2.355,"-",linewidth=3,color="k",label="vectorial")
    #f3_ax1.plot(pinhole_scale,PSF_im_fitting_arr*2.355,"-",linewidth=3,color="r",label="scalar")

    f3_ax1.plot(pinhole_scale,PSF_im_fitting_arr*2.355,"o",linewidth=3,markersize=markersize,color="r")
    f3_ax1.plot(pinhole_scale,PSF_em_fitting_arr*2.355,"--",linewidth=3,color="r")

    f3_ax1.plot(pinhole_scale,PSF_im_fitting_arr_vect*2.355,"o",linewidth=3,markersize=markersize,color="k",label="ISM PSF")
    f3_ax1.plot(pinhole_scale,PSF_em_fitting_arr_vect*2.355,"--",linewidth=3,color="k",label="WF PSF")


    f3_ax1.legend()

    f3_ax1.set_xlabel("swiping factor")
    f3_ax1.set_ylabel("FWHM of PSF (nm)")


    f3_ax1 = fig.add_subplot(gs[1, 1])
    f3_ax1.plot(pinhole_scale,PSF_em_fitting_arr/PSF_im_fitting_arr,"-",linewidth=3,markersize=10,color="r",label="scalar")
    f3_ax1.plot(pinhole_scale,PSF_em_fitting_arr_vect/PSF_im_fitting_arr_vect,"-",linewidth=3,markersize=10,color="k",label="vectorial")
    
    ind_vec=np.argmax(PSF_em_fitting_arr_vect/PSF_im_fitting_arr_vect)
    print("best sweep factor vec PSF:"+str(pinhole_scale[ind_vec]))
    
    ind_vec=np.argmax(PSF_em_fitting_arr/PSF_im_fitting_arr)
    print("best sweep factor SCALAR psf:"+str(pinhole_scale[ind_vec]))
    
    f3_ax1.plot(pinhole_scale,np.ones([31])*np.sqrt(2),"b-",linewidth=2,label=r"therotical expectation $\sqrt{2}$")
    #f3_ax1.set_ylim([1.,1.7])
    f3_ax1.set_ylabel("Improvement factor (ISM-WF)")
    f3_ax1.set_xlabel("swiping factor")
    f3_ax1.legend()

    f3_ax1 = fig.add_subplot(gs[1, 2])
    f3_ax1.plot(pinhole_scale,PSF_im_fitting_arr[0]/PSF_im_fitting_arr,"-",markersize=markersize,linewidth=3,color="r",label="scalar")
    f3_ax1.plot(pinhole_scale,PSF_im_fitting_arr_vect[0]/PSF_im_fitting_arr_vect,"-",markersize=markersize,linewidth=3,color="k",label="vectorial")
    #f3_ax1.set_ylim([1.,1.7])
    
    ind_vec=np.argmax(PSF_im_fitting_arr_vect[0]/PSF_im_fitting_arr_vect)
    print("best sweep factor vec PSF:"+str(pinhole_scale[ind_vec]))
    
    ind_vec=np.argmax(PSF_im_fitting_arr[0]/PSF_im_fitting_arr)
    print("best sweep factor SCALAR psf:"+str(pinhole_scale[ind_vec]))
    
    f3_ax1.plot(pinhole_scale,np.ones([31])*np.sqrt(2),"b-",linewidth=2,label=r"therotical expectation $\sqrt{2}$")
    f3_ax1.legend()
    #plt.legend()
    f3_ax1.set_xlabel("swiping factor")
    f3_ax1.set_ylabel("Improvement factor (ISM-confocal)")


fn="C:/research/spin disk/figure new/figure2/new/overdepth1/PSFs_width_oversizeof pinhole.mat"
mat=loadmat(fn)
PSF_im_fitting_arr=abs(mat['PSF_im_fitting_arr'][0])*0.1
PSF_em_fitting_arr=abs(mat['PSF_em_fitting_arr'][0])*0.1
pinhole_scale=np.linspace(0,400,15)

fn="C:/research/spin disk/figure new/figure2/new/overdepth1/PSFs_width_oversizeof pinhole_vect.mat"
mat=loadmat(fn)
PSF_im_fitting_arr_vect=abs(mat['PSF_im_fitting_arr'][0])*0.1
PSF_em_fitting_arr_vect=abs(mat['PSF_em_fitting_arr'][0])*0.1
#pinhole_scale=mat['pinhole_scale'][0]


fn="C:/research/spin disk/figure new/figure2/new/overdepth1/confocal_PSFs_width_oversizeof pinhole_vect.mat";
mat=loadmat(fn)
PSF_conf_fitting_arr_vect=abs(mat['PSF_conf_fitting_arr'][0])*0.1
fn="C:/research/spin disk/figure new/figure2/new/overdepth1/confocal_PSFs_width_oversizeof pinhole_scalar.mat";
mat=loadmat(fn)
PSF_conf_fitting_arr_scalar=abs(mat['PSF_conf_fitting_arr'][0])*0.1


f3_ax1 = fig.add_subplot(gs[2, 0])
f3_ax1.plot(pinhole_scale,PSF_im_fitting_arr,"-",linewidth=3,color="k",label="vectorial")
f3_ax1.plot(pinhole_scale,PSF_im_fitting_arr,"-",linewidth=3,color="r",label="scalar")

f3_ax1.plot(pinhole_scale,PSF_im_fitting_arr,"o-",linewidth=3,markersize=10,color="r")
f3_ax1.plot(pinhole_scale,PSF_em_fitting_arr,"--",linewidth=3,color="r")
f3_ax1.plot(pinhole_scale,PSF_conf_fitting_arr_scalar,"s-",markersize=10,linewidth=3,color="r")

f3_ax1.plot(pinhole_scale,PSF_im_fitting_arr_vect,"o-",linewidth=3,markersize=10,color="k",label="ISM PSF")
f3_ax1.plot(pinhole_scale,PSF_em_fitting_arr_vect,"--",linewidth=3,color="k",label="WF PSF")

f3_ax1.plot(pinhole_scale,PSF_conf_fitting_arr_vect,"s-",linewidth=3,markersize=10,color="k",label="confocal PSF")

#f3_ax1.legend()

f3_ax1.set_xlabel("pinhole size (A.U.)")
f3_ax1.set_ylabel("FWHM of PSF (nm)")


f3_ax1 = fig.add_subplot(gs[2, 1])
f3_ax1.plot(pinhole_scale,PSF_em_fitting_arr/PSF_im_fitting_arr,"-",linewidth=3,markersize=10,color="r",label="scalar")
f3_ax1.plot(pinhole_scale,PSF_em_fitting_arr_vect/PSF_im_fitting_arr_vect,"-",linewidth=3,markersize=10,color="k",label="vectorial")
f3_ax1.plot(pinhole_scale,np.ones([15])*np.sqrt(2),"b-",linewidth=2,label=r"therotical expectation $\sqrt{2}$")
#f3_ax1.set_ylim([1.,1.7])
f3_ax1.set_ylabel("Improvement factor (ISM-WF)")
f3_ax1.set_xlabel("pinhole size (A.U.)")
#f3_ax1.legend()

f3_ax1 = fig.add_subplot(gs[2, 2])
f3_ax1.plot(pinhole_scale,PSF_conf_fitting_arr_scalar/PSF_im_fitting_arr,"-",markersize=10,linewidth=3,color="r",label="scalar")
f3_ax1.plot(pinhole_scale,PSF_conf_fitting_arr_vect/PSF_im_fitting_arr_vect,"-",markersize=10,linewidth=3,color="k",label="vectorial")
#f3_ax1.set_ylim([1.,1.7])

f3_ax1.plot(pinhole_scale,np.ones([15])*np.sqrt(2),"b-",linewidth=2,label=r"therotical expectation $\sqrt{2}$")
#f3_ax1.legend()
#plt.legend()
f3_ax1.set_xlabel("pinhole size (A.U.)")
f3_ax1.set_ylabel("Improvement factor (ISM-confocal)")

fn="C:/research/spin disk/figure new/figure2/overpinhole_improvement.svg"
#plt.savefig(fn)
plt.show()