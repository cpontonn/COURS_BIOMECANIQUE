# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 10:58:20 2022

@author: cpontonn
"""
import numpy as np
# import matplotlib.pyplot as plt
from scipy import signal
import ezc3d

def Rp(markers):
    #pelvis-> world rotation matrix


    
    zp=(markers[21,:]-markers[10,:])/np.linalg.norm(markers[21,:]-markers[10,:]) # RFWT-LFWT
    xp=(0.5*(markers[21,:]+markers[10,:])-0.5*(markers[22,:]+markers[11,:]))/np.linalg.norm(0.5*(markers[21,:]+markers[10,:])-0.5*(markers[22,:]+markers[11,:]))
    yp=np.cross(zp,xp)
    xp=np.cross(yp,zp) 
    Rp=np.zeros((3,3))
    Rp[0,:]=[xp[0],yp[0],zp[0]]
    Rp[1,:]=[xp[1],yp[1],zp[1]]
    Rp[2,:]=[xp[2],yp[2],zp[2]]
    return Rp

def rotz(q):
    #z axis rotation matrix of angle q
    y=np.zeros((3,3))

    y[0,:]=[np.cos(q),-np.sin(q),0]
    y[1,:]=[np.sin(q),np.cos(q),0]
    y[2,:]=[0,0,1]
    return y

def LP_filter(x,f_cut,f_sam):
    # low pass filter, no phase shift
    b, a = signal.butter(4, f_cut, btype='low', analog=False, output='ba', fs=f_sam)
    y = signal.filtfilt(b, a, x, padlen=150)
    return y

def LP_filter_vec3(x,f_cut,f_sam):
    y=np.zeros((3,len(x[0,:])))
    y[0,:]=np.transpose(LP_filter(np.transpose(x[0,:]),f_cut,f_sam))
    y[1,:]=np.transpose(LP_filter(np.transpose(x[1,:]),f_cut,f_sam))
    y[2,:]=np.transpose(LP_filter(np.transpose(x[2,:]),f_cut,f_sam))
    
    return y

def diff2(x,dt):
    #2nd order finite difference for derivation
    y=np.zeros(len(x))
    n=len(x)-1
    y[0]=(-x[2]+4*x[1]-3*x[0])/(2*dt) 
    y[n]=(3*x[n]-4*x[n-1]+x[n-2])/(2*dt) 
    for i in range(1,n-1):
        y[i]=(x[i+1]-x[i-1])/(2*dt) #différence centrée d'ordre 2
    
    return y

def diff_vec3(x,dt):
    y=np.zeros((3,len(x[0,:])))
    y[0,:]=diff2(x[0,:],dt)
    y[1,:]=diff2(x[1,:],dt)
    y[2,:]=diff2(x[2,:],dt)
    
    return y

def generate_GRF(pf):
  
    F=np.zeros((len(pf['force'][:,0]),len(pf['force'][0,:])))#force
    P=np.zeros((len(pf['force'][:,0]),len(pf['force'][0,:])))#centre of pressure
    for f in range(0,len(pf['force'][0,:])-1):

        F[:,f]=pf['force'][:,f]# force in mocap coordinate system

        P[:,f]=pf['center_of_pressure'][:,f]/1000 # CoP in mocap coordinate system
    return F,P

def replace_nan(sig,m1):
    
    for i in range(0,len(sig[0,:])):
        if np.isnan(sig[0,i]):
            sig[0,i]=m1
        if np.isinf(sig[0,i]):
            sig[0,i]=m1   
        if np.isnan(sig[1,i]):
            sig[1,i]=m1
        if np.isinf(sig[1,i]):
            sig[1,i]=m1
        if np.isnan(sig[2,i]):
            sig[2,i]=m1
        if np.isinf(sig[2,i]):
            sig[2,i]=m1
    return sig

def resample_vec3(sig,nf):
    sig2=np.zeros((3,nf))
    sig2[0,:]=signal.resample(sig[0,:],nf)
    sig2[1,:]=signal.resample(sig[1,:],nf)
    sig2[2,:]=signal.resample(sig[2,:],nf)
    return sig2