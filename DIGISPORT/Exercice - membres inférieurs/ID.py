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
    #matrices de rotation du pelvis


    #matrice de rotation pelvis-> monde
    zp=(markers[21,:]-markers[10,:])/np.linalg.norm(markers[21,:]-markers[10,:]) # RFWT-LFWT
    xp=(0.5*(markers[21,:]+markers[10,:])-0.5*(markers[22,:]+markers[11,:]))/np.linalg.norm(0.5*(markers[21,:]+markers[10,:])-0.5*(markers[22,:]+markers[11,:]))
    yp=np.cross(zp,xp)
    xp=np.cross(yp,zp) # attention il faut bien s'assurer d'avoir une base orthonormée
    Rp=np.zeros((3,3))
    Rp[0,:]=[xp[0],yp[0],zp[0]]
    Rp[1,:]=[xp[1],yp[1],zp[1]]
    Rp[2,:]=[xp[2],yp[2],zp[2]]
    return Rp

def rotz(q):
    #matrices de rotation autour de z
    y=np.zeros((3,3))

    y[0,:]=[np.cos(q),-np.sin(q),0]
    y[1,:]=[np.sin(q),np.cos(q),0]
    y[2,:]=[0,0,1]
    return y

def LP_filter(x,f_coupure,f_echantillonnage):
    # filtre passe-bas, no phase shift
    b, a = signal.butter(4, f_coupure, btype='low', analog=False, output='ba', fs=f_echantillonnage)
    y = signal.filtfilt(b, a, x, padlen=150)
    return y

def LP_filter_vec3(x,f_coupure,f_echantillonnage):
    y=np.zeros((len(x[:,0]),3))
    y[:,0]=LP_filter(x[:,0],f_coupure,f_echantillonnage)
    y[:,1]=LP_filter(x[:,1],f_coupure,f_echantillonnage)
    y[:,2]=LP_filter(x[:,2],f_coupure,f_echantillonnage)
    
    return y

def diff2(x,dt):
    
    y=np.zeros(len(x))
    n=len(x)-1
    y[0]=(-x[2]+4*x[1]-3*x[0])/(2*dt) #différence avant d'ordre 2
    y[n]=(3*x[n]-4*x[n-1]+x[n-2])/(2*dt) #différence arrière d'ordre 2
    for i in range(1,n-1):
        y[i]=(x[i+1]-x[i-1])/(2*dt) #différence centrée d'ordre 2
    
    return y

def diff_vec3(x,dt):
    y=np.zeros((len(x[:,0]),3))
    y[:,0]=diff2(x[:,0],dt)
    y[:,1]=diff2(x[:,1],dt)
    y[:,2]=diff2(x[:,2],dt)
    
    return y

def generate_GRF(pf):
# pf_0['unit_force']          # Units of forces
# pf_0['unit_moment']         # Units of moments
       # Units of center of pressure

# pf_0['cal_matrix']          # Calibration matrix

# pf_0['origin']              # Position of the origin

# pf_0['force'][dim,frame]             # Force data
# pf_0['moment']              # Moment data
# pf_0['center_of_pressure']  # Center of pressure data
# pf_0['Tz']                  # Moment at center of pressure data

#     origin = 0.25*(pf['corners'][:,3]+pf['corners'][:,2]+pf['corners'][:,1]+pf['corners'][:,0])/1000
#     x = ((pf['corners'][:,3] + pf['corners'][:,0])/2) - ((pf['corners'][:,2] + pf['corners'][:,1])/2)
#     x = x/np.linalg.norm(x)
#     y = ((pf['corners'][:,0] + pf['corners'][:,1])/2) - ((pf['corners'][:,2] + pf['corners'][:,3])/2)
#     y = y/np.linalg.norm(y)
#     z = np.cross(x,y)
#     z = z/np.linalg.norm(z)
#     y = np.cross(z,x)
#     Rplatform=np.zeros((3,3))
#     Rplatform[0,:]=[x[0],y[0],z[0]]
#     Rplatform[1,:]=[x[1],y[1],z[1]]
#     Rplatform[2,:]=[x[2],y[2],z[2]]
        
    F=np.zeros((len(pf['force'][:,0]),len(pf['force'][0,:])))#effort
    P=np.zeros((len(pf['force'][:,0]),len(pf['force'][0,:])))#centre de pression
    for f in range(0,len(pf['force'][0,:])-1):
        #F[:,f]=Rplatform @ pf['force'][:,f]# effort dans la base mocap
        F[:,f]=pf['force'][:,f]# effort dans la base mocap
        #P[:,f]=origin+Rplatform @ pf['center_of_pressure'][:,f]/1000# effort dans la base mocap
        P[:,f]=pf['center_of_pressure'][:,f]/1000
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
    sig2=np.zeros((len(sig[:,0]),nf))
    sig2[0,:]=signal.resample(sig[0,:],nf)
    sig2[1,:]=signal.resample(sig[1,:],nf)
    sig2[2,:]=signal.resample(sig[2,:],nf)
    return sig2