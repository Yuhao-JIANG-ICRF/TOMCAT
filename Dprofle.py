#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 14:05:35 2025

@author: YJ281217
"""
# %%
## cleaning in TOMCAT
get_ipython().magic('reset -f')
get_ipython().magic('clear')
##
globals().clear()
# %%
'import Library or Package and release ram & close fig' 
import numpy as np
import matplotlib.pyplot as plt
import gc as gc

gc.collect()    #release ram
plt.close('all')     #close all figures

# %%


namecase='WEST'
N0=6.02e19
N1=2.25e19
ap=0.43
R0=2.5
B0=3.657
xkap=1.4
del0=0.05
delt=0.8
T0e=1.609e3
T0i=1.37e3
T1=0.25e3
exponn=1.0
expont=1.5

theta = np.linspace(0, 2*np.pi-0.001, 200)
rho=np.linspace(0, ap, 20)
R = np.zeros((20,200))
Z = np.zeros((20,200))

for k in range (0,20):   
    
    shaf = del0 *(1-(rho[k]/ap)**2)   
    
    for it in range (0,200):
        R[k,it]  = R0 + shaf + rho[k]*np.cos(theta[it])- delt*rho[k]**2*np.sin(theta[it])**2
        Z[k,it]  = xkap * rho[k] * np.sin(theta[it])
        
Rsep = R0 + ap*np.cos(theta+delt*ap*np.sin(theta))
Zsep = xkap * ap * np.sin(theta)
        
fig = plt.figure(1,figsize=(5, 8))
for k in range (0,20):
    plt.plot(R[k,:],Z[k,:],'--')

plt.plot(Rsep,Zsep,'k:')
plt.scatter(R0,0)

#fig = plt.figure(1)
#fig.axis('equal')
fig.gca().axis('equal')

r = np.linspace(0, ap, 200)
def nt(A0, A1, expon):
    return  (A0-A1)*(1-(r/ap)**2)**expon + A1

Nprf = nt(N0,N1,exponn)
TprfE = nt(T0e,T1,expont)
TprfI = nt(T0i,T1,expont)

plt.figure(2,figsize=(10, 5))
plt.subplots_adjust(left=0.1,right=0.9,
                    wspace=0.3)

plt.subplot(1,2,1)
plt.plot(r,Nprf,'b')

plt.subplot(1,2,2)
plt.plot(r,TprfE,'b',label='$T_e$')
plt.plot(r,TprfI,'r',label='$T_i$')




