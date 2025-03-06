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
import device
name = 'WEST'
device.get_device_input(name)
from device import N0, N1, ap, R0, B0, xkap, del0, delt, T0e, T0i, T1, exponn, expont


# ZZZ0, RRR0, Rsep0,rho0, Nprf0, TprfE0,TprfI0
# %%
nt = 200
nr = 10
na = 200
theta = np.linspace(0, 2*np.pi-0.001, nt)
rho=np.linspace(0, ap, nr)
r = np.linspace(0, ap, na)
R = np.zeros((nr,nt))
Z = np.zeros((nr,nt))

for k in range (0,nr):      
    shaf = del0 *(1-(rho[k]/ap)**2)   
    
    for it in range (0,200):
        R[k,it]  = R0 + shaf + rho[k]*np.cos(theta[it])- delt*rho[k]**2*np.sin(theta[it])**2
        Z[k,it]  = xkap * rho[k] * np.sin(theta[it])
        
Rsep = R0 + ap*np.cos(theta+delt*ap*np.sin(theta))
Zsep = xkap * ap * np.sin(theta)
        
fig = plt.figure(1,figsize=(5, 8))


plt.plot(device.RRR0, device.ZZZ0, color = 'b', label = 'Ref',linewidth=3) 
plt.plot(device.Rsep0, device.Zsep0, 'x', color = 'b',linewidth=3) 
plt.plot(device.RRR1, device.ZZZ1, color = 'r', label = 'Ana-EVE') 
plt.plot(device.Rsep0, device.Zsep0, 'x', color = 'r') 

for k in range (0,nr-1):
    plt.plot(R[k,:],Z[k,:],'k--')

plt.plot(R[nr-1,:],Z[nr-1,:],'k', label = 'Ana-TOMCAT')
plt.scatter(R[0,0],0,color='red')

plt.xlabel(r'$\mathregular{R\ [m]}$') 
plt.ylabel(r'$\mathregular{Z\ [m]}$')
plt.legend()
plt.grid(True)
plt.tight_layout()

fig.gca().axis('equal')

#%%
def nt(A0, A1, expon):
    return  (A0-A1)*(1-(r/ap)**2)**expon + A1

Nprf = nt(N0,N1,exponn)
TprfE = nt(T0e,T1,expont)
TprfI = nt(T0i,T1,expont)

plt.figure(2,figsize=(10, 5))
plt.subplots_adjust(left=0.1,right=0.9,
                    wspace=0.3)

plt.subplot(1,2,1)
plt.plot(r,Nprf,'g',label='$n_e$')

YL = plt.ylim()
plt.ylim(0,YL[1])
plt.xlabel(r'$\rho$')
plt.ylabel(r'$\mathregular{n_e\ [m^{-3}]}$')
plt.legend()
plt.grid(True)
plt.tight_layout()

plt.subplot(1,2,2)
plt.plot(r,TprfE/1000,'g',label='$T_e$')
plt.plot(r,TprfI/1000,'g--',label='$T_i$')

YL = plt.ylim()
plt.ylim(0,YL[1])
plt.xlabel(r'$\rho$')
plt.ylabel(r'$\mathregular{T\ [eV]}$')
plt.legend()
plt.grid(True)
plt.tight_layout()
#%%
fig = plt.figure(figsize = (6.4, 2.92))
ax = fig.subplots(1, 1)
ax.plot(device.rho0, device.Nprf0, 'b', label = 'Ref',linewidth=3) 
ax.plot(device.rho1, device.Nprf1, 'r', label = 'Ana-EVE')
ax.set_xlabel(r'$s_p$')
ax.set_ylabel(r'$\mathregular{n_e\ [m^{-3}]}$')
ax.legend()
ax.grid(True)
plt.tight_layout()

fig = plt.figure(figsize = (6.4, 2.92))
ax = fig.subplots(1, 1)
ax.plot(device.rho0, device.TprfE0, 'b', label = 'Ref $T_e$',linewidth=3) 
ax.plot(device.rho1, device.TprfE1, 'r', label = 'Ana-EVE $T_e$')

ax.plot(device.rho0, device.TprfI0, 'b--', label = 'Ref $T_i$',linewidth=3) 
ax.plot(device.rho1, device.TprfI1, 'r--', label = 'Ana-EVE $T_i$')
ax.set_xlabel(r'$s_p$')
ax.set_ylabel(r'$\mathregular{T\ [eV]}$')
ax.legend()
ax.grid(True)
plt.tight_layout()
