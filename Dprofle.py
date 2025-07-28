#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 14:05:35 2025

@author: YJ281217
"""
# %%
## cleaning in SPYDER
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
name = 'CFEDR'
N_line = 1    #how many line you want to compare as reference
FWall  = 1   #first wall

device.get_device_input(name)
from device import N0, N1, ap, R0, B0, xkap, del0, delt, T0e, T0i, T1, exponn, expont, Zdel


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
    
    for it in range (0,nt):
        R[k,it]  = R0 + shaf + rho[k]*np.cos(theta[it])- delt*rho[k]**2*np.sin(theta[it])**2
        Z[k,it]  = xkap * rho[k] * np.sin(theta[it])
        
Rsep = R0 + ap*np.cos(theta+delt*ap*np.sin(theta))
Zsep = xkap * ap * np.sin(theta)
      
  
fig = plt.figure(figsize=(5, 8))

plt.plot(device.RRR0, device.ZZZ0, color = 'b', label = 'Ref',linewidth=3) 
plt.plot(device.Rsep0, device.Zsep0, 'x', color = 'b',linewidth=3) 

"""##second line##"""
if N_line == 2:
    plt.plot(device.RRR1, device.ZZZ1, color = 'r', label = 'Ana-EVE') 
    plt.plot(device.Rsep1, device.Zsep1, 'x', color = 'r') 
"""##second line##"""

for k in range (0,nr-1):
    plt.plot(R[k,:],Z[k,:]+Zdel,'g--')

plt.plot(R[nr-1,:],Z[nr-1,:]+Zdel,'g', label = 'Ana-TOMCAT')
plt.scatter(R[0,0],Zdel,color='green')

if FWall == 1:
    plt.plot(device.FWR, device.FWZ, color = 'k') 

plt.xlabel(r'$\mathregular{R\ [m]}$') 
plt.ylabel(r'$\mathregular{Z\ [m]}$')
plt.legend()
plt.grid(True)
plt.tight_layout()

fig.gca().axis('equal')

#%%
def build_profile(A0, A1, expon):
    return  (A0-A1)*(1-(r/ap)**2)**expon + A1

Nprf = build_profile(N0,N1,exponn)
TprfE = build_profile(T0e,T1,expont)
TprfI = build_profile(T0i,T1,expont)


#%%  plot density and temperature
fig = plt.figure(figsize = (6.4, 5.5))
ax = fig.add_subplot(2,1,1)


ax.plot(device.rho0, device.Nprf0, 'b', label = 'Ref',linewidth=3) 
"""##second line##"""
if N_line == 2:
    ax.plot(device.rho1, device.Nprf1, 'r--', label = 'Ana-EVE')
"""##second line##"""
ax.plot(r/ap,Nprf,'g:',label='$n_e$')

YL = ax.get_ylim()
ax.set_ylim(0,YL[1])
ax.set_xlim(0,1)
ax.set_ylabel(r'$\mathregular{n_e\ [m^{-3}]}$')
ax.legend()
ax.grid(True)
plt.tight_layout()

#fig = plt.figure(figsize = (6.4, 2.92))
ax = fig.add_subplot(2, 1, 2)


ax.plot(device.rho0, device.TprfE0/1000, 'b', label = 'Ref $T_e$',linewidth=3)
ax.plot(device.rho0, device.TprfI0/1000, 'b--', label = 'Ref $T_i$',linewidth=3) 
"""##second line##""" 
if N_line == 2:
    ax.plot(device.rho1, device.TprfE1/1000, 'r', label = 'Ana-EVE $T_e$')
    ax.plot(device.rho1, device.TprfI1/1000, 'r--', label = 'Ana-EVE $T_i$')
"""##second line##""" 
ax.plot(r/ap,TprfE/1000,'g',label='$T_e$')
ax.plot(r/ap,TprfI/1000,'g--',label='$T_i$')

YL = ax.get_ylim()
ax.set_ylim(0,YL[1])
ax.set_xlim(0,1)
ax.set_xlabel(r'$\rho$')
ax.set_ylabel(r'$\mathregular{T\ [keV]}$')
ax.legend()
ax.grid(True)
plt.tight_layout()

#%%
R_ml = np.zeros(na)
R_mh = np.zeros(na)
for k in range (0,na):      
    shaf = del0 *(1-(r[k]/ap)**2)   
    
    R_ml[k]  = R0 + shaf + r[k]*1
    R_mh[k]  = R0 + shaf + r[k]*(-1)

fig = plt.figure(figsize=(6.4, 5.5))
ax = fig.add_subplot(2, 1, 1)

ax.plot(R_ml, Nprf, 'g', label=r'$n_e$ left')
ax.plot(R_mh, Nprf, 'g', label=r'$n_e$ right')
ax.set_ylabel(r'$n_e\ [m^{-3}]$')
ax.set_title(r'$n_e(R)$ at Midplane ($Z = 0$)')
ax.grid(True)
ax.legend()

ax = fig.add_subplot(2, 1, 2)
ax.plot(R_ml, TprfE/1000, 'g', label=r'$T_e$')
ax.plot(R_mh, TprfE/1000, 'g')
ax.plot(R_ml, TprfI/1000, 'g--', label=r'$T_i$')
ax.plot(R_mh, TprfI/1000, 'g--')

ax.set_xlabel('R [m]')
ax.set_ylabel(r'$\mathregular{T\ [keV]}$')
ax.set_title(r'$T(R)$ at Midplane ($Z = 0$)')
ax.grid(True)
ax.legend()

fig.tight_layout()
