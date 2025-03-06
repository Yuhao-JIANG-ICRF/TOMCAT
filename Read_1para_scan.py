#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  6 09:48:23 2025

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
from device import get_model_data
import gc as gc
import os
gc.collect()    #release ram
plt.close('all')     #close all figures

# %%
'Basic parameters'

pre = 'output/1para_scan'
UNDERSCORE='_'       #For MAC

if os.path.exists(pre):
    subfolders = [item for item in os.listdir(pre)
                  if os.path.isdir(os.path.join(pre, item))]
    
Ncase = len(subfolders) 

con   = []
B0    = []
freq  = []
Ntor  = []
P_sum = []
   
for i in range(Ncase):
    pre0 = 'output/1para_scan'+f'/{subfolders[i]}/'
    
    filename4 = pre0 + 'fort.31' #Poynting Flux and total power 
    aux31 = np.loadtxt(filename4)
    filename5 = pre0 + 'fort.30' #Power per species 
    aux30 = np.loadtxt(filename5)
    filename1 = pre0 + 'fort.80' + UNDERSCORE  #Basic parameters
    aux80 = np.loadtxt(filename1)
    
    name, Ns, name_species, conc, a, R0, Rmag = get_model_data(pre0)
    P_sumn = np.zeros((Ns+1))  
        

    
    B0n   = aux80[3]   #Magnetic field [T]
    freqn = aux80[2]/1e6   #ICRF frequency [MHz]
    Ntorn = aux80[1]   #

    xP      = aux31[:,0];
    Ptot    = aux31[:,1];
    Flux_re = aux31[:,3];
    Flux_im = aux31[:,4];

    xPow  = aux30[:,0];
    Pabs  = aux30[:,1:2*Ns:2];   #Ns species

    xd = np.diff(xPow)                   # computes differences; length is Nx-1
    xd = np.concatenate((xd, [xd[-1]]))  # append last element (equal to the previous difference)
    Pabs_sum = np.empty_like(Pabs)
    for k in range(Ns):
        Pabs_sum[:, k] = np.cumsum(Pabs[:, k] * xd)*100
        P_sumn[k+1] = Pabs_sum[-1, k]
    Ptot_sum = np.sum(Pabs_sum, axis=1, keepdims=True)
    P_sumn[0] = Ptot_sum[-1,0]
    
    con.append(conc)
    B0.append(B0n)
    freq.append(freqn)
    Ntor.append(Ntorn)
    P_sum.append(P_sumn)
P_sum = np.array(P_sum)
#%%
lists = {
    'B_0'             : (B0,   '[T]'),
    'f_{ICRF}'       : (freq, '[MHz]'),
    'N_{tor}'         : (Ntor, ''),
    name_species[-1]  : (con,  '[%]') 
    }
ftitle = f'{name}:'
ftitle1 = ''
ftitle2 = ''
for key, lst in lists.items():
    arr = np.array(lst[0])
    if np.all(arr == arr.flat[0]):
        if key != name_species[-1]:
            ftitle = ftitle+f'${key}$ = {lst[0][0]} {lst[1]}, '
        else:
            for i in range(Ns-1):
                if i==0:
                    ftitle1 = ftitle1 + '${%s} $' % name_species[i+1]
                    ftitle2 = ftitle2 + f'{conc[i]} %'
                else:
                    ftitle1 = ftitle1 + ': ${%s}$' % name_species[i+1]
                    ftitle2 = ftitle2 + f': {conc[i]}%'
            ftitle = ftitle + '\n' + ftitle1 +' = '+ftitle2 
    else:
        xlabel = 'X[${%s}$]' % key + f'{lst[1]}'
        xx = arr[:,-1]
 
    
'PLOT FIGURE 4----power deposition'
plt.rcParams.update({
    'lines.linewidth':2,
    'xtick.direction':'in',
    'ytick.direction':'in',
    'xtick.top':True,
    'xtick.bottom':True,
    'ytick.left':True,
    'ytick.right':True,
    'axes.grid':True
    })

colors = ['r', 'b', 'g', 'purple']


plt.figure(figsize=(6, 5))
plt.suptitle(ftitle)  
   
'Power absorption'
plt.plot(xx,P_sum[:,0],'k--',label='Total')
for i in range(Ns):  
        plt.plot(xx,P_sum[:,i+1],color=colors[i],label=f'${name_species[i]}$')
YL = plt.ylim()
plt.ylim(0,YL[1])
plt.xlim(np.min(xx),np.max(xx))
plt.ylabel('Absorbed Power [W/m]')
plt.xlabel(xlabel)
plt.legend()


Mtot = np.argmax(P_sum[:,0])
Mm = np.argmax(P_sum[:,-1])
Sumi = np.sum(P_sum[:,2:],axis=1)
Mi = np.argmax(Sumi)


print(f'The maximum of the total power absorption is Case {Mtot}:')
print(f'{xx[Mtot]}, $P_{{tot}}$ = {P_sum[Mtot,0]:.1f}%')
print(f'The maximum of the ions power absorption is Case {Mi}:')
print(f'{xx[Mi]}, $P_{{ions}}$ = {Sumi[Mi]:.1f}%')