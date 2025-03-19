#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  6 09:48:23 2025

@author: YJ281217
"""

# %%
## cleaning in SPYDER
get_ipython().magic('reset -f')
get_ipython().magic('clear')
##

# %%
globals().clear()
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

pre = 'output/2para_scan'
pref = 'output/'     #Save figures
UNDERSCORE='_'       #For MAC

con_scan = 0

if os.path.exists(pre):
    subfolders1 = [item for item in os.listdir(pre)
                  if os.path.isdir(os.path.join(pre, item))]
if os.path.exists(pre+'/1'):
    subfolders2 = [item for item in os.listdir(pre+'/1')
                  if os.path.isdir(os.path.join(pre+'/1', item))]
Ncase1 = len(subfolders1) 
Ncase2 = len(subfolders2)

con   = []
B0    = []
freq  = []
Ntor  = []
P_sum = []

for i in range(Ncase1):
    prem = pre+f'/{subfolders1[i]}/'
    for k in range(Ncase2):
        pre0 = prem+f'/{subfolders2[k]}/'
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
        
def cod(para):
    para1 = np.array(para)
    para2 = np.atleast_2d(para1)   
    if para2.shape[0]==1:
        para2 = para2.T
        paraf = para2[:,-1].reshape(Ncase1,Ncase2)
        return paraf   
    else:
        paraf1 = para2[:,-2].reshape(Ncase1,Ncase2)
        paraf2 = para2[:,-1].reshape(Ncase1,Ncase2)
        return paraf1,paraf2
        
   
    
        


confy, confx = cod(con)
B0f = cod(B0)
freqf = cod(freq)
Ntorf = cod(Ntor)

Psf = np.empty((Ns+1,Ncase1,Ncase2))        
P_sum = np.array(P_sum)       
for i in range(Ns+1): 
    Ps = P_sum[:,i].reshape(Ncase1,Ncase2)
    Psf[i] = Ps

#%%
def labelf(B0,freq,Ntor,conx,cony=None):
    lists = {
        'B_0'             : (B0,   '[T]', 'B0'),
        'f_{ICRF}'       : (freq, '[MHz]','freq'),
        'N_{tor}'         : (Ntor, '','Ntor'),
        f'X[{name_species[-1]}]'  : (conx,  '[%]','Mi-concentr')
        }
    if cony is not None:
        lists[f'X[{name_species[-2]}]']=(cony,  '[%]','Ma-concentr')
    return lists
def labeln(lists):    
    for key, lst in lists.items():
        arr = np.array(lst[0])
        if np.all(arr == arr[0]) != True:
            label = '$%s$' % key + lst[1]
            delk=key
            xx = arr
            minarr = np.min(arr)
            maxarr = np.max(arr)
            figname = lst[2]+f'_scan:{minarr}-{maxarr}'
    return label, xx, delk ,figname      
    
if con_scan==0:
    lists=labelf(B0f[:,0],freqf[:,0],Ntorf[:,0],confx[:,0])
    ylabel,yaxis,del1,figname1 = labeln(lists)
    lists=labelf(B0f[0,:],freqf[0,:],Ntorf[0,:],confx[0,:])
    xlabel,xaxis,del2,figname2 = labeln(lists)
else:
    lists=labelf(B0f[:,0],freqf[:,0],Ntorf[:,0],confx[:,0],confy[:,0])
    ylabel,yaxis,del1,figname1 = labeln(lists)
    lists=labelf(B0f[0,:],freqf[0,:],Ntorf[0,:],confx[0,:],confy[:,0])
    xlabel,xaxis,del2,figname2 = labeln(lists)
figname= figname1+'&'+figname2    
del lists[del1]
del lists[del2]
            
ftitle = f'{name}:'
ftitle1 = ''
ftitle2 = ''
for key, lst in lists.items():
    if key != 'X[{name_species[-1]}]' :
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


plt.figure(figsize=(6, 5))
plt.suptitle(ftitle+'\n Total Power absorption [%]')  
levels=np.linspace(0,100,50)
cp1 = plt.contourf(xaxis,yaxis,Psf[0],levels=levels)
plt.colorbar(cp1)
plt.ylabel(ylabel)
plt.xlabel(xlabel)


for i in range(Ns):  
    plt.figure(figsize=(6, 5))
    plt.suptitle(ftitle+'\n'+'${%s} $' % name_species[i] + ' Power absorption [%]')  
    levels=np.linspace(0,100,50)
    cp1 = plt.contourf(xaxis,yaxis,Psf[i+1],levels=levels)
    plt.colorbar(cp1)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)


Mtot = np.unravel_index(np.argmax(Psf[0]),Psf[0].shape)
nMtot = tuple(x+1 for x in Mtot)
Sumi = np.sum(Psf[2:],axis=0)
Mi = np.unravel_index(np.argmax(Sumi),Sumi.shape)
nMi = tuple(x+1 for x in Mi)



print(f'The maximum of the total power absorption is Case {nMtot}:')
print(f'{xaxis[Mtot[1]]},{yaxis[Mtot[0]]}, $P_{{tot}}$ = {Psf[0,Mtot[0],Mtot[1]]:.1f}%')
print(f'The maximum of the ions power absorption is Case {nMi}:')
print(f'{xaxis[Mi[1]]},{yaxis[Mi[0]]}, $P_{{ions}}$ = {Sumi[Mi[0],Mi[1]]:.1f}%')

# %% saving figures
import os

folder = pref+name+'/2para_scan/'
for i in range(Ns-2):
    folder = folder + name_species[i+1]
    
folder = folder+'('+name_species[-1]+')/'
print(folder)
    
if not os.path.exists(folder):
    os.makedirs(folder)   # create the folder if it does't exis


for i in plt.get_fignums():
    plt.figure(i)
    if i==1:       
        filename = f'Total-ABS-{figname}.png'
    else:
        filename = f'{name_species[i-2]}-ABS-{figname}.png'
    filepath = os.path.join(folder, filename)
    plt.savefig(filepath)


   