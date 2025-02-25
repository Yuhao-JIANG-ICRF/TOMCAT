#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 10:33:13 2025

@author: YJ281217
"""
# %%
globals().clear()
get_ipython().magic('reset -f')
get_ipython().magic('clear')

# %%
'import Library or Package and release ram & close fig' 
import numpy as np
import matplotlib.pyplot as plt
from device import get_model_data
import gc as gc
gc.collect()    #release ram
plt.close('all')     #close all figures

# %%
'Basic parameters'

pre = 'output/'
UNDERSCORE='_'       #For MAC

'''
name    = name of device,
Ns      = number of species,
name_species = name of species (including electrons),
conc    = concentrations (only ions),
a       = minor radius, 
R0      = major radius, 
Rmag    = position magnectic axis
'''
name, Ns, name_species, conc, a, R0, Rmag = get_model_data(pre)

'''
if iuseR = 1 , the central is at 0,
if iuseR = 1 , the central is at R0,
if isueR = 2 , the central is at Rmag
'''
iuseR=0    

print("Device:", name)
print("Number of species:", Ns)
print("Species:", name_species)
print("Concentrations:", conc,'%')

# %%
'READ DATA'
filename1 = pre + 'fort.80' + UNDERSCORE  #Basic parameters
aux80 = np.loadtxt(filename1)
filename2 = pre + 'fort.25' #Density and temp
aux25 = np.loadtxt(filename2)
filename3 = pre + 'fort.29' #Electric field
aux29 = np.loadtxt(filename3)
filename4 = pre + 'fort.31' #Poynting Flux and total power 
aux31 = np.loadtxt(filename4)
filename5 = pre + 'fort.30' #Power per species 
aux30 = np.loadtxt(filename5)
filename6 = pre + 'fort.72' + UNDERSCORE #Power fraction
aux72 = np.loadtxt(filename6)
filename7 = pre + 'fort.27' #Some extra quantities
aux27 = np.loadtxt(filename7)


B0   = aux80[3]   #Magnetic field [T]
freq = aux80[2]/1e6   #ICRF frequency [MHz]
Ntor = aux80[1]   #
conc = aux80[0]   #concentration %

xn = aux25[:,0];
n  = aux25[:,3:3+Ns]/1e19;       #Density for n species
T  = aux25[:,3+Ns:3+2*Ns]/1e3;   #Temp for n species 

xE        = aux29[:,0];
Eplus_re  = aux29[:,1];
Eplus_im  = aux29[:,2];
Eminus_re = aux29[:,3];
Eminus_im = aux29[:,4];
Epar_re   = aux29[:,5];
Epar_im   = aux29[:,6];

xP      = aux31[:,0];
Ptot    = aux31[:,1];
Flux_re = aux31[:,3];
Flux_im = aux31[:,4];

xPow  = aux30[:,0];
Pabs  = aux30[:,1:2*Ns:2];   #Ns species

xd = np.diff(xPow)                 # computes differences; length is Nx-1
xd = np.concatenate((xd, [xd[-1]]))  # append last element (equal to the previous difference)
Pabs_sum = np.empty_like(Pabs)
for i in range(Ns):
    Pabs_sum[:, i] = np.cumsum(Pabs[:, i] * xd)*100
Ptot_sum = np.sum(Pabs_sum, axis=1, keepdims=True)




powerfrac  = aux72[1:Ns+1];    #Ns species
powround   = np.round(100*powerfrac);
xf    = aux27[:,0];
disp_roots = aux27[:,1:8:2];
 
B1 = np.min(xf);      #Bonundary at HFS
B2 = np.max(xf);      #Bonundary at LFS
DoublePass=Flux_re[-1]-Flux_re[0];


# %%

'PLOT FIGURE'
colors = ['r', 'b', 'g', 'k']

plt.figure()
plt.subplots_adjust(left=0.1,right=0.9,
                    wspace=0.3)
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
plt.suptitle(f'''
             {name}:
             $B_0$={B0:.1f} [T],
             $f_{{ICRF}}$={freq:.0f} [MHz],
             $N_{{tor}}$={Ntor:.0f}
             '''.replace('\n',''))
'Density'
plt.subplot(2,3,1)
for i in range(Ns):
    if i==1:
        plt.plot(xn,n[:,i],color=colors[i])
    elif n[0,i]==n[0,i-1]:
        continue
    else:
        plt.plot(xn,n[:,i],color=colors[i])        
YL = plt.ylim()
XL = np.max(xn)
plt.plot([XL,XL],[0,YL[1]],'k--')

plt.xlim(0,a)
plt.ylim(0,YL[1])

plt.title('')
plt.xlabel('r')
plt.ylabel('Density [1e19$m^{-3}$]')

'Temp'
plt.subplot(2,3,4)
for i in range(Ns):
    if i==1:    
        plt.plot(xn,T[:,i],color=colors[i])
    elif T[0,i]==T[0,i-1]:
        continue
    else: 
        plt.plot(xn,T[:,i],color=colors[i])
    
YL = plt.ylim()
XL = np.max(xn)
plt.plot([XL,XL],[0,YL[1]],'k--')

plt.xlim(0,a)
plt.ylim(0,YL[1])

plt.legend(['$T_{e}$','$T_{i}$'])
plt.title('')
plt.xlabel('r')
plt.ylabel('Temp [keV]')

'Disp root'
plt.subplot(2,3,2)
plt.plot(xP,disp_roots)
YL = plt.ylim()
XL = np.max(xn)
#plt.plot([XL,XL],[0,YL[1]],'k--')

plt.xlim(-a,a)
plt.ylim(YL[0],YL[1])

plt.title('')
plt.xlabel('r')
plt.ylabel('Disp roots')

'Poynting Flux'
plt.subplot(2,3,3)
plt.plot(xP,Flux_re,'b')
plt.plot(xP,Flux_im,'r')
YL = plt.ylim()
XL = np.max(xn)
#plt.plot([XL,XL],[0,YL[1]],'k--')

plt.xlim(-a,a)
plt.ylim(YL[0],YL[1])

plt.legend(['Real','Image'])
plt.title('')
plt.xlabel('r')
plt.ylabel('Poynting Flux')

'Electric field'
plt.subplot(2,3,5)
plt.plot(xP,Eplus_re,'b')
plt.plot(xP,Eplus_im,'r')
plt.plot(xP,np.sqrt(Eplus_re**2+Eplus_im**2),'k--')
YL = plt.ylim()
XL = np.max(xn)
#plt.plot([XL,XL],[0,YL[1]],'k--')

plt.xlim(-a,a)
plt.ylim(YL[0],YL[1])

plt.legend(['Real','Image'])
plt.title('')
plt.xlabel('r')
plt.ylabel('|$E^{+}$|')

'Power absorption'
plt.subplot(2,3,6)
plt.plot(xP,Pabs)
plt.plot(xP,Ptot,'k--')
YL = plt.ylim()
XL = np.max(xn)
#plt.plot([XL,XL],[0,YL[1]],'k--')

plt.xlim(-a,a)
plt.ylim(0,YL[1])

plt.legend(['$P_{tot}$'])
plt.title('')
plt.xlabel('r')
plt.ylabel('Absorbed Power [W/m]')

plt.show

# %%
'PLOT FIGURE'
colors = ['r', 'b', 'g', 'k']

plt.figure()
plt.subplots_adjust(left=0.1,right=0.9,
                    wspace=0.3)
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
plt.suptitle(f'''
             {name}:
             $B_0$={B0:.1f} [T],
             $f_{{ICRF}}$={freq:.0f} [MHz],
             $N_{{tor}}$={Ntor:.0f}
             '''.replace('\n',''))
'Density'
plt.subplot(1,2,1)
for i in range(Ns):
    if i==0:
        plt.plot(xn,n[:,i],color=colors[i],label=f'n$_{{{name_species[i]}}}$')
    elif n[0,i]==n[0,i-1]:
        continue
    else:
        plt.plot(xn,n[:,i],color=colors[i],label=f'n$_i$')  
        

YL = plt.ylim()
XL = np.max(xn)
plt.plot([XL,XL],[0,YL[1]],'k--')

plt.xlim(0,a)
plt.ylim(0,YL[1])

plt.legend()
plt.title('')
plt.xlabel('r')
plt.ylabel('Density [1e19$m^{-3}$]')

'Temp'
plt.subplot(1,2,2)
for i in range(Ns):
    if i==0:    
        plt.plot(xn,T[:,i],color=colors[i])
    elif T[0,i]==T[0,i-1]:
        continue
    else: 
        plt.plot(xn,T[:,i],color=colors[i])
    

YL = plt.ylim()
XL = np.max(xn)
plt.plot([XL,XL],[0,YL[1]],'k--')

plt.xlim(0,a)
plt.ylim(0,YL[1])

plt.legend(['$T_{e}$','$T_{i}$'])
plt.title('')
plt.xlabel('r')
plt.ylabel('Temp [keV]')


plt.show
# %%
'PLOT FIGURE----Poyting flux and Electric field'
colors = ['r', 'b', 'g', 'k']

plt.figure()
plt.subplots_adjust(left=0.1,right=0.9,
                    wspace=0.3)
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
plt.suptitle(f'''
             {name}:
             $B_0$={B0:.1f} [T],
             $f_{{ICRF}}$={freq:.0f} [MHz],
             $N_{{tor}}$={Ntor:.0f}
             '''.replace('\n',''))

'Poynting Flux'
plt.subplot(1,2,1)
plt.plot(xP,Flux_re,'b')
plt.plot(xP,Flux_im,'r')
YL = plt.ylim()
XL = np.max(xn)
#plt.plot([XL,XL],[0,YL[1]],'k--')

plt.xlim(-a,a)
plt.ylim(YL[0],YL[1])

plt.legend(['Real','Image'])
plt.title('')
plt.xlabel('r')
plt.ylabel('Poynting Flux')

'Electric field'
plt.subplot(1,2,2)
plt.plot(xP,Eplus_re,'b')
plt.plot(xP,Eplus_im,'r')
plt.plot(xP,np.sqrt(Eplus_re**2+Eplus_im**2),'k--')
YL = plt.ylim()
XL = np.max(xn)
#plt.plot([XL,XL],[0,YL[1]],'k--')

plt.xlim(-a,a)
plt.ylim(YL[0],YL[1])

plt.legend(['Real','Image'])
plt.title('')
plt.xlabel('r')
plt.ylabel('|$E^{+}$|')

plt.show


# %%
'PLOT FIGURE'
colors = ['r', 'b', 'g', 'k']

plt.figure()
plt.subplots_adjust(left=0.1,right=0.9,
                    wspace=0.3)
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
plt.suptitle(f'''
             {name}:
             $B_0$={B0:.1f} [T],
             $f_{{ICRF}}$={freq:.0f} [MHz],
             $N_{{tor}}$={Ntor:.0f}
             '''.replace('\n',''))

'Power absorption'
plt.plot(xP,Ptot,'k--',label=f'$P_{{tot}}$:{Ptot_sum[-1,0]:.0f}%')
for i in range(Ns):  
        plt.plot(xP,Pabs[:,i],color=colors[i],label=f'{name_species[i]}:{Pabs_sum[-1,i]:.0f}%')

YL = plt.ylim()
XL = np.max(xn)
#plt.plot([XL,XL],[0,YL[1]],'k--')

plt.xlim(-a,a)
plt.ylim(0,YL[1])

plt.legend()
plt.title('')
plt.xlabel('r')
plt.ylabel('Absorbed Power [W/m]')

plt.show




