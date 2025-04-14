#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 10:33:13 2025

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
from device import get_model_data
import gc as gc
gc.collect()    #release ram
plt.close('all')     #close all figures

# %%
'Basic parameters'

pre = 'output/'
#pre = 'output/1para_scan/22/'
#pre = 'output/2para_scan/80/10/'

pref = 'output/'     #Save figures
UNDERSCORE='_'       #For MAC
'''
if iuseR = 0 , the central is at 0,
if iuseR = 1 , the central is at R0
'''
iuseR=0     #switch to choose the x axis
Bon = 0     #swithc to have the boundry, 0--off; 1--on
MA = 0     #swithc to plot magnetic axis, 0--off; 1--on
RA = 0     #swithc to plot major radius, 0--off; 1--on

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


print("Device:", name)
print("Number of species:", Ns)
print("Species:", name_species)
print("Concentrations:", conc,'%')
if MA==1:
    print('Magnetic Axis: ON')
if RA==1:
    print('Major Radius: ON')
if Bon==1:
    print('Bondary: ON')

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
filename7 = pre + 'fort.27' #
aux27 = np.loadtxt(filename7)


B0   = aux80[3]   #Magnetic field [T]
freq = aux80[2]/1e6   #ICRF frequency [MHz]
Ntor = aux80[1]   #

xn = aux25[:,0]/a;
n  = aux25[:,3:3+Ns];       #Density for n species
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

xd = np.diff(xPow)                   # computes differences; length is Nx-1
xd = np.concatenate((xd, [xd[-1]]))  # append last element (equal to the previous difference)
Pabs_sum = np.empty_like(Pabs)
for i in range(Ns):
    Pabs_sum[:, i] = np.cumsum(Pabs[:, i] * xd)*100
Ptot_sum = np.sum(Pabs_sum, axis=1, keepdims=True)




powerfrac  = aux72[1:Ns+1];    #Ns species
powround   = np.round(100*powerfrac);
xf    = aux27[:,0];
disp_roots = aux27[:,1];
 

DoublePass=Flux_re[-1]-Flux_re[0];
# %%



# %%
'PLOT FIGURE'
colors = ['r', 'b', 'g', 'purple']
ftitle = f'''{name}: 
$B_0$={B0:.1f} [T], 
$f_{{ICRF}}$={freq:.1f} [MHz], 
$N_{{tor}}$={Ntor:.0f}
'''.replace('\n','')
ftitle1 = ''
ftitle2 = ''
for i in range(Ns-1):
    if i==0:
        ftitle1 = ftitle1 + '${%s} $' % name_species[i+1]
        ftitle2 = ftitle2 + f'{conc[i]} %'
    else:
        ftitle1 = ftitle1 + ': ${%s}$' % name_species[i+1]
        ftitle2 = ftitle2 + f': {conc[i]}%'
    
ftitlef = ftitle + '\n' + ftitle1 +' = '+ftitle2

if iuseR == 1:
    R = R0
    xlab = 'R [m]'
if iuseR == 0:
    R = 0
    xlab = 'x [m]'
    

XL = np.max(xn)
xE = xE+R
xf = xf+R
xP = xP+R
xPow = xPow+R
xmin = R-a
xmax = R+a
B1 = np.min(xf)      #Bonundary at HFS
B2 = np.max(xf)      #Bonundary at LFS

plt.rcParams.update({
    'lines.linewidth':2,
    'xtick.direction':'in',
    'ytick.direction':'in',
    'xtick.top':True,
    'xtick.bottom':True,
    'ytick.left':True,
    'ytick.right':True,
    'axes.grid':True,
    'axes.labelsize':16,
    'xtick.labelsize':14,
    'ytick.labelsize':14,
    'legend.fontsize':13,
    'figure.titlesize':15
    })
# %%
plt.figure(1,figsize=(11, 7))
plt.subplots_adjust(left=0.1,right=0.9,
                    wspace=0.3)

plt.suptitle(ftitlef)

'Density'
plt.subplot(2,3,1)
for i in range(Ns):
    if i==0:
        plt.plot(xn,n[:,i],color=colors[i],label='$n_{%s}$' % name_species[i])
    elif n[0,i]==n[0,i-1]:
        continue
    else:
        plt.plot(xn,n[:,i],color=colors[i],label='$n_{%s}$' % name_species[i])    
        
YL = plt.ylim()
plt.xlim(0,1)
plt.ylim(0,YL[1])
plt.title('')
plt.ylabel('Density [$m^{-3}$]')

'Temp'
plt.subplot(2,3,4)
for i in range(Ns):
    if i==0:    
        plt.plot(xn,T[:,i],color=colors[i],label='$T_{e}$')
    elif T[0,i]==T[0,i-1]:
        continue
    else: 
        plt.plot(xn,T[:,i],color=colors[i],label='$T_{i}$')
    
YL = plt.ylim()
plt.xlim(0,1)
plt.ylim(0,YL[1])
plt.title('')
plt.xlabel('$\\rho$')
plt.ylabel('Temp [keV]')

'Disp root'
plt.subplot(2,3,2)
plt.plot(xP,disp_roots,label='Fast Wave')
YL = plt.ylim()
XL = np.max(xn)
#plt.plot([XL,XL],[0,YL[1]],'k--')

plt.xlim(xmin,xmax)
plt.ylim(YL[0],YL[1])

plt.title('')
plt.ylabel('Disp roots')

'Poynting Flux'
plt.subplot(2,3,3)
plt.plot(xP,Flux_re,'b',label='Real')
plt.plot(xP,Flux_im,'r',label='Image')

YL = plt.ylim()
plt.ylim(YL[0],YL[1])
plt.xlim(xmin,xmax)
plt.title('')
plt.ylabel('Poynting Flux')

'Electric field'
plt.subplot(2,3,5)
plt.plot(xP,np.sqrt(Eplus_re**2+Eplus_im**2),'g--',label='|E|')
plt.plot(xP,Eplus_re,'b',label='Real')
plt.plot(xP,Eplus_im,'r',label='Image')

YL = plt.ylim()
plt.xlim(xmin,xmax)
plt.ylim(YL[0],YL[1])
plt.title('')
plt.xlabel(xlab)
plt.ylabel('|$E^{+}$|')

'Power absorption'
plt.subplot(2,3,6)
plt.plot(xP,Ptot,'k--',label=f'$P_{{total}}$:{Ptot_sum[-1,0]:.0f}%')
for i in range(Ns):  
        plt.plot(xP,Pabs[:,i],color=colors[i],label=f'${name_species[i]}$:{Pabs_sum[-1,i]:.0f}%')

YL = plt.ylim()
plt.ylim(0,YL[1])
plt.xlabel(xlab)
plt.title('')
plt.ylabel('Absorbed Power [W/m]')

fig = plt.figure(1)
for ax in fig.axes:
    ax.legend()
# %%
'PLOT FIGURE 2'

plt.figure(2,figsize=(10, 5))
plt.subplots_adjust(left=0.1,right=0.9,
                    wspace=0.3)

'Density'
plt.subplot(1,2,1)
for i in range(Ns):
    if i==0:
        plt.plot(xn,n[:,i],color=colors[i],label='$n_{%s}$' % name_species[i])
    elif n[0,i]==n[0,i-1]:
        continue
    else:
        plt.plot(xn,n[:,i],color=colors[i],label='$n_{%s}$' % name_species[i])  
        

plt.title('')
plt.ylabel('Density [$m^{-3}$]')

'Temp'
plt.subplot(1,2,2)
for i in range(Ns):
    if i==0:    
        plt.plot(xn,T[:,i],color=colors[i],label='$T_{e}$')
    elif T[0,i]==T[0,i-1]:
        continue
    else: 
        plt.plot(xn,T[:,i],color=colors[i],label='$T_{i}$')
    

plt.title('')
plt.ylabel('Temp [keV]')

# %%
'PLOT FIGURE----Poyting flux and Electric field'

plt.figure(3,figsize=(10, 5))
plt.subplots_adjust(left=0.1,right=0.9,
                    wspace=0.3)

'Poynting Flux'
plt.subplot(1,2,1)

plt.plot(xP,Flux_re,'b',label='Real')
plt.plot(xP,Flux_im,'r',label='Image')

plt.title('')
plt.ylabel('Poynting Flux')

'Electric field'
plt.subplot(1,2,2)
plt.plot(xP,np.sqrt(Eplus_re**2+Eplus_im**2),'g--',label='|E|')
plt.plot(xP,Eplus_re,'b',label='Real')
plt.plot(xP,Eplus_im,'r',label='Image')


plt.title('')
plt.ylabel('|$E^{+}$|')

# %%
'PLOT FIGURE 4----power deposition'

plt.figure(4,figsize=(6, 5.5))
plt.subplots_adjust(left=0.1,right=0.9,
                    wspace=0.3)

'Power absorption'
plt.plot(xP,Ptot,'k--',label=f'$Total$:{Ptot_sum[-1,0]:.0f}%')
for i in range(Ns):  
        plt.plot(xP,Pabs[:,i],color=colors[i],label=f'${name_species[i]}$:{Pabs_sum[-1,i]:.0f}%')

plt.title('')
plt.ylabel('Absorbed Power [W/m]')

#%%
for i in plt.get_fignums():
    fig = plt.figure(i)
    plt.suptitle(ftitlef)
    if i==1:
        continue
    for ax in fig.axes:
        YL = ax.get_ylim()
        if i==3:
            ax.set_ylim(YL[0],YL[1])
        else:
            ax.set_ylim(0,YL[1])
        if i==2:
            ax.set_xlim(0,1)
            ax.set_xlabel('$\\rho$')
        else:
            ax.set_xlim(xmin,xmax)  
            ax.set_xlabel(xlab)
        
            if MA==1:
                ax.axvline(x=xmax-a+Rmag-R0, color='k',linestyle='--',label='magnetic axis',linewidth=1.5)
            if RA==1:
                ax.axvline(x=xmax-a, color='k',linestyle=':',linewidth=1.5)
            if Bon==1:
                ax.axvline(x=B1, color='k',linestyle=':',linewidth=1.5)  
                ax.axvline(x=B2, color='k',linestyle=':',linewidth=1.5)
        ax.legend()

    

        
    
plt.show
#%% SAVE figures
import os


folder = pref+name+'/'
for i in range(Ns-2):
    folder = folder + name_species[i+1]
    
folder = folder+'('+name_species[-1]+')/'
folder = folder + f'Ntor-{Ntor:.0f}_Fre-{freq:.1f}_Con-{conc[-1]:.0f}'
print(folder)
    
if not os.path.exists(folder):
    os.makedirs(folder)   # create the folder if it does't exis

figname = ['Sum','Dens-Temp','Poyn-Elec','Pabs'] 
for i in plt.get_fignums():
    plt.figure(i)
    filename = f'{figname[i-1]}.pdf'
    filepath = os.path.join(folder, filename)
    plt.savefig(filepath,dpi=300)
    
    
    
