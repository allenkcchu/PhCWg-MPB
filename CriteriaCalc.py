#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 11:47:29 2022

@author: allenchu
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.constants import c, pi
from scipy import interpolate

plt.rcParams.update({'font.size':10})
plt.close('all')

FigPath = os.path.join(os.getcwd(),'s1shiftSweep','Figures')
DataPath = os.path.join(os.getcwd(),'s1shiftSweep','Data')

a = 396
r = 130
s2x = 0
s2y = 0
s1x = 85
s1y = -155

index = 3.1

nSiO2 = 1.44
nSRN = index

L = 10e-6                      # m
gamma0 = 550                   # 1/W/m
alpha0 = 0.6                   # dB/mm
alpha0 = alpha0/4.31*1e3       # 1/m

fig, ax = plt.subplots(figsize=(12,4),ncols=3)

dfb = pd.read_csv(os.path.join(DataPath,f'USRN_band_a{a}_r{r}_s1x{s1x}_s1y{s1y}_s2x{s2x}_s2y{s2y}_n{int(index*100)}.dat'))
dfb = dfb.apply(pd.to_numeric, errors='coerce')
dfb.dropna(axis=1,inplace=True)

selectedBand = 12
k = dfb[' k1'].values
# kn = np.linspace(k[0],k[-1],1001)
LL_SRN = k/nSRN
LL_SiO2 = k/nSiO2
band = dfb[f' te band {selectedBand}'].values
# f = interpolate.interp1d(k, band,kind='quadratic')
# bandn = f(kn)
ax[0].plot(k,LL_SRN,'k:')
ax[0].plot(k,LL_SiO2,'k')
radiationMode = (band>=LL_SiO2)
band[radiationMode] = np.nan
ax[0].plot(k,band,label=f'Band {selectedBand}')
# ax[0].plot(kn,bandn,'--')   
        
ax[0].set_xlim(0.35,0.5)
ax[0].set_ylim(a/1600,a/1500)
ax[0].set_xlabel('Wave Vector ($ka/2\pi$)')
ax[0].set_ylabel('Frequency ($\omega a/2\pi c$)')
ax[0].grid(True)


# dk = np.diff(kn)
# dw = np.diff(bandn)
# vg = (dk/dw)**-1
# wp = 0.5*(bandn[:-1]+bandn[1:])
# bandp = 0.5*(bandn[:-1]+bandn[1:])
# ax[1].plot(a/bandp,np.abs(1/vg))

dfv = pd.read_csv(os.path.join(DataPath,f'USRN_gvd_a{a}_r{r}_s1x{s1x}_s1y{s1y}_s2x{s2x}_s2y{s2y}_n{int(index*100)}.dat'),delimiter='Vector3',engine='python',
                  names=dfb.columns[5:].values)
for N, col in enumerate(dfv.columns):
        dfv[col] = dfv[col].str.extract(pat='<(.*),.*,.*>.*')
        
dfv = dfv.apply(pd.to_numeric, errors='coerce')
vg = dfv[f' te band {selectedBand}'].values
vg[radiationMode] = np.nan
ng = abs(1/vg)

ax[1].plot(a/band,ng)
ax[1].set_xlim(1530,1580)
ax[1].set_ylim(0,100)
ax[1].grid(True)
ax[1].set_xlabel('Wavelength (nm)')
ax[1].set_ylabel('Group index ($n_g$)')

f = c/(a/band*1e-9)
vginv = 1/(c*vg)
w = 2*pi*f
beta2 = np.diff(vginv)/np.diff(w)   # s^2/m
beta2 = beta2*1e24 # ps^2/m
D2 = beta2*(-0.7619) # ps/nm/m
wave = 0.5*(c/f[1:]+c/f[:-1])*1e9
ax[2].plot(wave,D2)
# ax[1].plot(a/band,1/vg,'--')
    
ax[2].set_xlim(1550,1560)
ax[2].set_ylim(-10000,5000)
ax[2].grid(True)
ax[2].set_xlabel('Wavelength (nm)')
ax[2].set_ylabel('Dispersion (ps/nm/m)')

S = (np.interp(1560,a/band[~np.isnan(band)],ng[~np.isnan(band)])/index)
gamma = gamma0*S**2
alpha = alpha0*S**2
GVM = D2*L
intRange = (wave<=1560) & (wave>=1550)
GVM = np.trapz(GVM[intRange],wave[intRange])
loss = alpha*S**2*L
print(f'L={L*1e6}um')
print(f'alpha={alpha0*1e-3:.3f}(1/mm)')
print(f'S=ng/n={S:.3f}')
print(f'alpha*L={alpha*L:.5f}')
print(f'gamma={gamma:.3f}(1/W/m)')
print(f'GVM={GVM:.3f}(ps)')
print(f'nonlinear/dispersion={gamma*L/GVM:.3f}(1/W/ps)')
print(f'loss/dispersion={alpha*L/GVM:.3f}(1/ps)')

fig.suptitle(f'a={a}nm, r={r}nm, nUSRN={index}\ns1x={s1x}nm, s1y={s1y}nm, s2x={s2x}nm, s2y={s2y}nm')
fig.tight_layout()

fig.savefig(f'Criteria_a{a}_r{r}_nCore{int(100*index)}_s1x{s1x}_s1y{s1y}_s2x{s2x}_s2y{s2y}.png',transparent=True)

# ax[1].grid(True)
# ax[1].set_xlabel('Wavelength (nm)')
# ax[1].set_ylabel('Group index ($n_g$)')

# ax[0].legend(ncol=3,fontsize=8)
# fig.tight_layout()
    
# # fig.savefig(f'USRNbandCalc_a{a}_r{r}_s1x{s1x}_s1y{s1y}_s2x{s2x}_s2y{s2y}_n{int(index*100)}.png',transparent=True)
# # fig.savefig(os.path.join(FigPath,f'USRNs2Shift_a{a}_r{r}_n{int(index*100)}_s1x{s1x}_s1y{s1y}_s2y{s2y}.png'),transparent=True)
    