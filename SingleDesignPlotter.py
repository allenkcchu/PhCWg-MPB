#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 10:09:59 2022

@author: allenchu
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.constants import c, pi

plt.rcParams.update({'font.size':12})
# plt.close('all')

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

selectedBand = range(1,13)
# selectedBand = [11]
C = plt.cm.jet(np.linspace(0,1,len(selectedBand)))

fig, ax = plt.subplots(figsize=(12,4),ncols=2)

dfb = pd.read_csv(os.path.join(DataPath,f'USRN_band_a{a}_r{r}_s1x{s1x}_s1y{s1y}_s2x{s2x}_s2y{s2y}_n{int(index*100)}.dat'))
dfb = dfb.apply(pd.to_numeric, errors='coerce')
dfb.dropna(axis=1,inplace=True)

ax[0].plot(dfb[' k1'],dfb[' k1']/nSRN,'k--')
ax[0].plot(dfb[' k1'],dfb[' k1']/nSiO2,'k')
for N in range(np.size(dfb,1)-5):
    dfb[f' te band {N+1}'].iloc[dfb[f' te band {N+1}']>=dfb[' k1']/nSiO2] = np.nan
    if N in range(len(selectedBand)):
        ax[0].plot(dfb[' k1'],dfb[f' te band {N+1}'],color=C[N],label=f'Band {selectedBand[N]}')
       
        
ax[0].set_xlim(0.35,0.5)
ax[0].set_ylim(a/1600,a/1500)
ax[0].set_xlabel('Wave Vector ($ka/2\pi$)')
ax[0].set_ylabel('Frequency ($\omega a/2\pi c$)')

dfv = pd.read_csv(os.path.join(DataPath,f'USRN_gvd_a{a}_r{r}_s1x{s1x}_s1y{s1y}_s2x{s2x}_s2y{s2y}_n{int(index*100)}.dat'),delimiter='Vector3',engine='python',
                  names=dfb.columns[5:].values)
for N, col in enumerate(dfv.columns):
        dfv[col] = dfv[col].str.extract(pat='<(.*),.*,.*>.*')
        
dfv = dfv.apply(pd.to_numeric, errors='coerce')
dfv.dropna(axis=1,inplace=True)

# fig, ax = plt.subplots(figsize=(12,8))
for N, col in enumerate(dfv.columns):
    if N in range(len(selectedBand)):
        ax[1].plot(a/dfb[col],np.abs(1/dfv[col]),color=C[N])
    
ax[1].set_xlim(1500,1600)
ax[1].set_ylim(0,100)
ax[1].grid(True)
ax[1].set_xlabel('Wavelength (nm)')
ax[1].set_ylabel('Group index ($n_g$)')

ax[0].legend(ncol=3,fontsize=8)
fig.tight_layout()
    
# fig.savefig(f'USRNbandCalc_a{a}_r{r}_s1x{s1x}_s1y{s1y}_s2x{s2x}_s2y{s2y}_n{int(index*100)}.png',transparent=True)
# fig.savefig(os.path.join(FigPath,f'USRNs2Shift_a{a}_r{r}_n{int(index*100)}_s1x{s1x}_s1y{s1y}_s2y{s2y}.png'),transparent=True)
    