#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 06:27:00 2022

@author: allenchu
"""

import math
import meep as mp
import numpy as np
from meep import mpb
import matplotlib.pyplot as plt
import argparse
import os
import pickle

efields = []
dfields = []

def get_efields(ms, band):
    efields.append(ms.get_efield(band, bloch_phase=True))
    
def get_dfields(ms, band):
    dfields.append(ms.get_dfield(band, bloch_phase=True))

# parser = argparse.ArgumentParser()
# parser.add_argument("-a", help="", type=float)
# parser.add_argument("-r", help="", type=float)
# parser.add_argument("-s1x", help="", type=float)
# parser.add_argument("-s1y", help="", type=float)
# parser.add_argument("-s2x", help="", type=float)
# parser.add_argument("-s2y", help="", type=float)
# parser.add_argument("-nSRN", help="", type=float)
a_nm = 396
r_nm = 130
s1x_nm = 85
s1y_nm = -155
s2x_nm = 0
s2y_nm = 0
n_SRN = 3.1
# args = parser.parse_args()

plt.close('all')

# Used for calculation of line defect bands
# For Tau-K direction Waveguide

a = 1               # although not used, can be defined as a symbol
# h = mp.inf        # for 2D

h_nm = 300
h = h_nm / a_nm     # for 3D
r = r_nm / a_nm
s1x = s1x_nm / a_nm
s1y = s1y_nm / a_nm
s2x = s2x_nm / a_nm
s2y = s2y_nm / a_nm

# n_SRN = 2.72
n_SiO2 = 1.44
SRN = mp.Medium(epsilon=n_SRN**2)
Glass = mp.Medium(epsilon=n_SiO2**2)

slab_si = mp.Block(center=mp.Vector3(0,0,0), material=SRN, size=mp.Vector3(mp.inf, mp.inf, h))
cyl_air = mp.Cylinder(center=mp.Vector3(0,0,0), height=h, radius=r, material=Glass)

# Define the positions of the system
wg_width = 1                        # 1 means w1 waveguide
vec_1 = mp.Vector3(0.5 , 0.5*3**0.5)          # vector of tri-lattice
vec_2 = mp.Vector3(0.5 ,-0.5*3**0.5)          

# Define a position shift, because it won't work when the structure is sym.
# pos_shift = mp.Vector3(0.1, 0, 0)

# Define positions of holes, u means up
vec_h1u = mp.Vector3(0, 0.5*3**0.5*wg_width, 0)
vec_h1u2 = vec_h1u + mp.Vector3(s1x,s1y,0)
vec_h2u = vec_h1u + vec_1 + mp.Vector3(s2x,s2y,0)
vec_h3u = vec_h1u - vec_2 + mp.Vector3(s2x,s2y,0)
vec_h4u = vec_h1u - vec_2 + vec_1
vec_h5u = vec_h4u + vec_1
vec_h6u = vec_h4u - vec_2
vec_h7u = vec_h6u + vec_1
vec_h1d = mp.Vector3(0, -0.5*3**0.5*wg_width, 0)
vec_h1d2 = vec_h1d + mp.Vector3(s1x,-s1y,0)
vec_h2d = vec_h1d + vec_2 + mp.Vector3(s2x,-s2y,0)
vec_h3d = vec_h1d - vec_1 + mp.Vector3(s2x,-s2y,0)
vec_h4d = vec_h1d - vec_1 + vec_2
vec_h5d = vec_h4d + vec_2
vec_h6d = vec_h4d - vec_1
vec_h7d = vec_h6d + vec_2

# define geometry before run
default_material = Glass
geometry_lattice = mp.Lattice(size=mp.Vector3(1,vec_h7u.y-vec_h7d.y))

geometry = list()
geometry.append(slab_si)
geometry.append(mp.geometric_object_duplicates(shift_vector=vec_h1u2,min_multiple=1,max_multiple=1,go=cyl_air)[0])
geometry.append(mp.geometric_object_duplicates(shift_vector=vec_h2u,min_multiple=1,max_multiple=1,go=cyl_air)[0])
geometry.append(mp.geometric_object_duplicates(shift_vector=vec_h3u,min_multiple=1,max_multiple=1,go=cyl_air)[0])
geometry.append(mp.geometric_object_duplicates(shift_vector=vec_h4u,min_multiple=1,max_multiple=1,go=cyl_air)[0])
geometry.append(mp.geometric_object_duplicates(shift_vector=vec_h5u,min_multiple=1,max_multiple=1,go=cyl_air)[0])
geometry.append(mp.geometric_object_duplicates(shift_vector=vec_h5u,min_multiple=1,max_multiple=1,go=cyl_air)[0])
geometry.append(mp.geometric_object_duplicates(shift_vector=vec_h6u,min_multiple=1,max_multiple=1,go=cyl_air)[0])
geometry.append(mp.geometric_object_duplicates(shift_vector=vec_h7u,min_multiple=1,max_multiple=1,go=cyl_air)[0])
geometry.append(mp.geometric_object_duplicates(shift_vector=vec_h1d2,min_multiple=1,max_multiple=1,go=cyl_air)[0])
geometry.append(mp.geometric_object_duplicates(shift_vector=vec_h2d,min_multiple=1,max_multiple=1,go=cyl_air)[0])
geometry.append(mp.geometric_object_duplicates(shift_vector=vec_h3d,min_multiple=1,max_multiple=1,go=cyl_air)[0])
geometry.append(mp.geometric_object_duplicates(shift_vector=vec_h4d,min_multiple=1,max_multiple=1,go=cyl_air)[0])
geometry.append(mp.geometric_object_duplicates(shift_vector=vec_h5d,min_multiple=1,max_multiple=1,go=cyl_air)[0])
geometry.append(mp.geometric_object_duplicates(shift_vector=vec_h5d,min_multiple=1,max_multiple=1,go=cyl_air)[0])
geometry.append(mp.geometric_object_duplicates(shift_vector=vec_h6d,min_multiple=1,max_multiple=1,go=cyl_air)[0])
geometry.append(mp.geometric_object_duplicates(shift_vector=vec_h7d,min_multiple=1,max_multiple=1,go=cyl_air)[0])


# Define parameters to run
resolution = mp.Vector3(50,50,20)
num_bands = 12
mesh_size = 3

k_points = [mp.Vector3(0.3,0,0),
            mp.Vector3(0.5,0,0)]
k_interp = 39
k_points = mp.interpolate(k_interp, k_points)

print(f"a-nm: {a_nm}") 
print(f"r: {r}, h: {h}, wg-width: {wg_width}")
print(f"resolution: {resolution}")

# run the simulations
ms = mpb.ModeSolver(
    default_material=default_material,
    geometry_lattice=geometry_lattice,
    geometry=geometry,
    k_points=k_points,
    num_bands=num_bands,
    resolution=resolution,
    mesh_size=mesh_size
)





kList = [0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44,0.46,0.48,0.5]
# kList = [0.3,0.4,0.5]

for kNum, k in enumerate(kList):
    ms.run_te(mpb.output_at_kpoint(mp.Vector3(k,0,0), mpb.fix_efield_phase, get_efields, get_dfields))
    
    
md = mpb.MPBData(rectify=True, resolution=32, periods=1)
converted = []
for modeNum, f in enumerate(efields):
    # Get just the z component of the efields
    if modeNum in np.array([12])*(np.arange(len(kList))+1)-1:
        # f = f[:, :, 0, 0]
        f = np.squeeze(np.sum(np.abs(f)**2,axis=3)**0.5)
        converted.append(md.convert(f))

eps = ms.get_epsilon()
converted_eps = md.convert(eps)

fig, ax = plt.subplots(nrows=1,ncols=len(kList),figsize=(12,3))
for i, f in enumerate(converted):
    ax[i].contour(converted_eps.T, cmap='binary')
    ax[i].imshow(np.abs(f).T, interpolation='spline36', cmap='jet', alpha=0.9) # RdBu
    # ax[i].set_ylim(int(np.size(f,1)*0.35),int(np.size(f,1)*0.65))
    ax[i].axis('off')
    ax[i].set_title(f'k = {kList[i]}')

fig.suptitle(f'a={a_nm}nm, r={r_nm}nm, nUSRN={n_SRN}, nSiO2={n_SiO2}\ns1x={s1x_nm}nm, s1y={s1y_nm}nm, s2x={s2x_nm}nm, s2y={s2y_nm}nm')
fig.tight_layout()
data = [dfields,efields,converted_eps,a_nm,r_nm]
with open(os.path.join('LossCalc','fieldTest'),'wb') as fp:
    pickle.dump(data,fp)
# fig.savefig(f'Efield_a{a_nm}_r{r_nm}nm_nCore{int(100*n_SRN)}_nCladding{int(100*n_SiO2)}_s1x{s1x_nm}_s1y{s1y_nm}_s2x{s2x_nm}_s2y{s2y_nm}.png',transparent=True)

# plt.figure();plt.imshow(converted_eps.T,interpolation='spline36',cmap='binary')

