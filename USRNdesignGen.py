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


parser = argparse.ArgumentParser()
parser.add_argument("-a", help="", type=float)
parser.add_argument("-r", help="", type=float)
parser.add_argument("-s1x", help="", type=float)
parser.add_argument("-s1y", help="", type=float)
parser.add_argument("-s2x", help="", type=float)
parser.add_argument("-s2y", help="", type=float)
parser.add_argument("-nSRN", help="", type=float)
args = parser.parse_args()

plt.close('all')

# Used for calculation of line defect bands
# For Tau-K direction Waveguide

a = 1               # although not used, can be defined as a symbol
# h = mp.inf        # for 2D
a_nm = args.a
h_nm = 300
h = h_nm / a_nm     # for 3D
r_nm = args.r
r = r_nm / a_nm
s1x_nm = args.s1x
s1y_nm = args.s1y
s2x_nm = args.s2x
s2y_nm = args.s2y
s1x = s1x_nm / a_nm
s1y = s1y_nm / a_nm
s2x = s2x_nm / a_nm
s2y = s2y_nm / a_nm

# n_SRN = 2.72
n_SRN = args.nSRN
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
k_interp = 41
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

    

ms.run_te(mpb.display_group_velocities)


md = mpb.MPBData(rectify=True)
eps = ms.get_epsilon()
converted_eps = md.convert(eps)
plt.figure();plt.imshow(converted_eps.T,interpolation='spline36',cmap='binary')

