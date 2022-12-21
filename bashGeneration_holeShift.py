#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 00:06:11 2022

@author: allenchu
"""

import numpy as np
import os

def genScript(path,a,r,s1x,s1y,s2x,s2y,nSRN):
    line1 = f'python USRNdesignGen.py -a {a} -r {r} -s1x {s1x} -s1y {s1y} -s2x {s2x} -s2y {s2y} -nSRN {nSRN} > {path}/USRN_a{a}_r{r}_s1x{s1x}_s1y{s1y}_s2x{s2x}_s2y{s2y}_n{int(100*nSRN)}.out;\n'
    line2 = f'grep freqs {path}/USRN_a{a}_r{r}_s1x{s1x}_s1y{s1y}_s2x{s2x}_s2y{s2y}_n{int(100*nSRN)}.out > {path}/USRN_band_a{a}_r{r}_s1x{s1x}_s1y{s1y}_s2x{s2x}_s2y{s2y}_n{int(100*nSRN)}.dat;\n'
    line3 = f'grep velocity {path}/USRN_a{a}_r{r}_s1x{s1x}_s1y{s1y}_s2x{s2x}_s2y{s2y}_n{int(100*nSRN)}.out > {path}/USRN_gvd_a{a}_r{r}_s1x{s1x}_s1y{s1y}_s2x{s2x}_s2y{s2y}_n{int(100*nSRN)}.dat;\n'

    return line1+line2+line3
path = os.path.join('s2shiftSweep','Data')
filename = 's2shiftSweep'
a = 410
r = 120
s1x = 200
s1y = -40
s2xList = np.linspace(0,200,11)
s2yList = np.linspace(140,140,1)
nSRN = 2.7

with open(f'{filename}.sh','w') as test:
    test.write('#!/bin/bash\n')
    for s2y in s2yList:
        for s2x in s2xList:
            s2x = int(s2x)
            s2y = int(s2y)
            line = genScript(path,a,r,s1x,s1y,s2x,s2y,nSRN)
            test.write(line)