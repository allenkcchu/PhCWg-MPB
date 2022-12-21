#!/bin/bash
python USRNdesignGen.py -a 410 -r 120 -s1x 200 -s1y -40 -s2x 0 -s2y 140 -nSRN 2.7 > aSweep/Data/USRN_a410_r120_s1x200_s1y-40_s2x0_s2y140_n270.out;
grep freqs aSweep/Data/USRN_a410_r120_s1x200_s1y-40_s2x0_s2y140_n270.out > aSweep/Data/USRN_band_a410_r120_s1x200_s1y-40_s2x0_s2y140_n270.dat;
grep velocity aSweep/Data/USRN_a410_r120_s1x200_s1y-40_s2x0_s2y140_n270.out > aSweep/Data/USRN_gvd_a410_r120_s1x200_s1y-40_s2x0_s2y140_n270.dat;
