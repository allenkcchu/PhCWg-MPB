#!/bin/bash
python USRNdesignGen.py -a 410 -r 120 -s1x 0 -s1y 0 -s2x 0 -s2y 0 -nSRN 3.1 > SingleTest/Data/USRN_a410_r120_s1x0_s1y0_s2x0_s2y0_n310.out;
grep freqs SingleTest/Data/USRN_a410_r120_s1x0_s1y0_s2x0_s2y0_n310.out > SingleTest/Data/USRN_band_a410_r120_s1x0_s1y0_s2x0_s2y0_n310.dat;
grep velocity SingleTest/Data/USRN_a410_r120_s1x0_s1y0_s2x0_s2y0_n310.out > SingleTest/Data/USRN_gvd_a410_r120_s1x0_s1y0_s2x0_s2y0_n310.dat;
