#!/bin/sh

# example for plotting H and C shifts of RNA/DNA
# copy this script to a working directory, make local changes as needed

# first, extract files with atomname, exp calc from rdb files:

row < xxxx.rdb atomname mat H | row atomname ne H3 | row atomname ne H1 | \
   jointbl res atomname xxxx_exp.rdb| column atomname exp calc | \
   headchg -del > H.dat

row < xxxx.rdb atomname mat C \
   jointbl res atomname xxxx_exp.rdb| column atomname exp calc | \
   headchg -del > C.dat

#  second, use these files as input to plotshifts.py:

$AFNMRHOME/src/afnmr/plotshifts.py -i H.dat  C.dat \
              -t '$^1$H Shifts' '$^{13}$C Shifts' \
              --linear-fit  --no-stats \
              -x "exp shifts (ppm)" \
              -y "calc shifts (ppm)" \
              --scale-factor 2.0 --legend-scale-factor 1.3 \
              --r-label 0.4 0.1 --rrmse --legend-columns=3 -o calc_vs_exp.pdf

