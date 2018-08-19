#!/bin/sh

#  Sample script to show how to run a single structure through
#     rism3d.snglpnt.  Substitute "xxxx" with some labe for your system.

#  100 mM NaCl in cSPCE water:  (change as needed)
solvent=$AMBERHOME/dat/rism1d/cSPCE_KH_NaCl_0.1M.xvv

rism3d.snglpnt --pdb xxxx.amber.pdb --prmtop xxxx.parm7 \
               --xvv $solvent --tolerance 1e-6 --verbose 2 \
               --guv  xxxx  --mdiis_nvec 10 --buffer 20.0  > xxxx.r3d
