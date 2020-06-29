#!/bin/sh 

#  Sample script to show how to run a single structure through
#     sander, invoking RISM.  Substitute "xxxx" with some labe for your system.

cat > mdin.rism.$1 <<EOF
  Single point calculation, invoking 3D-RISM
 &cntrl
    ntx=1, ntpr=1, ntwx=0, ntwr=0, ioutfm=1,
    imin=0, nstlim=0, 
    ntb=1,
    irism=1,
    cut=8.0,
 /
 &rism
    periodic='pme'
    closure='kh'
    buffer=1., grdspc=0.35,0.35,0.35,
    solvcut=8.0,
    verbose=2,
    npropagate=5,
    mdiis_del=0.7, mdiis_nvec=5, tolerance=1e-6,
    apply_rism_force=0,
    volfmt='dx',
    centering=0,
    ntwrism=1,
 /
EOF

sander -O -i mdin.rism.$1 -o xxxx.out \
    -p xxxx.parm7 -c xxxx.rst7 \
    -xvv $AMBERHOME/dat/rism1d/cSPCE_KH_NaCl_0.1M.xvv \
    -electronMap xxxx

/bin/rm mdin.rism.$1 restrt mdinfo
