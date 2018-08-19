#! /usr/bin/env python
from __future__ import print_function
import sys,os
import numpy as np
import argparse

sys.path.append( os.environ['AMBERHOME'] + '/AmberTools/src/xtalutil/Analysis')
import ReadAmberFiles as raf
print(raf.__file__)

#======================================================================#
#                                                                      #
# Convert cpptraj normal mode file (evecs) to nmd file for viewing of  #
# normal modes in VMD. Requires VMD 1.9.1 or higher.                   #
#                                                                      #
# Arguments:                                                           #
#    -n number of normal modes to include in the nmd file (default=20) #
#    -ipdb name of pdb file (atoms must correspond to the normal modes)#
#    -ievecs evecs normal mode file created by cpptraj                 #
#    -ofile name of output nmd file (default out.nmd)                  #
#                                                                      #
# http://www.csb.pitt.edu/ProDy/index.html#nmwiz                       #
# v1: pawel janowski                                                   #
#                                                                      #
#======================================================================#


def run(nmodes, ipdb_name, ievecs_name, ofile_name):
    # Get info from pdb
    ipdb = raf.pdb(ipdb_name)
    name = ipdb_name.rstrip('.pdb')
    atomnames = ipdb.Get_AtomNames()
    resnames = ipdb.Get_ResidueNames()
    resids = ipdb.Get_ResidueNumbers()
    bfactors = ipdb.Get_Bfactors()
    coordinates = ipdb.Get_Coords()

    # Write info from pdb to nmd file
    ofile = open(ofile_name, 'w')
    ofile.write('nmwiz_load %s\n' % ofile_name)
    print('name %s' % name, file=ofile)
    print('atomnames ', end="", file=ofile)
    print(*atomnames, file=ofile)
    print('resnames ', end="", file=ofile)
    print(*resnames, file=ofile)
    print('resids ', end="", file=ofile)
    print(*resids, file=ofile)
    print('bfactors ', end="", file=ofile)
    print(*bfactors, file=ofile)
    print('coordinates ', end="", file=ofile)
    print(*["%8.3f " % xyz for atom in coordinates for xyz in atom], end="", file=ofile)
    # This would've been the python2.x way of doing it:
    #~ for atom in coordinates:
    #~ for coordinate in atom:
    #~ ofile.write("%8.3f " %coordinate)
    print('', file=ofile)

    # Get modes from evecs file and write to nmd file
    vecfile = open(ievecs_name, 'r')
    line = ''
    while line[0:5] != ' ****':
        line = vecfile.readline()
    for i in range(nmodes):
        mode_number, freq = vecfile.readline().strip().split()
        assert int(mode_number) == i + 1
        mode = []
        line = vecfile.readline()
        while line[0:5] != ' ****':
            for i in line.strip().split():
                mode.append(float(i))
            line = vecfile.readline()
        amplitude = (1 / float(freq))
        print('mode %2d %12.10f ' % (int(mode_number), amplitude), end="", file=ofile)
        print(*["%12.5f " % fluc for fluc in mode], end="", file=ofile)
        print('', file=ofile)
    vecfile.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-n",
        help="number of normal modes to include (default=20)",
        default=20)
    parser.add_argument(
        "-ipdb",
        help="name of pdb file (atoms must correspond to normal modes")
    parser.add_argument(
        "-ievecs",
        help="name of evecs file (output of cpptraj diagmatrix)")
    parser.add_argument(
        "-ofile",
        help="name of nmd file to output (default=out.nmd)",
        default="out.nmd")
    args = parser.parse_args()
    run(int(args.n), args.ipdb, args.ievecs, args.ofile)
