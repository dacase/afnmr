/* Mark a cubic lattice as to solvent accessibility.
    Copyright (c) 1993--1995 by Donald Bashford and Tony You.

    This source code file is part of the MEAD (Macroscopic
    Electrostatics with Atomic Detail) package of objects and
    programs.  MEAD is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 1, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; see the file COPYING.  If not, write to
    the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA
    02139, USA.

    Donald Bashford can be contacted by electronic mail by the address,
    bashford@scripps.edu, or by paper mail at Department of Molecular
    Biology, The Scripps Research Institute, 10666 North Torrey Pines
    Road, La Jolla, California 92037.
$Id: SAVcalc_cuberep.cc,v 2.3 2001/04/22 01:32:03 bashford Exp $
*/
#include "MEAD/SolvAccVol.h"
#include "MEAD/Shell.h"
#include "MEAD/Cone.h"
#include "MEAD/Sausage.h"
#include "MEAD/globals.h"
#include "MEAD/CubeLatSpec.h"
#include "MEAD/AccTag.h"
#include "MEAD/CLShell.h"

#include <iostream>

void
SolvAccVolRep::calc_cuberep(const CubeLatSpec& cls, AccTag* acc_array)
{
  if( !is_calculated ) {
    blab1 << "SolvAccVolRep::calc_cuberep called but Anal rep not calculated"
      << "\nSo call anal_calc" << endl;
    anal_calc();
  }

  int grid_dim = cls.get_grid_dim();
  int nsq = grid_dim*grid_dim;
  int ncube = nsq*grid_dim;

  // Mark everthing exterior to start..
  int h;
  for(h = 0; h<ncube; ++h) {
    acc_array[h] = exterior;
  }

  if (sph_count != 0) {  // if no atoms, all is exterior.

    // Create an array of CLShells to do the Shell and Cone work.
    int cs_count = 0;
    CLShell *cl_sh = new CLShell [sph_count];
    for(int at=0; at<sph_count; ++at) {
      CLShell cs(atsph[at],cls);
      if (cs.in_lattice()) {
	cl_sh[cs_count] = cs;
	++cs_count;
      }
    }


    // On first pass through shells, mark atom interiors as interior.
    int ci;
    for (ci=0; ci<cs_count; ++ci) {
      if (cl_sh[ci].in_lattice())
	cl_sh[ci].mark_by_radii(acc_array);
    }

    if (blab3pt != &cnull) {
      // Do some unnecessary work to report statistics.
      int nodeinshell = 0;
      int nodeinatom = 0;
      for(h=0; h<ncube; h++ ) {
	AccTag ea = acc_array[h];
	if(ea == undecided )
	  nodeinshell++;
	else if (ea == interior)
	  nodeinatom++;
      }
      blab3 << "Total number of nodes inside atoms is " << nodeinatom << endl;
      blab3 << "Total number of nodes inside shells is " << nodeinshell
	<< endl;
    }

    // On second pass, cone-free regions of shells will be marked exterior
    // and cone-covered regions of shells will be marked undecided.
    for (ci=0; ci<cs_count; ++ci) {
      if (cl_sh[ci].in_lattice())
	cl_sh[ci].mark_by_cones(acc_array);
    }

    delete [] cl_sh;

    if (blab3pt != &cnull) {
      // Do some unnecessary work to report statistics.
      int nodeincone = 0;
      for(h=0; h<ncube; h++ ) {
	if( acc_array[h] == undecided )
	  nodeincone++;
      }
      blab3 << "Total number of nodes inside cones is " << nodeincone << endl;
    }

    // Now pass through tube neighborhoods.
    for( int s1=0; s1<sausage_count; s1++ ) {
      sglist[s1].mark_cubelat(cls, acc_array);
    }
    blab2 << "Finishing checking sausages " << endl;

    // Mark still-undecided nodes as interior and in-tube nodes as exterior
    int in_sphere_count = 0;
    int in_cone_count = 0;
    int in_cone_acc = 0;
    for(h=0; h<ncube; h++ ) {
      if( acc_array[h] == interior ) {
	in_sphere_count++;
      }
      else
	if( acc_array[h] == undecided  ) {
          in_cone_count++;
          acc_array[h] = interior;
	}
	else
	  if( acc_array[h] == in_tube) {
            in_cone_acc++;
            acc_array[h] = exterior;
	  }
    }

    blab2 << "Number of points inside atoms and water inacc. is "
      << in_sphere_count;
    blab2 << "Number of points inside cones and water inacc. is "
      << in_cone_count << endl;
    blab2<< "Number of points inside cones but water acc. becuase in tubes is "
      << in_cone_acc << endl;
    blab2 << "Total number of nodes water inaccessible is "
      << in_sphere_count + in_cone_count << endl;
  }
  else {
    blab2 << "No Atoms present.  All node points accessible." << endl;
  }
}

// SAVcalc_cuberep.cc ends here
