/* Mark an arbitrary set of points as to solvent accessibility.
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

$Id: SAVtag_points.cc,v 2.4 2001/04/22 01:32:03 bashford Exp $
*/
#include "MEAD/SolvAccVol.h"
#include "MEAD/Shell.h"
#include "MEAD/TPShell.h"
#include "MEAD/Cone.h"
#include "MEAD/Sausage.h"
#include "MEAD/globals.h"
#include "MEAD/AccTag.h"

#include <iostream>

// FIXME
// Setting up the box and maybe also a series of sub-boxes could be
// done in the ctor.

void
SolvAccVolRep::tag_points(int npts, const Coord* pt, AccTag* acc_array)
{
  blab3 << "SolvAccVolRep::tag_points entered" << endl;
  if( !is_calculated ) {
    blab1 << "SolvAccVolRep::tag_points called but Anal rep not calculated"
      << "\nSo call anal_calc" << endl;
    anal_calc();
  }

  // Mark everthing exterior to start..
  int h;
  for(h = 0; h<npts; ++h) {
    acc_array[h] = exterior;
  }
  if (sph_count==0) {
    blab2 << "No Atoms present.  All points accessible." << endl;
    return;
  }

// Only points in a box around the molecule need be examined.
  int *idxmolbox = new int[npts];
  int in_molbox = 0;
  Coord udag(1,1,1);
  Coord mol_top_corner = atsph[0].get_coord()
    + atsph[0].get_outer_rad() * udag;
  Coord mol_bottom_corner = atsph[0].get_coord()
    - atsph[0].get_outer_rad() * udag;
  int at;
  for (at=1; at<sph_count; ++at) {
    Coord top_at = atsph[at].get_coord() + atsph[at].get_outer_rad() * udag;
    Coord bot_at = atsph[at].get_coord() - atsph[at].get_outer_rad() * udag;
    if (top_at.x > mol_top_corner.x) mol_top_corner.x = top_at.x;
    if (top_at.y > mol_top_corner.y) mol_top_corner.y = top_at.y;
    if (top_at.z > mol_top_corner.z) mol_top_corner.z = top_at.z;
    if (bot_at.x < mol_bottom_corner.x) mol_bottom_corner.x = bot_at.x;
    if (bot_at.y < mol_bottom_corner.y) mol_bottom_corner.y = bot_at.y;
    if (bot_at.z < mol_bottom_corner.z) mol_bottom_corner.z = bot_at.z;
  }
  blab3 << "SolvAccVolRep::tag_points finished pass through shells\n"
    << "mol_top_corner = " << mol_top_corner << ", mol_bottom_corner = "
      << mol_bottom_corner << "." << endl;

  for (h=0; h<npts; ++h) {
    if (pt[h] < mol_top_corner && pt[h] > mol_bottom_corner)
      idxmolbox[in_molbox++] = h;
  }
  blab3 << "SolvAccVolRep::tag_points: " << in_molbox
    << " points found within molecule corners." << endl;

  // Create an array of TPShells to do the Shell and Cone work.

  // Scan list of points to see how many are inside the sphere.
  int cs_count = 0;
  TPShell *tp_sh = new TPShell [sph_count];
  int *idxwork = new int[npts];
  for(at=0; at<sph_count; ++at) {
    // Define corners of atom covering box to be scan pts against.
    Coord coord = atsph[at].get_coord();
    float outer_rad = atsph[at].get_outer_rad();
    Coord diag(outer_rad,outer_rad,outer_rad);
    Coord bottom_corner = coord - diag;
    Coord top_corner = coord + diag;
    // Scan points and put indices of near ones in idxwork
    int num_inbox=0;
    for (int ih=0; ih<in_molbox; ++ih) {
      int h = idxmolbox[ih];
      if(pt[h].x >= bottom_corner.x)
	if(pt[h].x <= top_corner.x)
	  if(pt[h].y >= bottom_corner.y)
	    if(pt[h].y <= top_corner.y)
	      if(pt[h].z >= bottom_corner.z)
		if(pt[h].z <= top_corner.z) {
		  idxwork[num_inbox]=h;
		  ++num_inbox;
		}
    }
    if (num_inbox) {
      TPShell cs(atsph[at], num_inbox, idxwork);
      tp_sh[cs_count] = cs;
      ++cs_count;
    }
  }
  delete [] idxwork;

/*
  int cs_count = 0;
  TPShell *tp_sh = new TPShell [sph_count];
  for(int at=0; at<sph_count; ++at) {
    TPShell cs(atsph[at], npts, pt);
    if (cs.has_pts()) {
      tp_sh[cs_count] = cs;
      ++cs_count;
    }
    blab3 << "at = " << at << endl;
  }
*/

  // On first pass through shells, mark atom interiors as interior.
  int ci;
  for (ci=0; ci<cs_count; ++ci)
    tp_sh[ci].mark_by_radii(pt, acc_array);

  if (blab3pt != &cnull) {
    // Do some unnecessary work to report statistics.
    int nodeinshell = 0;
    int nodeinatom = 0;
    for(h=0; h<npts; h++ ) {
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
  for (ci=0; ci<cs_count; ++ci)
    tp_sh[ci].mark_by_cones(pt, acc_array);

  delete [] tp_sh;

  if (blab3pt != &cnull) {
    // Do some unnecessary work to report statistics.
    int nodeincone = 0;
    for(h=0; h<npts; h++ ) {
      if( acc_array[h] == undecided )
	nodeincone++;
    }
    blab3 << "Total number of nodes inside cones is " << nodeincone << endl;
  }

  // Now pass through tube neighborhoods.
  for (h=0; h<npts; ++h)
    if (acc_array[h] == undecided)
      for(int s1=0; s1<sausage_count; s1++ )
	if (sglist[s1].pt_inside(pt[h])) {
	  acc_array[h] = in_tube;
	  break;
	}

  blab2 << "Finishing checking sausages " << endl;

  // Mark still-undecided nodes as interior and in-tube nodes as exterior
  int in_sphere_count = 0;
  int in_cone_count = 0;
  int in_cone_acc = 0;
  for(h=0; h<npts; h++ ) {
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
  blab2 << "Number of points inside cones but water acc. becuase in tubes is "
    << in_cone_acc << endl;
  blab2 << "Total number of nodes water inaccessible is "
    << in_sphere_count + in_cone_count << endl;

}


int
SolvAccVolRep::accessible(const Coord& r)
{
  if( !is_calculated ) {
    blab1 << "SolvAccVolRep::tag_points called but Anal rep not calculated"
      << "\nSo call anal_calc" << endl;
    anal_calc();
  }

  AccTag acc = exterior;

  // FIXME.  The use of TPShell to do these one-point functions below is
  // awkward and slow.  There should be one-point functions of Shell.

  // Create an array of TPShells to do the Shell and Cone work.
  // Kluge around lack of one-point functions using arrays of one element.
  Coord pt[1] = {r};
  AccTag acc_array[1];
  acc_array[0] = acc;
  int cs_count = 0;
  TPShell *tp_sh = new TPShell [sph_count];
  for(int at=0; at<sph_count; ++at) {
    TPShell cs(atsph[at], 1, pt);
    if (cs.has_pts()) {
      tp_sh[cs_count] = cs;
      ++cs_count;
    }
  }

  // On first pass through shells, mark atom interiors as interior.
  int ci;
  for (ci=0; ci<cs_count; ++ci)
    tp_sh[ci].mark_by_radii(pt, acc_array);

  // On second pass, cone-free regions of shells will be marked exterior
  // and cone-covered regions of shells will be marked undecided.
  for (ci=0; ci<cs_count; ++ci)
    tp_sh[ci].mark_by_cones(pt, acc_array);

  delete [] tp_sh;
  acc = acc_array[0];

  // Now pass through tube neighborhoods.
  if (acc == undecided)
    for(int s1=0; s1<sausage_count; s1++ )
      if (sglist[s1].pt_inside(r)) {
	acc = in_tube;
	break;
      }
  if (acc == undecided)
    acc = interior;
  else if (acc == in_tube)
    acc = exterior;
  return (acc == exterior);
}

// SAVtag_points.cc ends here
