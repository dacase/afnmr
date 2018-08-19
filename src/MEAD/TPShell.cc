/* Variant of Shell for dealing with arbitrary sets of points to mark.
    Copyright (c) 1993--1995 by Donald Bashford.

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

$Id: TPShell.cc,v 1.6 2004/12/06 17:58:42 bashford Exp $
*/

#include "MEAD/TPShell.h"

#include "MEAD/AccTag.h"
#include "MEAD/Shell.h"
#include "MEAD/Cone.h"
#include "MEAD/globals.h"

// This Ctor installs a pre-defined index
TPShell::TPShell(const Shell& sph, int nidx, const int* idx)
: coord(sph.get_coord()), inner_rad(sph.get_inner_rad()),
  outer_rad(sph.get_outer_rad()), outer_rad_sq(sph.get_outer_rad_sq()),
  num_inbox(nidx)
{
  if (num_inbox==0)
    return;

  inbox_pt_index = new int[num_inbox];
  for (int k=0; k<num_inbox; ++k)
    inbox_pt_index[k] = idx[k];

  if (sph.is_buried()) flag=buried;
  else if (sph.is_partially_buried()) flag = partially_buried;
  else if (sph.is_free()) flag = free;
  else ::error("ERROR CLShell ctor could not determine flag value");

  numcones = sph.get_numcones();
  if (numcones) {
    cone = new Cone [numcones];
    int ic = 0;
    for (ShellConeIterator sci(sph); sci.valid(); sci.next(), ++ic) {
      const Cone* oc = sci.cur_cone_pt();
      cone[ic].cos_ang1_sq = oc->cos_ang1_sq;
      cone[ic].unit_axis = oc->unit_axis;
    }
    if (ic != numcones)
      ::error ("TPShell constructor: cone number mismatch");
  }
  else
    cone = 0;
}



TPShell::TPShell(const Shell& sph, int npts, const Coord* pts)
: coord(sph.get_coord()), inner_rad(sph.get_inner_rad()),
  outer_rad(sph.get_outer_rad()), outer_rad_sq(sph.get_outer_rad_sq())
{
 // Zero some things at first, in case all pts are far away.
  numcones=0;
  inbox_pt_index = 0;
  cone = 0;

  if (sph.is_buried()) flag=buried;
  else if (sph.is_partially_buried()) flag = partially_buried;
  else if (sph.is_free()) flag = free;
  else ::error("ERROR CLShell ctor could not determine flag value");

  // Define corners of atom covering box to be scan pts against.
  Coord diag(outer_rad,outer_rad,outer_rad);
  Coord bottom_corner = coord - diag;
  Coord top_corner = coord + diag;

  // Scan list of points to see how many are inside the sphere.
  num_inbox=0;
  int *temparr = new int[npts]; // temparr contents will become inbox_pt_index
  for (int h=0; h<npts; ++h) {
    if(pts[h].x >= bottom_corner.x)
      if(pts[h].x <= top_corner.x)
	if(pts[h].y >= bottom_corner.y)
	  if(pts[h].y <= top_corner.y)
	    if(pts[h].z >= bottom_corner.z)
	      if(pts[h].z <= top_corner.z) {
		temparr[num_inbox]=h;
		++num_inbox;
	      }
  }
  if (num_inbox==0) {
    delete [] temparr;
    return;
  }

  // Copy the temp array into inbox_pt_index;
  inbox_pt_index = new int[num_inbox];
  for (int k=0; k<num_inbox; ++k)
    inbox_pt_index[k] = temparr[k];
  delete [] temparr;
  // Add the cones
  numcones = sph.get_numcones();
  if (numcones) {
    cone = new Cone [numcones];
    int ic = 0;
    for (ShellConeIterator sci(sph); sci.valid(); sci.next(), ++ic) {
      const Cone* oc = sci.cur_cone_pt();
      cone[ic].cos_ang1_sq = oc->cos_ang1_sq;
      cone[ic].unit_axis = oc->unit_axis;
    }
    if (ic != numcones)
      ::error ("TPShell constructor: cone number mismatch");
  }
  else
    cone = 0;
}

TPShell::TPShell(const TPShell& old)
: coord(old.coord), inner_rad(old.inner_rad),
  outer_rad(old.outer_rad), outer_rad_sq(old.outer_rad_sq),
  flag(old.flag), num_inbox(old.num_inbox), numcones(old.numcones)
{
  if (num_inbox>0) {
    inbox_pt_index = new int[num_inbox];
    for(int h=0; h<num_inbox; ++h)
      inbox_pt_index[h] = old.inbox_pt_index[h];
  }
  else
    inbox_pt_index = 0;
  if (numcones>0) {
    cone = new Cone[numcones];
    for (int i=0; i<numcones; ++i)
      cone[i] = old.cone[i];
  }
  else
    cone = 0;
}

TPShell&
TPShell::operator= (const TPShell& old)
{
  coord = old.coord;
  inner_rad = old.inner_rad;
  outer_rad = old.outer_rad;
  outer_rad_sq = old.outer_rad_sq;
  flag = old.flag;

  if (num_inbox) delete [] inbox_pt_index;
  num_inbox = old.num_inbox;
  if (num_inbox>0) {
    inbox_pt_index = new int[num_inbox];
    for(int h=0; h<num_inbox; ++h)
      inbox_pt_index[h] = old.inbox_pt_index[h];
  }
  else
    inbox_pt_index = 0;

  if (numcones) delete [] cone;
  numcones = old.numcones;
  if (numcones) {
    cone = new Cone[numcones];
    for (int i=0; i<numcones; ++i)
      cone[i] = old.cone[i];
  }
  else
    cone = 0;
  return *this;
}

TPShell::~TPShell()
{
  if (numcones && cone) delete [] cone;
  if (num_inbox) delete [] inbox_pt_index;
}

void TPShell::mark_by_radii(const Coord *pt, AccTag *acc_array) const
{
  if (num_inbox==0) return;
  float inner_rsq = inner_rad*inner_rad;
  if (flag==free) {
    for (int h=0; h<num_inbox; ++h) {
      int idx = inbox_pt_index[h];
      if( acc_array[idx] != interior ) {
	Coord r = pt[idx]-coord;
	float rsq = r*r;
	if (rsq<inner_rsq)
	  acc_array[idx] = interior;
	else if (rsq<outer_rad_sq)
	  acc_array[idx] = exterior;
      }
    }
  }
  else {
    for (int h=0; h<num_inbox; ++h) {
      int idx = inbox_pt_index[h];
      if( acc_array[idx] != interior ) {
	Coord r = pt[idx]-coord;
	float rsq = r*r;
	if (rsq<inner_rsq)
	  acc_array[idx] = interior;
	else if (rsq<outer_rad_sq)
	  acc_array[idx] = undecided;
      }
    }
  }
}




void TPShell::mark_by_cones(const Coord *pt, AccTag *acc_array) const
{
  if (num_inbox==0) return;
  if (flag == buried || flag == free) return;  // buried or free atom case

  for (int h=0; h<num_inbox; ++h) {
    int idx = inbox_pt_index[h];
    if (acc_array[idx] == undecided) {
      Coord r = pt[idx]-coord;
      float rsq = r*r;
      if (rsq<=outer_rad_sq) {
	int inside_cone = 0;
	for( int nc=0; nc<numcones; nc++) {
	  float dot_prod = r * cone[nc].unit_axis;
	  float cone_cos_ang1_sq = cone[nc].cos_ang1_sq;
	  if( cone_cos_ang1_sq >= 0 ) {
	    if( dot_prod > 0 ) {
	      float pt_cos_sq = dot_prod*dot_prod / rsq;
	      if (pt_cos_sq > cone_cos_ang1_sq) {
		inside_cone = 1;
		break;
	      }
	    }
	  }
	  else {  // A negative cone cosine sq means ang > pi/2
	    if( dot_prod < 0 ) {
	      float pt_cos_sq = dot_prod*dot_prod / rsq;
	      if (pt_cos_sq < -cone_cos_ang1_sq) {
		inside_cone = 1;
		break;
	      }
	    }
	    else {
	      inside_cone = 1;
	      break;
	    }
	  }
	}
	if (inside_cone==0)
	  acc_array[idx] = exterior;
      }
    }
  }

}

// TPShell.cc ends here
