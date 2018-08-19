/* A varient of the Shell class specialized for cubic lattices.
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

$Id: CLShell.cc,v 1.14 2004/12/06 17:58:06 bashford Exp $
*/

#include "MEAD/CLShell.h"

#include "MEAD/AccTag.h"
#include "MEAD/CubeLatSpec.h"
#include "MEAD/CLShell.h"
#include "MEAD/Shell.h"
#include "MEAD/Cone.h"


CLShell::CLShell (const Shell& sph, const CubeLatSpec& cls)
{
  // Some of these prelim variables could be static members of the
  // class to save space and/or time.  But this would require measure to
  // insure only one set of CLShells is in use (existence) at a time.

  numcones = 0;  // Destructor-safe initial value in case out of lattice

  float spacing = cls.get_spacing();
  Coord grid_center_in_space = cls.get_center();
  grid_dim = cls.get_grid_dim();
  float grlen = (float) (grid_dim - 1);
  float halfgrlen = grlen/2;
  Coord grid_center_in_grid (halfgrlen, halfgrlen, halfgrlen);

  if (sph.is_buried()) flag=buried;
  else if (sph.is_partially_buried()) flag = partially_buried;
  else if (sph.is_free()) flag = free;
  else ::error("ERROR CLShell ctor could not determine flag value");
  inner_rad = sph.get_inner_rad() / spacing;
  outer_rad = sph.get_outer_rad() / spacing;
  outer_rad_sq = sph.get_outer_rad_sq() / (spacing*spacing);
  coord = (sph.get_coord() - grid_center_in_space)/spacing
    + grid_center_in_grid;

  // Now define corners of atom covering box to be scanned in grid-base coords.

  Coord diag(outer_rad,outer_rad,outer_rad);
  Coord bottom_corner = coord - diag;
  Coord top_corner = coord + diag;

  is_in = 1;
  if( bottom_corner.x>=0 && bottom_corner.x<grid_dim)
    i1 = (int) bottom_corner.x;
  else
    if( bottom_corner.x < 0 )
      i1 = 0;
    else
      {is_in = 0; return;}
  if( bottom_corner.y>=0 && bottom_corner.y<grid_dim)
    j1 = (int) bottom_corner.y;
  else
    if( bottom_corner.y < 0 )
      j1 = 0;
    else
      {is_in = 0; return;}
  if( bottom_corner.z>=0 && bottom_corner.z<grid_dim)
    k1 = (int) bottom_corner.z;
  else
    if( bottom_corner.z < 0 )
      k1 = 0;
    else
      {is_in = 0; return;}
  if( top_corner.x>=0 && top_corner.x<grid_dim-1 )
    i2 = (int) top_corner.x+1;
  else
    if( top_corner.x >= grid_dim-1 )
      i2 = grid_dim-1;
    else
      {is_in = 0; return;}
  if( top_corner.y>=0 && top_corner.y<grid_dim-1 )
    j2 = (int) top_corner.y+1;
  else
    if( top_corner.y >= grid_dim-1 )
      j2 = grid_dim-1;
    else
      {is_in = 0; return;}
  if( top_corner.z>=0 && top_corner.z<grid_dim-1 )
    k2 = (int) top_corner.z+1;
  else
    if( top_corner.z >= grid_dim-1 )
      k2 = grid_dim-1;
    else
      {is_in = 0; return;}

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
      ::error ("CLShell constructor: cone number mismatch");
  }
  else
    cone = 0;
}

CLShell::CLShell(const CLShell& old)
{
  coord = old.coord;
  inner_rad = old.inner_rad;
  outer_rad = old.outer_rad;
  outer_rad_sq = old.outer_rad_sq;
  flag = old.flag;
  numcones = old.numcones;

  grid_dim = old.grid_dim;
  is_in = old.is_in;
  i1 = old.i1;
  j1 = old.j1;
  k1 = old.k1;
  i2 = old.i2;
  j2 = old.j2;
  k2 = old.k2;
  if (numcones) {
    cone = new Cone[numcones];
    for (int i=0; i<numcones; ++i)
      cone[i] = old.cone[i];
  }
  else
    cone = 0;
}


CLShell&
CLShell::operator= (const CLShell& old)
{
  coord = old.coord;
  inner_rad = old.inner_rad;
  outer_rad = old.outer_rad;
  outer_rad_sq = old.outer_rad_sq;
  flag = old.flag;
  grid_dim = old.grid_dim;
  is_in = old.is_in;
  i1 = old.i1;
  j1 = old.j1;
  k1 = old.k1;
  i2 = old.i2;
  j2 = old.j2;
  k2 = old.k2;
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

CLShell::~CLShell()
{
  if (numcones && cone) delete [] cone;
}


void CLShell::mark_by_radii(AccTag *acc_array)
{
  if (!is_in) return;
  int nsq = grid_dim*grid_dim;
  float rsq = inner_rad*inner_rad;
  if (flag == free) {
    for(int i=i1; i<=i2; i++ ) {
      int insq = i*nsq;
      float dxsq = (i-coord.x)*(i-coord.x);
      for(int j=j1; j<=j2; j++ ) {
	int jline = j*grid_dim;
	float dysq = (j-coord.y)*(j-coord.y);
	for(int k=k1; k<=k2; k++ ) {
	  int h = insq + jline + k;
	  if( acc_array[h] != interior ) {
	    float dsq = (k-coord.z)*(k-coord.z);
	    dsq = dsq + dxsq + dysq;
	    if( dsq < rsq ) {
	      acc_array[h] = interior;
	    }
	    else {
	      if( dsq < outer_rad_sq) {
		acc_array[h] = exterior;
	      }
	    }
	  }
	}
      }
    }
  }
  else {
    for(int i=i1; i<=i2; i++ ) {
      int insq = i*nsq;
      float dxsq = (i-coord.x)*(i-coord.x);
      for(int j=j1; j<=j2; j++ ) {
	int jline = j*grid_dim;
	float dysq = (j-coord.y)*(j-coord.y);
	for(int k=k1; k<=k2; k++ ) {
	  int h = insq + jline + k;
	  if( acc_array[h] != interior ) {
	    float dsq = (k-coord.z)*(k-coord.z);
	    dsq = dsq + dxsq + dysq;
	    if( dsq < rsq ) {
	      acc_array[h] = interior;
	    }
	    else {
	      if( dsq < outer_rad_sq) {
		acc_array[h] = undecided;
	      }
	    }
	  }
	}
      }
    }
  }
}


void CLShell::mark_by_cones(AccTag *acc_array)
{
  if (!is_in) return;
  if (flag == buried || flag == free) return;  // buried or free atom case
  int nsq = grid_dim*grid_dim;
  Coord cent_to_grid_pt;
  for(int i=i1; i<=i2; i++ ) {
    int insq = i*nsq;
    cent_to_grid_pt.x = i-coord.x;
    float dxsq = cent_to_grid_pt.x * cent_to_grid_pt.x;
    for(int j=j1; j<=j2; j++ ) {
      int jline = j*grid_dim;
      cent_to_grid_pt.y = j-coord.y;
      float dysq = cent_to_grid_pt.y * cent_to_grid_pt.y;
      for(int k=k1; k<=k2; k++ ) {
	int h = insq + jline + k;
	if( acc_array[h] == undecided ) {
	  cent_to_grid_pt.z = k-coord.z;
	  float dsq = cent_to_grid_pt.z * cent_to_grid_pt.z;
	  dsq = dsq + dxsq + dysq;
	  if( dsq <= outer_rad_sq ) {
	    int inside_cone = 0;
	    for( int nc=0; nc<numcones; nc++) {
	      float dot_prod = cent_to_grid_pt * cone[nc].unit_axis;
	      float cone_cos_ang1_sq = cone[nc].cos_ang1_sq;
	      if( cone_cos_ang1_sq >= 0 ) {
		if( dot_prod > 0 ) {
		  float pt_cos_sq = dot_prod*dot_prod / dsq;
		  if (pt_cos_sq > cone_cos_ang1_sq) {
		    inside_cone = 1;
		    break;
		  }
		}
	      }
	      else {  // A negative cone cosine sq means ang > pi/2
		if( dot_prod < 0 ) {
		  float pt_cos_sq = dot_prod*dot_prod / dsq;
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
	    if( inside_cone == 0 ) {
	      acc_array[h] = exterior;
	    }
	  }
	}
      }
    }
  }


}

// CLShell.cc ends here
