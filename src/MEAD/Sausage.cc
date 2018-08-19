/* Sausage is the volume swept by probe rolling between spheres.
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

$Id: Sausage.cc,v 2.6 2007/05/28 01:26:43 bashford Exp $
*/


#include "MEAD/Sausage.h"
#include "MEAD/CubeLatSpec.h"
#include "MEAD/AccTag.h"
#include "MEAD/VertElem.h"
#include "MEAD/Vertex.h"
#include "MEAD/Pair.h"
#include "MEAD/Shell.h"

#include <string>
#include <iostream>
#include <math.h>

// FIXME.  This code is hard to understand.  Maybe Donut and Sausage
// should be subclasses of a Tube base class.  Also beter variable names
// should be used.

// This ctor makes a normal sausage in which the vertices become ends.
Sausage::Sausage (const Pair* pair, const VertElem *v1, const VertElem *v2,
		  float probe_radius)
{
  coord_1 = v1->vpt->coord;
  coord_2 = v2->vpt->coord;
  coord_3 = pair->get_overlap_center();
  angle = v2->angle - v1->angle;
  if (angle < 0)
    ::error("ERROR: Sausage constructor got negative angle from VertElems");
// FIXME! why not do this from pair->get_overlap_center() ?
  rad = (coord_1 - coord_3) * (coord_1 - coord_3);
  Rprob = probe_radius;
}

// This ctor makes a donut.
Sausage::Sausage (const Pair* pair, const Shell *s1, const Shell *s2,
		  float probe_radius)
{
  coord_1 = s1->get_coord();
  coord_2 = s2->get_coord();
  coord_3 = pair->get_overlap_center();
// FIXME rad is really the square of a radius.  Bad name.
  rad = pair->get_overlap_rad_sq();
  angle = -10; // A silly code for donut.
  Rprob = probe_radius;
}

// FIXME.
// pt_inside is ineffecient and needs a substantial re-working of the
// Sausage class.  Sausage and Donut should be subclasses of Tube.
// Many quantities should be calculated by the constructor including,
// axis unit vector, ring radius, etc.  SqAngle class should be used
// as a replacement for the crazy manipulations of angles here and
// it should be possible to avoid any calls to trig functions.
// Similar comments apply to mark_cubelat, below.

// See if a point is inside the tube.
int
Sausage::pt_inside(const Coord& r) const
{
  // FIXME!  Inefficient to calculate these here.
  // Should be Sausage data elements?
  float ring_radius = sqrt(rad);
  float ring_minus_probe = ring_radius - Rprob;
  float ring_plus_probe = ring_radius + Rprob;
  float ring_radius_sq = rad;
  float ring_minus_probe_sq = ring_minus_probe*ring_minus_probe;
  float ring_plus_probe_sq = ring_plus_probe*ring_plus_probe;
  float Rprobsq = Rprob*Rprob;

  Coord cent_to_r = r - coord_3;
  float dis_to_cent_sq = cent_to_r*cent_to_r;

  int inside_saus = 0;

  if (dis_to_cent_sq > ring_plus_probe_sq) return 0;

  if( angle == -10.0 ) {                                     //   Donut case
    if(dis_to_cent_sq>ring_minus_probe_sq
       && dis_to_cent_sq<ring_plus_probe_sq) {
      Coord donut_axis = coord_2 - coord_1;
      float inv_axis_lensq = 1/(donut_axis*donut_axis);
      float axis_dot_veccent = donut_axis*cent_to_r;
      float dis_to_plan_sq = inv_axis_lensq*axis_dot_veccent*axis_dot_veccent;
      if(dis_to_plan_sq <= Rprobsq) {
	float sin_ang_sq = dis_to_plan_sq/dis_to_cent_sq;
	float dis_to_arc_sq;
	//should have two distances, only short one is needed in donut cases
	if( sin_ang_sq <= 1 )
	  dis_to_arc_sq = dis_to_cent_sq + ring_radius_sq
	    - 2*sqrt(dis_to_cent_sq*ring_radius_sq*(1-sin_ang_sq));
	else {
	  if( sin_ang_sq-1 < 0.01 )
	    dis_to_arc_sq = dis_to_cent_sq + ring_radius_sq;
	  else {
	    dis_to_arc_sq = 0;  // prevent compiler warning
	    cerr << "Sausage::pt_inside: Angular domain error.  " <<
		      " Sin^2(x) =  " << sin_ang_sq <<  endl;
	    ::error();
	  }
	}
	if (dis_to_arc_sq < Rprobsq)
	  inside_saus = 1;
      }
    }
  }
  else {  // Sausage case
    if( dis_to_cent_sq>ring_minus_probe_sq
       && dis_to_cent_sq<ring_plus_probe_sq ) {
      // Substitute some sane names for crazy member varible names (FIXME!)
      const Coord& center = coord_3;
      const Coord& head = coord_1;
      const Coord& tail = coord_2;
      // Construct unit vector perpendicular to the plane
      // FIXME!! It is slow to do this here.  Should be data member.
      Coord uperp = cross(head-center, tail-center);
      uperp = uperp / sqrt(uperp*uperp);
      float dis_to_plan = uperp * cent_to_r;
      float dis_to_plan_sq = dis_to_plan*dis_to_plan;
      if( dis_to_plan_sq <= Rprobsq) {
	// See if r is inside sphere sitting at either head or tail.
	Coord head_to_r = r - head;
	float dis_to_sg_head_sq = head_to_r * head_to_r;
	if( dis_to_sg_head_sq < Rprobsq )
	  inside_saus = 1;
	else {
	  Coord tail_to_r = r - tail;
	  float dis_to_sg_tail_sq = tail_to_r * tail_to_r;
	  if( dis_to_sg_tail_sq < Rprobsq)
	    inside_saus = 1;
	  else { // End spheres did not decide it so calc dist to arc, etc.
	    float dis_to_arc1_sq,dis_to_arc2_sq;
	    float sin_ang_sq = dis_to_plan_sq/dis_to_cent_sq;
	    if( sin_ang_sq < 1 ) {
	      float d = 2*sqrt(dis_to_cent_sq*ring_radius_sq*(1-sin_ang_sq));
	      dis_to_arc1_sq = dis_to_cent_sq + ring_radius_sq - d;
	      dis_to_arc2_sq = dis_to_cent_sq + ring_radius_sq + d;
	    }
	    else {
	      if( sin_ang_sq-1 < 0.001 )
		dis_to_arc1_sq=dis_to_arc2_sq=dis_to_cent_sq+ring_radius_sq;
	      else {
		// line below just to prevent compiler warnings
		dis_to_arc1_sq = dis_to_arc2_sq = 0;
		cerr << "Angular domain error in Sausage.cc : "<<
		  " Sin^2(x) =  " << sin_ang_sq << endl;
		::error();
	      }
	    }
	    // Check each dis_to_arc to see if it is valid and inside sausage.
	    if( dis_to_arc1_sq <= Rprobsq) {
	      //find vector from sg_cent to the proj of on sg_plane
	      Coord cen_to_proj = cent_to_r - uperp*dis_to_plan;
	      float denomin = sqrt((cen_to_proj*cen_to_proj)*ring_radius_sq);
	      float dot_prod = (head-center)*cen_to_proj;
	      float ang1;
	      if(fabs( dot_prod ) <= denomin)
		ang1 = acos( dot_prod/denomin );
	      else {
		if( fabs( dot_prod )-denomin < 0.01 )
		  ang1 = dot_prod > 0.0F ? 0.0F : 3.1415926F;
		else {
		  ang1 = 0; // prevent compiler warning
		  cout<<"Angular domain error in Sausage.cc"<<endl;
		  cout<<"Cos(x) > 1  "<< endl;
		}
	      }
	      dot_prod = (tail-center)*cen_to_proj;
	      float ang2;
	      if( fabs( dot_prod ) < denomin  )
		ang2 = acos( dot_prod/denomin );
	      else {
		if( fabs( dot_prod )-denomin < 0.01 )
		  ang2 = dot_prod > 0.0F ? 0.0F : 3.1415926F;
		else {
		  ang2 = 0; // prevent compiler warning
		  cout<<"Angular domain error in Sausage.cc"<<endl;
		  cout<<"Cos(x) > 1 "<<endl;
		}
	      }
	      if( fabs( ang1+ang2-angle ) < 0.01 )
		inside_saus = 1;
	      else {
		if( fabs(ang1+ang2+angle-6.283185) < 0.01 ) {
		  if(dis_to_arc2_sq < Rprobsq)
		    inside_saus = 1;
		}
		else
		  if( angle > 3.1415926 ) {
		    float dang = fabs( ang1 - ang2 );
		    if( fabs( dang-(6.283185-angle) ) < 0.01 )
		      inside_saus = 1;
		  }
	      }
	    }
	  }
	}
      }
    }
  }
  return inside_saus;
}

void
Sausage::mark_cubelat(const CubeLatSpec& cls, AccTag acc_array[])
{
  // The Sausage internal data is in real space units, so we need
  // local variables in grid space units.

  float spacing = cls.get_spacing();
  float sp_coeff = 1/spacing;
  Coord grid_center_in_space = cls.get_center();
  int grid_dim = cls.get_grid_dim();
  int nsq = grid_dim*grid_dim;
  float grlen = (float) (grid_dim - 1);
  float halfgrlen = grlen/2;
  Coord grid_center_in_grid (halfgrlen, halfgrlen, halfgrlen);
  float sp_coeff_sq = sp_coeff*sp_coeff;
  float Rprob_in_grid = Rprob*sp_coeff;
  float Rprobsq_in_grid = Rprob_in_grid*Rprob_in_grid;

  Coord sg_center_in_grid = (coord_3 - grid_center_in_space)
    *sp_coeff + grid_center_in_grid;
  Coord sg_head_in_grid = (coord_1 - grid_center_in_space)
    *sp_coeff + grid_center_in_grid;
  Coord sg_tail_in_grid = (coord_2 - grid_center_in_space)
    *sp_coeff + grid_center_in_grid;

  float r = sqrt( rad );
  float r0 = r - Rprob;
  float r1 = r + Rprob;
  float r1_in_grid = r1*sp_coeff;
  float rsq_in_grid = rad*sp_coeff_sq;
  float r0sq_in_grid = r0*r0*sp_coeff_sq;
  float r1sq_in_grid = r1*r1*sp_coeff_sq;

  // Now define corners of sausage_covering box to be scanned
  // in grid-base coords.

  Coord top_corner, bottom_corner;

  if( angle < 3.1415926 && angle != -10.0 ) {
    if( sg_center_in_grid.x>sg_tail_in_grid.x
       && sg_center_in_grid.x>sg_head_in_grid.x ) {
      bottom_corner.x = sg_center_in_grid.x - r1_in_grid;
      if( sg_tail_in_grid.x > sg_head_in_grid.x )
	top_corner.x = sg_tail_in_grid.x + Rprob_in_grid;
      else
	top_corner.x = sg_head_in_grid.x + Rprob_in_grid;
    }
    else
      if( sg_center_in_grid.x<sg_tail_in_grid.x
	 && sg_center_in_grid.x<sg_head_in_grid.x ) {
	top_corner.x = sg_center_in_grid.x + r1_in_grid;
	if( sg_tail_in_grid.x < sg_head_in_grid.x )
	  bottom_corner.x = sg_tail_in_grid.x - Rprob_in_grid;
	else
	  bottom_corner.x = sg_head_in_grid.x - Rprob_in_grid;
      }
      else {
	if( sg_tail_in_grid.x < sg_head_in_grid.x ) {
	  top_corner.x = sg_head_in_grid.x + Rprob_in_grid;
	  bottom_corner.x = sg_tail_in_grid.x - Rprob_in_grid;
	}
	else {
	  top_corner.x = sg_tail_in_grid.x + Rprob_in_grid;
	  bottom_corner.x = sg_head_in_grid.x - Rprob_in_grid;
	}
      }
    if( sg_center_in_grid.y>sg_tail_in_grid.y
       && sg_center_in_grid.y>sg_head_in_grid.y ) {
      bottom_corner.y = sg_center_in_grid.y - r1_in_grid;
      if( sg_tail_in_grid.y > sg_head_in_grid.y )
	top_corner.y = sg_tail_in_grid.y + Rprob_in_grid;
      else
	top_corner.y = sg_head_in_grid.y + Rprob_in_grid;
    }
    else
      if( sg_center_in_grid.y<sg_tail_in_grid.y
	 && sg_center_in_grid.y<sg_head_in_grid.y ) {
	top_corner.y = sg_center_in_grid.y + r1_in_grid;
	if( sg_tail_in_grid.y < sg_head_in_grid.y )
	  bottom_corner.y = sg_tail_in_grid.y - Rprob_in_grid;
	else
	  bottom_corner.y = sg_head_in_grid.y - Rprob_in_grid;
      }
      else {
	if( sg_tail_in_grid.y < sg_head_in_grid.y ) {
	  top_corner.y = sg_head_in_grid.y + Rprob_in_grid;
	  bottom_corner.y = sg_tail_in_grid.y - Rprob_in_grid;
	}
	else {
	  top_corner.y = sg_tail_in_grid.y + Rprob_in_grid;
	  bottom_corner.y = sg_head_in_grid.y - Rprob_in_grid;
	}
      }
    if( sg_center_in_grid.z>sg_tail_in_grid.z
       && sg_center_in_grid.z>sg_head_in_grid.z ) {
      bottom_corner.z = sg_center_in_grid.z - r1_in_grid;
      if( sg_tail_in_grid.z > sg_head_in_grid.z )
	top_corner.z = sg_tail_in_grid.z + Rprob_in_grid;
      else
	top_corner.z = sg_head_in_grid.z + Rprob_in_grid;
    }
    else
      if( sg_center_in_grid.z<sg_tail_in_grid.z
	 && sg_center_in_grid.z<sg_head_in_grid.z ) {
	top_corner.z = sg_center_in_grid.z + r1_in_grid;
	if( sg_tail_in_grid.z < sg_head_in_grid.z )
	  bottom_corner.z = sg_tail_in_grid.z - Rprob_in_grid;
	else
	  bottom_corner.z = sg_head_in_grid.z - Rprob_in_grid;
      }
      else {
	if( sg_tail_in_grid.z < sg_head_in_grid.z ) {
	  top_corner.z = sg_head_in_grid.z + Rprob_in_grid;
	  bottom_corner.z = sg_tail_in_grid.z - Rprob_in_grid;
	}
	else {
	  top_corner.z = sg_tail_in_grid.z + Rprob_in_grid;
	  bottom_corner.z = sg_head_in_grid.z - Rprob_in_grid;
	}
      }
  }
  else {
    Coord diag(r1_in_grid, r1_in_grid, r1_in_grid);
    bottom_corner = sg_center_in_grid - diag;
    top_corner = sg_center_in_grid + diag;
  }

  int i1, i2, j1, j2, k1, k2;

  if( bottom_corner.x>=0 && bottom_corner.x<grid_dim)
    i1 = (int) bottom_corner.x;
  else
    if( bottom_corner.x < 0 )
      i1 = 0;
    else
      return;
  if( bottom_corner.y>=0 && bottom_corner.y<grid_dim)
    j1 = (int) bottom_corner.y;
  else
    if( bottom_corner.y < 0 )
      j1 = 0;
    else
      return;
  if( bottom_corner.z>=0 && bottom_corner.z<grid_dim)
    k1 = (int) bottom_corner.z;
  else
    if( bottom_corner.z < 0 )
      k1 = 0;
    else
      return;
  if( top_corner.x>=0 && top_corner.x<grid_dim-1 )
    i2 = (int) top_corner.x+1;
  else
    if( top_corner.x >= grid_dim-1 )
      i2 = grid_dim-1;
    else
      return;
  if( top_corner.y>=0 && top_corner.y<grid_dim-1 )
    j2 = (int) top_corner.y+1;
  else
    if( top_corner.y >= grid_dim-1 )
      j2 = grid_dim-1;
    else
      return;
  if( top_corner.z>=0 && top_corner.z<grid_dim-1 )
    k2 = (int) top_corner.z+1;
  else
    if( top_corner.z >= grid_dim-1 )
      k2 = grid_dim-1;
    else
      return;

  // check all nodes inside the box

  if( angle == -10.0 ) {                                     //   Donut case
    float coeff = 1/( ( sg_tail_in_grid-sg_head_in_grid )
		     * ( sg_tail_in_grid-sg_head_in_grid ) );
    for(int i=i1; i<=i2; i++ ) {
      int insq = i*nsq;
      float a = (sg_tail_in_grid.x-sg_head_in_grid.x)*(i-sg_center_in_grid.x);
      float dxsq = (i-sg_center_in_grid.x)*(i-sg_center_in_grid.x);
      for(int j=j1; j<=j2; j++ ) {
	int jline = j*grid_dim;
	float b = (sg_tail_in_grid.y-sg_head_in_grid.y)
	  *(j-sg_center_in_grid.y);
	float dysq = (j-sg_center_in_grid.y)*(j-sg_center_in_grid.y);
	for(int k=k1; k<=k2; k++ )  {
	  int h = insq + jline + k;
	  if( acc_array[h] == undecided ) {
	    Coord grid_pt_in_grid(static_cast<float>(i),
			static_cast<float>(j),
			static_cast<float>(k));
	    float dis_to_cent_sq = (k-sg_center_in_grid.z)*(k-sg_center_in_grid.z);
	    dis_to_cent_sq = dis_to_cent_sq + dxsq + dysq;
	    if( dis_to_cent_sq>r0sq_in_grid && dis_to_cent_sq<r1sq_in_grid ) {
	      float c = (sg_tail_in_grid.z-sg_head_in_grid.z)
		* (k-sg_center_in_grid.z);
	      float dis_to_plan_sq = coeff*(a+b+c)*(a+b+c);
	      if( dis_to_plan_sq <= Rprobsq_in_grid ) {
		float sin_ang_sq = dis_to_plan_sq/dis_to_cent_sq;
                float dis_to_arc;  //should have two distances, only the short
		if( sin_ang_sq <= 1 )  //one is needed in donut cases
		  dis_to_arc = dis_to_cent_sq + rsq_in_grid
		    - 2*sqrt(dis_to_cent_sq*rsq_in_grid*(1-sin_ang_sq));
		else {
		  if( sin_ang_sq-1 < 0.01 )
		    dis_to_arc = dis_to_cent_sq + rsq_in_grid;
		  else {
		    dis_to_arc = 0; // prevent compiler warnings
		    cout << "Angular domain error in Sausage.cc : " <<
		      " Sin(x) =  " << sqrt(sin_ang_sq) <<  endl;
		    ::error();
		  }
		}
		if( dis_to_arc < Rprobsq_in_grid )
		  acc_array[h] = in_tube;
	      }
	    }
	  }
	}
      }
    }
  }
  else {                                                  // Real sausage cases
    Coord v1 = sg_head_in_grid - sg_center_in_grid;
    Coord v2 = sg_tail_in_grid - sg_center_in_grid;
    Coord v = cross(v1,v2);
    float coeff = 1/sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
    float a = v.x * coeff;
    float b = v.y * coeff;  //a,b,c determine the normal of sausage plane
    float c = v.z * coeff;

    for(int i=i1; i<=i2; i++ ) {
      int insq = i*nsq;
      for(int j=j1; j<=j2; j++ ) {
	int jline = j*grid_dim;
	for(int k=k1; k<=k2; k++ ) {
	  int h = insq + jline + k;
	  if( acc_array[h] == undecided ) {
	    Coord grid_pt_in_grid(static_cast<float>(i),
			static_cast<float>(j),
			static_cast<float>(k));
	    float dis_to_cent_sq = (grid_pt_in_grid-sg_center_in_grid)
	      * (grid_pt_in_grid-sg_center_in_grid);
	    if( dis_to_cent_sq>r0sq_in_grid && dis_to_cent_sq<r1sq_in_grid ) {
	      float dis_to_plan = a*(grid_pt_in_grid.x-sg_center_in_grid.x)
		                + b*(grid_pt_in_grid.y-sg_center_in_grid.y)
		                + c*(grid_pt_in_grid.z-sg_center_in_grid.z);
	      float dis_to_plan_sq = dis_to_plan*dis_to_plan;
	      if( dis_to_plan_sq <= Rprobsq_in_grid ) {
		float inside_saus = 0;
		float dis_to_sg_head_sq = ( grid_pt_in_grid - sg_head_in_grid )
		  * ( grid_pt_in_grid - sg_head_in_grid );
		if( dis_to_sg_head_sq < Rprobsq_in_grid )
		  inside_saus = 1;
		else {
		  float dis_to_sg_tail_sq = ( grid_pt_in_grid - sg_tail_in_grid )
		    * ( grid_pt_in_grid - sg_tail_in_grid );
		  if( dis_to_sg_tail_sq < Rprobsq_in_grid )
		    inside_saus = 1;
		  else {
		    float dis_to_arc1_sq,dis_to_arc2_sq;
		    float sin_ang_sq = dis_to_plan_sq/dis_to_cent_sq;
		    if( sin_ang_sq < 1 ) {
		      float d = 2*sqrt(dis_to_cent_sq*rsq_in_grid*(1-sin_ang_sq));
		      dis_to_arc1_sq = dis_to_cent_sq + rsq_in_grid - d;
		      dis_to_arc2_sq = dis_to_cent_sq + rsq_in_grid + d;
		    }
		    else {
		      if( sin_ang_sq-1 < 0.001 )
			dis_to_arc1_sq=dis_to_arc2_sq=dis_to_cent_sq+rsq_in_grid;
		      else {
			// line below just to prevent compiler warnings
			dis_to_arc1_sq = dis_to_arc2_sq = 0;
			cerr << "Angular domain error in Sausage.cc : "<<
			  " Sin(x) =  " << sqrt(sin_ang_sq) << endl;
			::error();
		      }
		    }

                //check each dis_to_arc to see if it is valid and inside sausage

		    if( dis_to_arc1_sq <= Rprobsq_in_grid ) {
                      //find vector from sg_cent to the proj of (ijk) on sg_plane
		      float x1 = (i-sg_center_in_grid.x) - a*dis_to_plan;
		      float y1 = (j-sg_center_in_grid.y) - b*dis_to_plan;
		      float z1 = (k-sg_center_in_grid.z) - c*dis_to_plan;
		      float denomin = sqrt( (x1*x1+y1*y1+z1*z1)*rsq_in_grid );
		      float dot_prod = x1*(sg_head_in_grid.x-sg_center_in_grid.x)
		                     + y1*(sg_head_in_grid.y-sg_center_in_grid.y)
		      	             + z1*(sg_head_in_grid.z-sg_center_in_grid.z);
		      float ang1;
		      if( fabs( dot_prod ) <= denomin  )
			ang1 = acos( dot_prod/denomin );
		      else {
			if( fabs( dot_prod )-denomin < 0.01 )
			  ang1 = dot_prod > 0.0F ? 0.0F : 3.1415926F;
			else {
			  ang1 = 0; // prevent compler warnings
			  cout<<"Angular domain error in Sausage.cc"<<endl;
			  cout<<"Cos(x) > 1  "<< endl;
			}
		      }
		      dot_prod = x1*(sg_tail_in_grid.x-sg_center_in_grid.x)
			       + y1*(sg_tail_in_grid.y-sg_center_in_grid.y)
			       + z1*(sg_tail_in_grid.z-sg_center_in_grid.z);
		      float ang2;
		      if( fabs( dot_prod ) < denomin  )
			ang2 = acos( dot_prod/denomin );
		      else {
			if( fabs( dot_prod )-denomin < 0.01 )
			  ang2 = dot_prod > 0.0F ? 0.0F : 3.1415926F;
			else {
			  ang2 = 0; // prevent compiler warnings
			  cout<<"Angular domain error in Sausage.cc"<<endl;
			  cout<<"Cos(x) > 1 "<<endl;
			}
		      }
		      if( fabs( ang1+ang2-angle ) < 0.01 )
			inside_saus = 1;

		      else {
			if( fabs(ang1+ang2+angle-6.283185) < 0.01 ) {
			  if( dis_to_arc2_sq < Rprobsq_in_grid )
			    inside_saus = 1;
			}
			else
			  if( angle > 3.1415926 ) {
			    float dang = fabs( ang1 - ang2 );
			    if( fabs( dang-(6.283185-angle) ) < 0.01 )
			      inside_saus = 1;
			  }
		      }
		    }
		  }
		}
		if( inside_saus )
		  acc_array[h] = in_tube;
	      }
	    }
	  }
	}
      }
    }
  }
}

ostream& Sausage::write_top_in_binary(ostream& bfile)
{
 float top_data[12];
 if( !bfile )
   {
    cerr<<"Sausage::write_top_in_binary : can not open the file"<<endl;
    return bfile;
   }
 else
   {
    top_data[0] = rad;
    top_data[1] = angle;
    top_data[2] = Rprob;
    top_data[3] = coord_1.x;
    top_data[4] = coord_1.y;
    top_data[5] = coord_1.z;
    top_data[6] = coord_2.x;
    top_data[7] = coord_2.y;
    top_data[8] = coord_2.z;
    top_data[9] = coord_3.x;
    top_data[10] = coord_3.y;
    top_data[11] = coord_3.z;
    bfile.write((char *) top_data, sizeof(float)*12);
   }
 return bfile;
}

ostream& Sausage::write_top_in_ascii(ostream& afile)
{
 afile << rad << "  "  << angle << "  " << Rprob << endl;
 coord_1.print(afile);
 afile << "  ";
 coord_2.print(afile);
 afile << "  ";
 coord_3.print(afile);
 return afile;
}

int Sausage::read_top_in_binary(istream& bfile)
{
 float top_data[12];
 bfile.read((char *) top_data, sizeof(float)*12);
 if( !bfile )
   {
    cerr<<"ERROR : Sausage::read_top_in_binary "<<endl;
    return 0;
   }
 rad = top_data[0];
 angle = top_data[1];
 Rprob = top_data[2];
 coord_1.x = top_data[3];
 coord_1.y = top_data[4];
 coord_1.z = top_data[5];
 coord_2.x = top_data[6];
 coord_2.y = top_data[7];
 coord_2.z = top_data[8];
 coord_3.x = top_data[9];
 coord_3.y = top_data[10];
 coord_3.z = top_data[11];
 return 1;
}

int Sausage::read_top_in_ascii(istream& afile)
{
 string word;
 int sausage_numb;
 afile >> word;
 afile >> sausage_numb;
 afile >> rad;
 afile >> angle;
 afile >> Rprob;
 coord_1.read(afile);
 coord_2.read(afile);
 coord_3.read(afile);
 return 1;
}

// Sausage.cc ends here
