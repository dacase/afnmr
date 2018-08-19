/* calculate analytical rep. of solv. accessible vol from spheres and radii.
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
$Id: SAVanal_calc.cc,v 2.12 2012/04/09 20:15:39 agiammon Exp $
*/

#include <vector>
using std::vector;

#include "MEAD/SolvAccVol.h"
#include "MEAD/Shell.h"
#include "MEAD/Pair.h"
#include "MEAD/Cone.h"
#include "MEAD/Sausage.h"
#include "MEAD/globals.h"

#include <iostream>
#include <math.h>


void SolvAccVolRep::anal_calc()
{
  if (is_calculated) {
    cerr << "SolvAccVolRep::anal_calc: WARNING: SolvAccVolRep, "
      << ",\nhas already been calculated\n"
	"returning with no new caluclation" << endl;
    return;
  }

  blab1 << "SolvAccVolRep::anal_calc doing calculation" << endl;



  /*  Find all completely buried atoms and exclude them from atom list */
  int atcount = sph_count;
  sph_count = 0;
  int i;
  for(i=0; i<atcount; i++ )
    {
      if( atsph[i].is_free() )
	{
	  for(int j=i+1; j<atcount; j++ )
	    {
	      float diffradsq = atsph[i].get_inner_rad()
				- atsph[j].get_inner_rad();
	      diffradsq *= diffradsq;
	      Coord rij = atsph[j].get_coord() - atsph[i].get_coord();
	      float disq = rij*rij;
	      if (disq <= diffradsq)
		{
                  if( atsph[i].get_inner_rad() > atsph[j].get_inner_rad() )
		    atsph[j].bury();
                  else
                    {
                      atsph[i].bury();
                      break;
                    }
		}
	    }
	  if( atsph[i].is_free() )
	    sph_count++;
	}
    }

  if( sph_count <= 0 )
    ::error("SolvAccVolRep::anal_calc went from non-zero atoms to spheres <= 0!");

// Transfer the non-buried atoms to a new list just the right size

  Shell *atsph1 = new Shell [sph_count];
  int olsp=0, newsp=0;
  for(; olsp<atcount; ++olsp ) {
    if( atsph[olsp].is_free() ) {
      atsph1[newsp] = atsph[olsp];
      ++newsp;
    }
  }
  if (newsp != sph_count)
    ::error("SolvAccVolRep::anal_calc: sphere count went wrong");
  delete [] atsph;
  atsph = atsph1;
  blab3 << "Finishing checking all buried atoms\t# 0f sphere = "
    << newsp <<endl;

// Build the pair data structures

  Pair **pair_array = new Pair* [sph_count];
  int *pairlist_length_at = new int [sph_count];
  Pair *tmp_pair = new Pair [sph_count];
  int total_number_of_pairs = 0;

  for( i=0; i<sph_count; i++ )
    {
// Collect elements for the i-slice of the pair array in tmp_pair.
// (Also do some flagging of spheres).
      Pair *ij_pair = tmp_pair;
      int i_pair_count = 0;
      for(int j=i+1; j<sph_count; j++ )
	{
	  float sumradsq = atsph[i].get_outer_rad()+atsph[j].get_outer_rad();
	  sumradsq *= sumradsq;
	  Coord rij = atsph[j].get_coord() - atsph[i].get_coord();
	  float disq = rij*rij;
	  if (disq==0.0) {
	    // Spheres coincide.  They must have identical radii or
	    // else the smaller should have been removed above.
	    if (atsph[i].get_outer_rad() != atsph[j].get_outer_rad())
	      ::error("SAV::anal_calc: unexpected coinciding spheres",
		      " of different radii\n");
	    cerr << "SAV::anal_calc: WARNING: unusual case of spheres\n"
		 << "  with coinciding centers and identical radii.\n"
		 << "  Geometric algorithms could get shaky here!" << endl;
	    // For such sphere pairs, don't make a Pair or mark anything.
	    // Would this "missing" pair cause triplet problems?  FIXME?
	  }
	  if (disq < sumradsq)
	    {
	      Pair p(i, atsph[i], j, atsph[j]);
	      *ij_pair = p;
	      ++ij_pair;
	      ++i_pair_count;
	      ++total_number_of_pairs;

	      /* After this checking all atom should be flagged
		 buried.  These are initial flag values.  Later some
		 of these will become partially_buried. */
	      atsph[i].bury();
	      atsph[j].bury();
	    }
	}
      // Allocate an i-slice of the array and fill it from tmp_array
      pairlist_length_at[i] = i_pair_count;
      pair_array[i] = new Pair [i_pair_count];
      for (int jj=0; jj<i_pair_count; ++jj)
	pair_array[i][jj] = tmp_pair[jj];  // Pointers might be faster here.

      if( atsph[i].is_free() )
	blab3 << "atom " << i << " is a isolated atom" << endl;
    }
  delete [] tmp_pair;
  blab2 << "Number of pairs is "<< total_number_of_pairs << endl;

  /* Find and validate all geometrically allowed vertexes to identify
      those on surface */

#ifdef POPULATION_COUNT
  cout << "Before the triplet loop:" << endl;
  Vertex::population_report();
  VertElem::population_report();
#endif
  int at_triplets = 0;
  int sp_triplets = 0;

  int elim_vertex = 0;

  for(int at_i=0; at_i<sph_count-1; at_i++ )
    {
      Pair *ij_pair = pair_array[at_i];
      Pair *ij_pair_end = ij_pair + pairlist_length_at[at_i];
      for( ; ij_pair<ij_pair_end; ij_pair++ )
	{
	  int at_j = ij_pair->j_index();
	  Coord rij = atsph[at_j].get_coord() - atsph[at_i].get_coord();;
          float rij_sq = rij*rij;
	  float Ri_sq = atsph[at_i].get_outer_rad_sq();
	  float Rj_sq = atsph[at_j].get_outer_rad_sq();
	  float rij_dot_Ri = ( Ri_sq + rij_sq - Rj_sq )*0.5F;
	  for(Pair *ik_pair=ij_pair+1; ik_pair<ij_pair_end; ik_pair++ )
	    {
	      int at_k = ik_pair->j_index();
	      // We now have at_i, at_j and at_k.  Are they a triplex?
	      Pair *jk_pair = pair_array[at_j];
	      Pair *jk_pair_end = jk_pair + pairlist_length_at[at_j];

	      for( ; jk_pair<jk_pair_end && jk_pair->j_index() < at_k;
		  jk_pair++ )
		;
	      if (jk_pair >= jk_pair_end || jk_pair->j_index() != at_k)
		  continue;	// This is not a triplex.

	      // Otherwise, this IS a triplex, so process it...

	      at_triplets++;

	      /* this if statement eliminates the calculations and
		 validation of vertexes involving three buried
		 pairs. However, it depends on how fast all buried
		 pairs being flagged, and the number of vertexes
		 eliminated may be a very small portion of the total.
		 so it may not lead to a speed up but an opposite
		 effect (as happened in tRNA case) */

// FIXME this condition removed for now so that all tripets will be
// considered.  This maximizes pair burial
//	      if( jk_pair->flag!=Pair::buried || ik_pair->flag!=Pair::buried
//		 || ij_pair->flag!=Pair::buried )
	      if(1)
		{
		  float alpha, beta, gamma;
		  Coord rik = atsph[at_k].get_coord() - atsph[at_i].get_coord();
		  float rik_sq = rik*rik;
		  float Rk_sq = atsph[at_k].get_outer_rad_sq();
		  float rij_dot_rik = rij*rik;
		  float rik_dot_Ri = ( Ri_sq + rik_sq - Rk_sq )*0.5F;
		  float denom =  rij_sq*rik_sq - rij_dot_rik*rij_dot_rik;
		  if( denom < 0.00005 )
		    {
		      blab1 << "Colinear case for the atom triplet" << endl;
		      gamma = -1.0;
		      alpha = beta = 0;
		    }
		  else
		    {
		      denom = 1/denom;
		      alpha = (rik_sq*rij_dot_Ri - rij_dot_rik*rik_dot_Ri)*denom;
		      beta = (rij_sq*rik_dot_Ri - rij_dot_rik*rij_dot_Ri)*denom;
		      gamma = Ri_sq - alpha*alpha*rij_sq-beta*beta*rik_sq
			    - 2*alpha*beta*rij_dot_rik;
		    }

		  if( gamma>0 ) // They form two geometrically possible vertices
		    {
		      sp_triplets++;
		      Coord rij_cross_rik = cross(rij,rik);
		      float rij_cross_rik_sq = rij_cross_rik*rij_cross_rik;
		      gamma = sqrt(gamma/rij_cross_rik_sq);
		      Coord temp1 = atsph[at_i].get_coord()+alpha*rij+beta*rik;
		      Coord temp2 = gamma * rij_cross_rik;
		      Coord trip1 = temp1 - temp2;
		      Coord trip2 = temp1 + temp2;

		      VertElem vt1( trip1 );
		      VertElem vt2( trip2 );

		      //validate the vertexes

		      if(! ij_pair->is_fully_buried())
			ij_pair->Check_vertex_pair(vt1,vt2);

		      if(! ik_pair->is_fully_buried())
			ik_pair->Check_vertex_pair(vt2,vt1);

		      if(! jk_pair->is_fully_buried())
			jk_pair->Check_vertex_pair(vt1,vt2);
		    }
		  else
		    /* Check if any of the overlapping circles are inside
		       another atom.
		       */
		    {
		      ij_pair->quickcheck_burial_by_shell(atsph[at_k]);
		      ik_pair->quickcheck_burial_by_shell(atsph[at_j]);
		      jk_pair->quickcheck_burial_by_shell(atsph[at_i]);
		    }
		}
              else
                elim_vertex++;
	    }
	}
      if( (at_i%100) == 0 )
	blab3 << "Finishing atom " << at_i << endl;
    }
#ifdef POPULATION_COUNT
  cout << "After the triplet loop:" << endl;
  Vertex::population_report();
  VertElem::population_report();
#endif
  blab2 << "Number of atom triplets is " << at_triplets << endl;
  blab2 << "Number of geometrically allowed triplets is " << sp_triplets << endl;
  blab2 << "Number of eliminated triplets is " << elim_vertex << endl;


  int nburied_pairs = 0;
  int nburied_shells = 0;
  int ngeom_buried_pairs = 0;
  // Take a trip through the pair list to correctly flag spheres.
  for( i=0; i<sph_count; i++ )
    {
      Shell *at_i_pt = atsph + i;
      Pair *ij_pair = pair_array[i];
      Pair *ij_pair_end = ij_pair + pairlist_length_at[i];
      for( ; ij_pair<ij_pair_end; ij_pair++ )
	{
	  Shell *at_j_pt = atsph + ij_pair->j_index();
	  switch (ij_pair->flag)
	    {
	    case  Pair::buried:
	      ++nburied_pairs;
	      break;
	    case Pair::geometrically_buried:
	      ++ngeom_buried_pairs;
	      break;
	    default:
	      at_i_pt->partially_bury();
	      at_j_pt->partially_bury();
	      break;
	    }
	}
      if (at_i_pt->is_buried()) ++nburied_shells;
    }

  blab3 << "Number of geometrically geom_buried_pairs = "<<ngeom_buried_pairs
    << "\nNumber of otherwise buried pairs = " << nburied_pairs
      << "\nNumber of buried shells = " << nburied_shells
	<< endl;

  int numb = 0;
  int totcones = 0;
  int numb_donut = 0;
  vector<Sausage> sg_plex;
  for( i=0; i<sph_count; i++ )
    {
      Shell *at_i_pt = atsph + i;
      float rad_i_sq = at_i_pt->get_outer_rad_sq();
      Coord coord_i = at_i_pt->get_coord();
      int i_not_buried = ! at_i_pt->is_buried();
      Pair *ij_pair = pair_array[i];
      Pair *ij_pair_end = ij_pair + pairlist_length_at[i];
      for( ; ij_pair<ij_pair_end; ij_pair++ )
	{
	  Shell *at_j_pt = atsph + ij_pair->j_index();
	  Coord coord_j = at_j_pt->get_coord();
	  int j_not_buried = ! at_j_pt->is_buried();

	  if (ij_pair->flag != Pair::geometrically_buried
	      && (i_not_buried || j_not_buried))
	    {
	      // Make the corrsponding cones for this pair and add to spheres
	      Coord vec_ij = coord_j - coord_i;
	      float atom_dis_sq = vec_ij * vec_ij;
	      Coord unit_vec_ij = vec_ij / sqrt (atom_dis_sq);
	      float rad_j_sq = at_j_pt->get_outer_rad_sq();

	      float d1 = 0.25F*(atom_dis_sq + rad_i_sq - rad_j_sq )
		* (atom_dis_sq + rad_i_sq - rad_j_sq) / atom_dis_sq;
	      float arc_rad_sq = rad_i_sq - d1;
	      float d2 = rad_j_sq - arc_rad_sq;

	      if (i_not_buried) { // Make the cone for the i atom.
		Cone cn;
		cn.unit_axis = unit_vec_ij;
                if (d1 > d2)
		    cn.cos_ang1_sq = d1/rad_i_sq;
                else {
                   if(atom_dis_sq > d2)
                      cn.cos_ang1_sq = d1/rad_i_sq;
                   else
                      cn.cos_ang1_sq = -d1/rad_i_sq;
		 }
		at_i_pt->add_cone(cn);
		++totcones;
	      }
	      if (j_not_buried) { // Make the cone for the j atom.
		Cone cn;
		cn.unit_axis = -unit_vec_ij;
                if(d1 < d2)
                   cn.cos_ang1_sq = d2/rad_j_sq;
                else {
		   if (atom_dis_sq > d1)
		       cn.cos_ang1_sq = d2/rad_j_sq;
		   else
		       cn.cos_ang1_sq = -d2/rad_j_sq;
		 }
		at_j_pt->add_cone(cn);
		++totcones;
	      }
	    }
	  // Make the Sausages for this pair, if any.
	  if( ij_pair->flag == Pair::partially_buried )
	    {
	      VertElem *vert_el_ptr = ij_pair->get_vert_head();
	      Coord c0 = ij_pair->get_overlap_center();
	      if( vert_el_ptr == 0 )
		::error("ERROR: anal_calc: A partially_buried pair ",
			"has a null vertex had pointer");
	      while(vert_el_ptr)
		{
		  VertElem* vert1 = vert_el_ptr;
		  if (vert_el_ptr->vpt->count != 3) {
		    blab2 << "WARNING: SAVanal_calc: vertex found with count = "
		      << vert_el_ptr->vpt->count << "\n";
		  }
		  vert_el_ptr = vert_el_ptr->next;
		  if(vert_el_ptr == 0 )
		    ::error("ERROR: anal_calc: Odd number of elements ",
			    "in a vertex list");
		  VertElem* vert2 = vert_el_ptr;
		  if (vert_el_ptr->vpt->count != 3) {
		    blab2 << "WARNING: SAVanal_calc: vertex found with count = "
		      << vert_el_ptr->vpt->count << "\n";
		  }
		  Sausage newsausage(ij_pair, vert1, vert2, Rprob);
                  sg_plex.push_back( newsausage );
		  vert_el_ptr = vert_el_ptr->next;
		}
	      numb++;
	    }
	  // Otherwise make the donut.
	  else if( ij_pair->flag == Pair::free )
	    {
	      Sausage newsausage(ij_pair, at_i_pt, at_j_pt, Rprob);
	      sg_plex.push_back( newsausage );
	      numb_donut++;
	    }
	}
    }
  sausage_count = sg_plex.size();

  if( sausage_count != 0 )
    {
      sglist = new Sausage [sausage_count];
      std::vector<Sausage>::iterator p;
      for( i = 0, p = sg_plex.begin();
           i < sausage_count;
           i++, p++ )
	sglist[i] = *p;
    }

  int non_bur_sh = sph_count - nburied_shells;
  blab2 << totcones << " in " << non_bur_sh << " shells or ave. "
    << ((float) totcones) / non_bur_sh << " cones/shell " << endl;

  blab2 << "Number of pairs with nonempty linked list is " << numb << endl;
  blab2 << "Number of sausage generated is " << sausage_count << endl;
  blab2 << "Number of donuts generated is " << numb_donut << endl;

  for( i=0; i<sph_count; i++ )
      delete [] pair_array[i];
  delete [] pair_array;
  delete [] pairlist_length_at;

  is_calculated = 1;
#ifdef POPULATION_COUNT
  blab1 << "At end of SolvAccVolRep::anal_calc:" << endl;
  Vertex::population_report();
  VertElem::population_report();
#endif
}

// SAVanal_calc.cc ends here
