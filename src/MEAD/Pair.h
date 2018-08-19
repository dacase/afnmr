/* A pair of shells whose outer radii overlap.  -*- C++ -*-

    Copyright (c) 1993--1995 by Donald Bashford

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

   $Id: Pair.h,v 2.8 2007/05/28 01:26:42 bashford Exp $
*/
#ifndef _Pair_h
#define _Pair_h 1

#include "MEAD/VertElem.h"
#include "MEAD/Shell.h"
// On msvc the define below pulls in macro defs line M_PI
#define _USE_MATH_DEFINES 1
#include <math.h>

class Pair {
  public :
  //   constructors and assignment
  Pair () : vert_head( 0 ) { }
  Pair (int i_index, const Shell& shell_i,
	int j_index, const Shell& shell_j);
  Pair ( const Pair& a );
  Pair& operator = ( const Pair& a );

  //   other member functions
  inline int i_index() const {return i_ndx;}
  inline int j_index() const {return j_ndx;}
  void Check_vertex_pair(const VertElem& v1, const VertElem& v2);

  inline VertElem* get_vert_head( void ) { return vert_head; }
  inline float get_overlap_rad_sq() const {return overlap_rad_sq;}
  inline int quickcheck_burial_by_shell(const Shell&);
  inline Coord get_overlap_center( void ) const { return overlap_center; }
  inline int is_fully_buried() const
    {return flag == buried || flag == geometrically_buried;}
  ~Pair( );

// A publicly accessible flag (hmmm?)
  enum Flag {free=0, partially_buried=1, buried=2, geometrically_buried=3 };
  Flag flag;

  private :
  inline void delete_vertex_ring();
  inline float vertex_angle(const VertElem& v) const;
  int i_ndx, j_ndx;
  float overlap_rad_sq;
  Coord ref_vector;
  Coord overlap_center;
  Coord Rij;
  VertElem *vert_head;

};

inline
Pair::Pair (int i_index, const Shell& shell_i,
	    int j_index, const Shell& shell_j)
{
  flag = free;
  i_ndx = i_index;
  j_ndx = j_index;

  Rij = shell_j.get_coord() - shell_i.get_coord();
  float sij = Rij*Rij;
  float temp = 0.25F/sij;
  float si = sij + shell_i.get_outer_rad_sq() - shell_j.get_outer_rad_sq();
  si = si*si*temp;
  float sj = sij + shell_j.get_outer_rad_sq() - shell_i.get_outer_rad_sq();
  sj = sj*sj*temp;
  float temp2;
  if( sij<sj && shell_i.get_outer_rad_sq() < shell_j.get_outer_rad_sq() )
    temp2 = -sqrt(si/sij);
  else
    temp2 = sqrt(si/sij);
  overlap_center = shell_i.get_coord() + temp2*Rij;
  overlap_rad_sq = shell_i.get_outer_rad_sq() - si;

  vert_head = 0;
}

inline
Pair :: Pair( const Pair& a )
{
  if (a.vert_head)
    ::error("error: Pair copy constructor may not be used\n",
	    "       to copy a Pair that has a non-empty vertex ring.");
  flag = a.flag;
  i_ndx = a.i_ndx;
  j_ndx = a.j_ndx;
  overlap_center = a.overlap_center;
  overlap_rad_sq = a.overlap_rad_sq;
  ref_vector = a.ref_vector;
  Rij = a.Rij;
  vert_head = 0;
}

inline
Pair& Pair :: operator = (const Pair& a)
{
  if (a.vert_head || vert_head)
    ::error("error: Pair assignment operator may not be used to assign\n",
	    "       to or from a Pair that has a non-empty vertex ring.");
  flag = a.flag;
  i_ndx = a.i_ndx;
  j_ndx = a.j_ndx;
  overlap_center = a.overlap_center;
  overlap_rad_sq = a.overlap_rad_sq;
  ref_vector = a.ref_vector;
  Rij = a.Rij;
  return *this;
}

inline
Pair :: ~Pair( )
{
  delete_vertex_ring();
}


inline void
Pair::delete_vertex_ring()
{
  for (VertElem *current = vert_head; current; ) {
    VertElem *temp = current->next;
    delete current;
    current = temp;
  }
  vert_head = 0;
}

// This function will give the angle of a vertex on the pairs ring
// provided it is really on the ring and ref_vector is defined (no checking).

inline float
Pair::vertex_angle(const VertElem& v) const
{
  float ang;
  Coord r = v.vpt->coord - overlap_center;
  float cos_ang = r*ref_vector/overlap_rad_sq;
  if( fabs(cos_ang)-1 > 0 ) {
    if( cos_ang > 0 )
      ang = 0;
    else
      ang = 3.1415926F;
  }
  else {
    float  dot_and_cross_product = Rij*cross( r,ref_vector );
    if( dot_and_cross_product > 0 )
      ang = acos( cos_ang );
    else
      ang = 6.2831852F - acos( cos_ang );
  }
  return ang;
}

/* WARNING!  If this Pair is already buried this function doesn't bother
   to check whether the pair is also buried by other_shell, it just returns
   zero.  So it cannot be used for general purpose checks of geometric
   burial.  Also it depends for its validity on the spcor3 problem having
   no solutions.
*/

inline int
Pair::quickcheck_burial_by_shell(const Shell& other_shell)
{
  int retval = 0;
  if (flag!=geometrically_buried) {
    Coord cent_to_shell = overlap_center - other_shell.get_coord();
    float a = cent_to_shell*cent_to_shell;
    float b = other_shell.get_outer_rad_sq() - overlap_rad_sq;
    if (a<b) {
      flag = geometrically_buried;
      delete_vertex_ring();
      retval = 1;
    }
  }
  return retval;
}
#endif

// Pair.h ends here
