// This is -*- C++ -*-
#ifndef _TPShell_h
#define _TPShell_h 1

/*
$Id: TPShell.h,v 2.6 2007/05/28 01:26:43 bashford Exp $
*/

/*
    This source code file is part of the MEAD (Macroscopic
    Electrostatics with Atomic Detail) package of objects and
    programs, copyright (C) 1990 by Donald Bashford of the Department
    Molecular Biology of The Scripps Research Institute.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 1, or (at your option)
    any later version.

    The GNU General Public License should be in a file called COPYING.GNU.
    Some comments about copying and distribution by D. Bashford should
    be in a file called COPYING.DB

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

    Donald Bashford can be contacted by electronic mail by the address

       bashford@scripps.edu

    or by paper mail at

       Department of Molecular Biology
       The Scripps Research Institute
       10666 North Torrey Pines Road
       La Jolla, California  92037
*/

class Shell;
class CubeLatSpec;
struct Cone;
#include "MEAD/Coord.h"
#include "MEAD/AccTag.h"

class TPShell {
public:
  TPShell();
  TPShell(const Shell& sph, int npts, const Coord* coord);
  TPShell(const Shell& sph, int nidx, const int* idx);
  TPShell(const TPShell&);
  TPShell& operator = (const TPShell&);
  inline int has_pts() const  {return num_inbox>0;}
  void mark_by_radii(const Coord *pt, AccTag *acc_array) const ;
  void mark_by_cones(const Coord *pt, AccTag *acc_array) const ;
  ~TPShell();

private:
  // These same as Shell
  Coord coord;
  float inner_rad, outer_rad, outer_rad_sq;
  enum Flag {free=0, partially_buried=1, buried=2 };
  Flag flag;
  int numcones;

  // These new.
  int num_inbox;
  int *inbox_pt_index;
  Cone *cone;
};


inline TPShell::TPShell()
: inner_rad(0), outer_rad(0), outer_rad_sq(0),
  numcones(0), flag(free), num_inbox(0), cone(0), inbox_pt_index(0)
{
}

#endif

// TPShell.h ends here
