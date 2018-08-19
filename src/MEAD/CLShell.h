// This is -*- C++ -*-
#ifndef _CLShell_h
#define _CLShell_h 1

/*
$Id: CLShell.h,v 2.6 2007/05/28 01:26:42 bashford Exp $
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

class CLShell {
public:
  CLShell();
  CLShell(const Shell& sph, const CubeLatSpec& cls);
  CLShell(const CLShell&);
  CLShell& operator = (const CLShell&);
  inline int in_lattice() {return is_in;}
  void mark_by_radii(AccTag *acc_array);
  void mark_by_cones(AccTag *acc_array);
  ~CLShell();

private:
  // These same as Shell
  Coord coord;
  float inner_rad, outer_rad, outer_rad_sq;
  enum Flag {free=0, partially_buried=1, buried=2 };
  Flag flag;
  int numcones;

  // These new.
  int grid_dim;
  int i1, i2, j1, j2, k1, k2;  // grid box around center
  int is_in;
  Cone *cone;
};


inline CLShell::CLShell() 
{
  inner_rad = outer_rad = outer_rad_sq = 0;
  numcones = 0;
  flag = free;
  is_in = 0;
  i1 = i2 = j1 = j2 = k1 = k2 = 0;
  grid_dim = 0;
  cone = 0;
}

#endif

// CLShell.h ends here
