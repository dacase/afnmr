// This is -*- C++ -*-
#ifndef _Sausage_h
#define _Sausage_h 1

/*
$Id: Sausage.h,v 2.6 2004/12/06 17:58:36 bashford Exp $
*/

/*
    This source code file is part of the MEAD (Macroscopic
    Electrostatics with Sphereic Detail) package of objects and
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

#include "MEAD/Coord.h"
#include "MEAD/AccTag.h"

class CubeLatSpec;
class VertElem;
class Pair;
class Shell;

class Sausage {
public:

  Sausage() {}
//  Sausage(const Coord& c1, const Coord& c2, const Coord& c3,
//	  float ang, float r);
  Sausage (const Pair*, const VertElem *, const VertElem*, float probe_radius);
  Sausage (const Pair*, const Shell*, const Shell*, float probe_radius);
  Sausage(const Sausage& a);
  Sausage& operator= (const Sausage& a);
  void mark_cubelat(const CubeLatSpec&, AccTag []);
  int pt_inside(const Coord&) const ;
  inline Coord get_sausage_head(){return coord_1;}
  inline Coord get_sausage_tail(){return coord_2;}
  inline Coord get_sausage_center(){return coord_3;}
  inline float get_sausage_rad(){return rad;}
  inline float get_sausage_angle(){return angle;}
  ostream& write_top_in_binary(ostream&);
  ostream& write_top_in_ascii(ostream&);
  int read_top_in_binary(istream&);
  int read_top_in_ascii(istream&);
  ~Sausage( ) {}
private:
  Coord   coord_1, coord_2, coord_3;
  float   angle, rad;
  float Rprob;
};

inline
Sausage::Sausage(const Sausage& a)
{
  coord_1 = a.coord_1;
  coord_2 = a.coord_2;
  coord_3 = a.coord_3;
  rad = a.rad;
  angle = a.angle;
  Rprob = a.Rprob;
}

inline Sausage&
Sausage::operator= (const Sausage& a)
{
  coord_1 = a.coord_1;
  coord_2 = a.coord_2;
  coord_3 = a.coord_3;
  rad = a.rad;
  angle = a.angle;
  Rprob = a.Rprob;
  return *this;
}

#endif

// Sausage.h ends here
