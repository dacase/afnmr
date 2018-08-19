// This is -*- C++ -*-
/* class for calculations relating to solvent-accessible volumes of molecules
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

$Id: SolvAccVol.h,v 2.9 2004/12/06 17:58:40 bashford Exp $
*/
#ifndef _SolvAccVol_h
#define _SolvAccVol_h 1

#include <string>
using std::string;
#include <iostream>
using std::ostream;
using std::istream;

// Classes used for implementation and communication..

#include "MEAD/AccTag.h"
class CubeLatSpec;
class Sausage;
class Shell;
class AtomSet;
class Coord;

class SolvAccVolRep {
friend class SolvAccVol;
public:
  static int instances();
private:
  SolvAccVolRep(const AtomSet&, float probe_radius);
  ~SolvAccVolRep();
  void anal_calc();
  void calc_cuberep(const CubeLatSpec&, AccTag* );
  void tag_points(int npts, const Coord*, AccTag*);
  int accessible(const Coord&);
  int check_is_calculated() const ;

  // Filename derived from name data member.
  void  write_top_in_binary() const ;
  void  write_top_in_ascii() const ;
  int  read_top_in_binary();
  int  read_top_in_ascii();

  void  write_top_in_binary(ostream&) const ;
  void  write_top_in_ascii(ostream&) const ;
  int  read_top_in_binary(istream&);
  int  read_top_in_ascii(istream&);

  void  write_top_in_binary(const string& filename) const ;
  void  write_top_in_ascii(const string& filename) const ;
  int  read_top_in_binary(const string& filename);
  int  read_top_in_ascii(const string& filename);

  static int instanceCount;
  int referenceCount;
  const float Rprob;
  int  sph_count;
  int  sausage_count;
  Shell  *atsph;
  Sausage *sglist;
  string  name;
  int  is_calculated;
};

//!wrap!
class SolvAccVol {
public:
  SolvAccVol(const AtomSet&);
  SolvAccVol(const AtomSet&, float probe_radius);
  SolvAccVol(const SolvAccVol&);
  //!nowrap!
  SolvAccVol& operator= (const SolvAccVol&);
  ~SolvAccVol();
  void anal_calc();
  //!nowrap!+
  void calc_cuberep(const CubeLatSpec&, AccTag* ) ;
  void tag_points(int npts, const Coord*, AccTag*) ;
  //!nowrap!-
  int accessible(const Coord&) ;
  int check_is_calculated() const ;

  // Filename derived from name data member.
  void  write_top_in_binary() const ;
  void  write_top_in_ascii() const ;
  int  read_top_in_binary();
  int  read_top_in_ascii();

  //!nowrap!+
  void  write_top_in_binary(ostream&) const ;
  void  write_top_in_ascii(ostream&) const ;
  int  read_top_in_binary(istream&);
  int  read_top_in_ascii(istream&);
  //!nowrap!-

  void  write_top_in_binary(const string& filename) const ;
  void  write_top_in_ascii(const string& filename) const ;
  int  read_top_in_binary(const string& filename);
  int  read_top_in_ascii(const string& filename);
private:
  SolvAccVol() {}
  SolvAccVolRep *rep;
};

#if SWIGPP_LITERAL_INCLUDE

%addmethods SolvAccVol {
//
// Return a Numeric_3D_Array_INT of AccTag's

  Numeric_3D_Array_INT * get_cuberep(const CubeLatSpec& cls)
  {
    int ngr = cls.get_grid_dim();
    int nsq = ngr * ngr;
    int ncube = nsq * ngr;
    Numeric_3D_Array_INT *array_data = new Numeric_3D_Array_INT [ncube];
    if (array_data != 0) {
      self->calc_cuberep(cls, (AccTag *) array_data);
      // Swap from C array to Fortran array representation
      for (int h = 0, i = 0; i < ngr; i++) {
        for (int j = 0; j < ngr; j++) {
          for (int k = 0; k < ngr; j++) {
            int fortind = i + j*ngr + k*nsq;
            if (fortind != h) {
              int tmp = array_data[fortind];
              array_data[fortind] = array_data[h];
              array_data[h] = tmp;
            }
            ++h;
          }
        }
      }
    }
    return array_data;
  }

//
// Mark an arbitrary set of points as to solvent accessibility,
// returning a 1D Array

  Numeric_1D_Array_INT * tag_points(int npts, const Coord * pt)
  {
    Numeric_1D_Array_INT *array_data = new Numeric_1D_Array_INT [npts];
    if (array_data != 0) {
      self->tag_points(npts, pt, (AccTag *) array_data);
    }
    return array_data;
  }
};

// This wraps up the AccTag enum, although it creates an entirely
// new enum here. Take care to update this AccTag_enum class as well as
// the python wrapper code (AccTag.py) if the C++ AccTag enum is changed.

%inline %{

class AccTag_enum {
public:
  enum { interior , exterior , undecided , in_tube };
};

%}

// Python code to create a python AccTag instance
%pragma(python) include="AccTag.py"

#endif // SWIGPP_LITERAL_INCLUDE

#endif

// SolvAccVol.h ends here
