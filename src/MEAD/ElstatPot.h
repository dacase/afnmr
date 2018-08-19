/* Base class for electrostatic potential, envelope for subclasses -*- C++ -*-

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

$Id: ElstatPot.h,v 2.11 2004/12/06 17:58:18 bashford Exp $
*/
#ifndef _ElstatPot_h
#define _ElstatPot_h 1

// FIXME!  More use of const-ness is needed to insure that charges
// and dielectric environments do not "change underneath" the ElstatPots
// they are in.

#include "MEAD/Coord.h"
// #include "MEAD/FinDiffMethod.h"
// #include "MEAD/globals.h"

class FDGridLevel;
class CubeLatSpec;
class ChargeDist_lett;
class DielectricEnvironment_lett;
class ElectrolyteEnvironment_lett;
class ChargeDist;
class DielectricEnvironment;
class ElectrolyteEnvironment;
class ElstatPot_lett;
class FinDiffMethod;

//DON'T !wrap!
class ElstatPot {
  //!nowrap!+
friend class AnalyticEP;  // FIXME!  Is this really needed?
  //!nowrap!-
public:
  //!nowrap!
  ElstatPot(DielectricEnvironment, ChargeDist, ElectrolyteEnvironment);
  ElstatPot(FinDiffMethod, DielectricEnvironment, ChargeDist,
	    ElectrolyteEnvironment);
  //!nowrap!+
  explicit ElstatPot (ElstatPot_lett* epp);
  ElstatPot();
  ElstatPot(const ElstatPot&);
  ElstatPot& operator= (const ElstatPot&);
  virtual ~ElstatPot();
  //!nowrap!-
  virtual void solve();
  virtual float value(Coord x) const ;
  virtual Coord field(Coord x) const ;
  virtual Coord displacement(Coord x) const ;
  //!nowrap!+
//  virtual float operator * (const ChargeDist& c) const
//    {return (saferep())->operator* (c);}
  static int instances() {return instanceCount;}
  //!nowrap!-
private:
  ElstatPot_lett *rep;
  static int instanceCount;
};


//!wrap!
class ElstatPot_lett {
  //!nowrap!+
friend class ElstatPot;
  //!nowrap!-
public:
  //!nowrap!+
  ElstatPot_lett ();
  ElstatPot_lett (const ElstatPot_lett&);
  virtual ~ElstatPot_lett();
  ElstatPot_lett& operator= (const ElstatPot_lett&);
  ElstatPot_lett (DielectricEnvironment_lett* e,
		      ChargeDist_lett* r,
		      ElectrolyteEnvironment_lett* ely);
  virtual void solve() = 0;
  virtual float value(Coord c) const = 0;
  virtual Coord field(Coord x) const = 0;
  virtual Coord displacement(Coord x) const = 0;
  static int instances() {return instanceCount;}
  void operator delete (void *p, size_t s);
//  virtual float operator * (const ChargeDist& c) const
//    {return c.MultiplyElstatPot(*this);}
  //!nowrap!-
private:
  ChargeDist* rho_obj;
  DielectricEnvironment* eps_obj;
  ElectrolyteEnvironment* electrolyte_obj;
  static int instanceCount;
  short referenceCount;
};

inline void ElstatPot::solve()
{rep->solve();}

inline float ElstatPot::value(Coord x) const
{return rep->value(x);}

inline Coord ElstatPot::field(Coord x) const
{ return rep->field(x); }

inline Coord ElstatPot::displacement(Coord x) const
{return rep->displacement(x);}

// These work for all envelopes and letters, but more specialized versions
// would be faster.  FIXME.

float operator* (const ChargeDist& c, const ElstatPot& e);
inline float operator* (const ElstatPot& e, const ChargeDist& c)
{ return operator*(c,e); }

float operator* (const ChargeDist& c, const ElstatPot_lett& e);
inline float operator* (const ElstatPot_lett& e, const ChargeDist& c)
{ return operator*(c,e); }

float operator* (const ChargeDist_lett& c, const ElstatPot& e);
inline float operator* (const ElstatPot& e, const ChargeDist_lett& c)
{ return operator*(c,e); }

float operator* (const ChargeDist_lett& c, const ElstatPot_lett& e);
inline float operator* (const ElstatPot_lett& e, const ChargeDist_lett& c)
{ return operator*(c,e); }

//!wrap!
class AnalyticEP : public ElstatPot_lett {
public:
  //!nowrap!+
  AnalyticEP () : ElstatPot_lett () {}
  AnalyticEP (DielectricEnvironment_lett* e, ChargeDist_lett* r,
	      ElectrolyteEnvironment_lett* ely) : ElstatPot_lett (e,r,ely) {}
  //!nowrap!-
};

#if SWIGPP_LITERAL_INCLUDE

%addmethods ElstatPot_lett {
//
// Add ElstatPot's.
// This is used by all the ElstatPot_lett derived classes.
// Note that if any __add__ operations are defined in an ElstatPot_lett
// derived class then this operation will be overridden (masked) and will
// have to be defined explicitly there.
//
  ElstatPotCombination operator+ (ElstatPot_lett& epl1, ElstatPot_lett& epl2);

  // The __radd__ version of the following defined in ElstatPotCombination.
  ElstatPotCombination operator+ (ElstatPot_lett& epl, ElstatPotCombination& epc);
//
// Return a Python NumPy array of potential values given a CubeLatSpec.
// This is used by all the ElstatPot_lett derived classes, and works via
// an implicit downcasted self->value() call.
//
  Numeric_3D_Array_FLOAT * get_cuberep(const CubeLatSpec& cls)
  {
    int grid_dim = cls.get_grid_dim();
    int nsq = grid_dim * grid_dim;
    int ncube = nsq * grid_dim;
    Numeric_3D_Array_FLOAT *array_data = new Numeric_3D_Array_FLOAT [ncube];
    if (array_data != 0) {
      float grlen = (float) (grid_dim - 1);
      float spacing = cls.get_spacing();
      float halfgrlen = grlen/2;
      Coord gridcenter_in_grid(halfgrlen, halfgrlen, halfgrlen);
      Coord gridcenter_in_space = cls.get_center();
      Coord gridpoint_in_space;
      for (int i=0; i<grid_dim; ++i) {
        gridpoint_in_space.x =
          spacing * ((float) i - gridcenter_in_grid.x)
            + gridcenter_in_space.x;
        for (int j=0; j<grid_dim; ++j) {
          gridpoint_in_space.y = spacing * ((float) j - gridcenter_in_grid.y)
            + gridcenter_in_space.y;
          for (int k=0; k<grid_dim; ++k) {
            gridpoint_in_space.z = spacing * ((float) k - gridcenter_in_grid.z)
              + gridcenter_in_space.z;
            int fortind = i + j*grid_dim + k*nsq;
            array_data[fortind] = self->value(gridpoint_in_space);
          }
        }
      }
    }
    return array_data;
  }
};

#endif // SWIGPP_LITERAL_INCLUDE


#endif

// ElstatPot.h ends here
