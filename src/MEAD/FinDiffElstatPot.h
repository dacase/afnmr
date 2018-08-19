/* Electrostatic potential for a dielectric slab -*- C++ -*-

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

$Id: FinDiffElstatPot.h,v 2.17 2007/02/05 22:17:13 bashford Exp $
*/
#ifndef _FinDiffElstatPot_h
#define _FinDiffElstatPot_h 1

#include <string>
using std::string;

#include "MEAD/ElstatPot.h"
#include "MEAD/FinDiffMethod.h"
#include "MEAD/FDGridLevel.h"

//!wrap!
class FinDiffElstatPot : public ElstatPot_lett {
public:
// FIXME!  How to get in the AnalyticEP of choice?
  //!nowrap!
  FinDiffElstatPot (FinDiffMethod, DielectricEnvironment_lett *,
		    ChargeDist_lett *, ElectrolyteEnvironment_lett*,
		    AnalyticEP*);
  FinDiffElstatPot (FinDiffMethod, DielectricEnvironment_lett *, ChargeDist_lett *, ElectrolyteEnvironment_lett*);
  FinDiffElstatPot (DielectricEnvironment_lett *, ChargeDist_lett *, ElectrolyteEnvironment_lett*);
  virtual ~FinDiffElstatPot ();
  virtual void solve();
  // Use an AVS field for coarse init.
  virtual void solve_using_coarse_init(string fieldname = "");
  void write_coarse_field(const string& fieldname);
  CubeLatSpec coarse_lattice_spec () const;
  //!nowrap!+
  float *coarse_field_into_array (int *size) const;
  FinDiffMethod get_method() {return method;}
  FDGridLevel *get_fine_grid_pot() {return fine_grid_pot;}
  FDGridLevel *get_coarser(FDGridLevel* theLevel) {return theLevel->coarser;}
  void get_val_array(FDGridLevel* theLevel, float* array) {theLevel->get_val_array(array);}
  //!nowrap!-
  virtual float value(Coord) const ;
//  virtual float operator * (const ChargeDist& c) const
//    {return c.finDiffElstatPotMultiply(*this);}
  virtual Coord field(Coord) const ;
  virtual Coord displacement(Coord) const ;
protected:
  FinDiffMethod method;
  FDGridLevel *fine_grid_pot;
  ChargeDist_lett* rho;
  DielectricEnvironment_lett* eps;
  ElectrolyteEnvironment_lett* electrolyte;
  bool this_owns_analytic_approx; // Ugh!
  AnalyticEP* analytic_approx;
  int solved;
};

#if SWIGPP_LITERAL_INCLUDE

%addmethods FinDiffElstatPot {
//
// Multiply with a ChargeDist_lett derived class.
// (The __rmul__ version is defined in the ChargeDist_lett derived classes
//
  float operator* (const FinDiffElstatPot& fde, const ChargeDist_lett& cdl);
//
// Scale operations
//
  ElstatPotCombination operator* (FinDiffElstatPot& fde, float scale);
  ElstatPotCombination operator* (float scale, FinDiffElstatPot& fde);
  ElstatPotCombination operator/ (FinDiffElstatPot& fde, float scale);
//
// Need to override the base class implementation of get_cuberep
// to take advantage of the cubic lattice solutions already in hand.
// Note that the FDGridLevels are accessed from finest -> coarsest
// while the FinDiffMethod levels are accessed from coarsest -> finest

  Numeric_3D_Array_FLOAT * get_cuberep(const CubeLatSpec& cls)
  {
    FDGridLevel* finest = self->get_fine_grid_pot();
    list<FDGridLevel*> levels(0);
    int ngr = cls.get_grid_dim();
    int ncube = ngr * ngr * ngr;
    Numeric_3D_Array_FLOAT * data_array = 0;
    // Look for the cls in the grid levels of this FinDiffElstatPot.
    if (finest != 0) {
      while (finest) {
        levels.push_back(finest);
        finest = self->get_coarser(finest);
      }
      FinDiffMethod theMethod(self->get_method());
      CubeLatSpec *theLattice = const_cast<CubeLatSpec*>(&cls);
      CubeLatSpec *FDlattice = theMethod.get_coarsest();
      int level=0;
      while (FDlattice) {
        if (*theLattice == *FDlattice) {
          list<FDGridLevel*>::reverse_iterator rind=levels.rbegin();
          for (; rind != levels.rend() && level > 0; ++rind, --level);
          if (rind != levels.rend() && level == 0) {
            data_array = new Numeric_3D_Array_FLOAT [ncube];
            FDGridLevel* theLevel = *rind;
            self->get_val_array(theLevel, data_array);
          }
          break;
        }
        FDlattice = theMethod.get_finer();
        ++level;
      }
    }
    // If the cls level wasn't found then use the base class method
    if (data_array == 0) {
      data_array = ElstatPot_lett_get_cuberep(self, cls);
    }
    return data_array;
  }
};

#endif // SWIGPP_LITERAL_INCLUDE


#endif

// FinDiffElstatPot.h ends here
