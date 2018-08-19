/* Class for electrostatic potential combinations -*- C++ -*-

    Copyright (c) 1993--2001 by Donald Bashford

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

$Id: ElstatPotCombination.h,v 2.2 2004/12/06 17:58:17 bashford Exp $
*/
#ifndef _ElstatPotCombination_h
#define _ElstatPotCombination_h 1

#include <list>
using std::list;
using std::pair;

#include "MEAD/ElstatPot.h"
#include "MEAD/Coord.h"
#include "MEAD/globals.h"

// You can't make pair template with a reference to an ElstatPot_lett
// abstract base class.
// So use the ElstatPot envelope class to represent them.

//!wrap!
class ElstatPotCombination : public list<pair<ElstatPot, float> > {
public:
  ElstatPotCombination ();
  ElstatPotCombination (const ElstatPotCombination& epc);
  //!nowrap!+
  ElstatPotCombination (ElstatPot& ep, float scale = 1.0);
  ElstatPotCombination (ElstatPot& ep1, ElstatPot& ep2);
  ElstatPotCombination (ElstatPot& ep, ElstatPot_lett& epl);
  ElstatPotCombination (ElstatPot_lett& epl, ElstatPot& ep);
  ElstatPotCombination (ElstatPot_lett& epl, float scale = 1.0);
  ElstatPotCombination (ElstatPot_lett& epl1, ElstatPot_lett& epl2);
  ElstatPotCombination& operator= (const ElstatPotCombination& epc);
  //!nowrap!-
  void solve();
  float value(Coord c) const;
  Coord field(Coord x) const;
  Coord displacement(Coord x) const;
  ElstatPotCombination& operator+= (const ElstatPotCombination& epc);
  //!nowrap!
  ElstatPotCombination& operator+= (ElstatPot& ep);
  ElstatPotCombination& operator+= (ElstatPot_lett& epl);
  ElstatPotCombination& operator*= (float scale);
  ElstatPotCombination& operator/= (float scale);
  //!nowrap!
  float get_scale() const {return global_scale;}
private:
  float global_scale;
};

// Methods to add and scale ElstatPot's.
// These create a new type of object, an ElstatPotCombination, representing
// a list of ElstatPot's each with an associated scale factor (from
// multiplication or division). Adding ElstatPot's together simply adds them
// to the list. Taking the value of an ElstatPotCombination sums the values
// of the individual ElstatPot's in the list, and applies the scale factors.

ElstatPotCombination operator* (ElstatPot& ep, float scale);
ElstatPotCombination operator* (float scale, ElstatPot& ep);
ElstatPotCombination operator* (ElstatPot_lett& epl, float scale);
ElstatPotCombination operator* (float scale, ElstatPot_lett& epl);
ElstatPotCombination operator/ (ElstatPot& ep, float scale);
ElstatPotCombination operator/ (ElstatPot_lett& epl, float scale);
ElstatPotCombination operator+ (ElstatPot_lett& epl1, ElstatPot_lett& epl2);
ElstatPotCombination operator+ (ElstatPot_lett& epl, ElstatPot& ep);
ElstatPotCombination operator+ (ElstatPot& ep, ElstatPot_lett& epl);
ElstatPotCombination operator+ (ElstatPot& ep1, ElstatPot& ep2);

// You can add more ElstatPot's to an ElstatPotCombination.
// You can globally scale an ElstatPotCombination.

ElstatPotCombination operator+ (const ElstatPotCombination& ep1, const ElstatPotCombination& ep2);
ElstatPotCombination operator+ (const ElstatPotCombination& ep, ElstatPot_lett& epl);
ElstatPotCombination operator+ (ElstatPot_lett& epl, const ElstatPotCombination& ep);
ElstatPotCombination operator+ (const ElstatPotCombination& epc, ElstatPot& ep);
ElstatPotCombination operator+ (ElstatPot& ep, const ElstatPotCombination& epc);
ElstatPotCombination operator* (const ElstatPotCombination& ep, float scale);
ElstatPotCombination operator* (float scale, const ElstatPotCombination& ep);
ElstatPotCombination operator/ (const ElstatPotCombination& ep, float scale);
ElstatPotCombination operator/ (float scale, const ElstatPotCombination& ep);

// You can multiply by a ChargeDist

float operator* (const ChargeDist_lett& cd, const ElstatPotCombination& ep);
float operator* (const ElstatPotCombination& ep, const ChargeDist_lett& cd);

#if SWIGPP_LITERAL_INCLUDE

%addmethods ElstatPotCombination {
//
// Needed to define __del__ to clean up new objects created by python
  ~ElstatPotCombination () {}
//
  ElstatPotCombination operator+ (const ElstatPotCombination& ep1, const ElstatPotCombination& ep2);

  // The __radd__ version is defined in the ElstatPot_lett class.
  ElstatPotCombination operator+ (const ElstatPotCombination& ep, ElstatPot_lett& epl);

  ElstatPotCombination operator* (const ElstatPotCombination& ep, float scale);
  ElstatPotCombination operator* (float scale, const ElstatPotCombination& ep);
  ElstatPotCombination operator/ (const ElstatPotCombination& ep, float scale);

  // The __rmul__ versions are defined in ChargeDist_lett derived classes.
  float operator* (const ElstatPotCombination& epc, const ChargeDist_lett& cd);

/*
 * If you could get at the ElstatPot rep, like
 * ElstatPot_lett* ElstatPot::get_rep() {return rep;}
 * and figure out it's underlying type at runtime,
 * then you could:
 *
 * Implement all the Python list methods
 *
 * Implement the get_cuberep() (like the following untested method)
 *
  Numeric_3D_Array_FLOAT * get_cuberep(const CubeLatSpec& cls)
  {
    int ngr = cls.get_grid_dim();
    int ncube = ngr * ngr * ngr;
    Numeric_3D_Array_FLOAT * data_array = new Numeric_3D_Array_FLOAT [ncube];
    Numeric_3D_Array_FLOAT * tmp_array = 0;
    bzero((void *) data_array, ncube * sizeof(Numeric_3D_Array_FLOAT));
    for (ElstatPotCombination::const_iterator ind = self->begin(); ind != self->end(); ++ind) {
      pair<ElstatPot, float> p = *ind;
      tmp_array = ElstatPot_lett_get_cuberep(p.first.get_rep(), cls);
      if (tmp_array != 0) {
        for (int i=0; i<ncube; ++i) {
          array_data[i] += p.second * tmp_array[i];
        }
        delete [] tmp_array;
        tmp_array = 0;
      }
    }
    float global_scale = self->get_scale();
    for (int i=0; i<ncube; ++i) {
      array_data[i] *= global_scale;
    }
    return array_data;
  }
 *
 */
};

#endif // SWIGPP_LITERAL_INCLUDE


#endif

// ElstatPotCombination.h ends here
