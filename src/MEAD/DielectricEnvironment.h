/* Base class for Diel. Evir. and envelope for subclasses -*- C++ -*-

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

$Id: DielectricEnvironment.h,v 2.10 2004/12/06 17:58:13 bashford Exp $
*/
#ifndef _DielectricEnvironment_h
#define _DielectricEnvironment_h 1


#include "MEAD/CubeLatSpec.h"
#include "MEAD/DielCubeRep.h"
#include "MEAD/Coord.h"

class DielectricEnvironment_lett;
class ElstatPot;

//DON'T !wrap!
class DielectricEnvironment {
// LANGUAGE BUG!!  The following friend declaration is OK with gcc
// but not cfront.  According to the ARM Sec 11.4 p. 251, cfront is right!
//    friend ElstatPot::ElstatPot (DielectricEnvironment, ChargeDist,
//                                 ElectrolyteEnvironment);
// So must get "overly friendly"...
  //!nowrap!+
  friend class ElstatPot;
  //!nowrap!-
public:
  // Canonical stuff
  //!nowrap!+
  DielectricEnvironment();
  DielectricEnvironment(const DielectricEnvironment& d);
  //!nowrap!-
  DielectricEnvironment(DielectricEnvironment_lett*);
  //!nowrap!+
  DielectricEnvironment& operator= (const DielectricEnvironment&);
  virtual ~DielectricEnvironment();

  virtual DielCubeRep get_cuberep(const CubeLatSpec& cls);
  virtual float epsext_value() const ;
  static int instances() {return instanceCount;}
  //!nowrap!-
private:
  DielectricEnvironment_lett *rep;
  static int instanceCount;
};

//!wrap!
class DielectricEnvironment_lett {
  //!nowrap!+
friend class DielectricEnvironment;
  //!nowrap!-
public:
  //!nowrap!+
  DielectricEnvironment_lett();
  DielectricEnvironment_lett(const DielectricEnvironment_lett&);
  DielectricEnvironment_lett& operator= (const DielectricEnvironment_lett&);
  virtual ~DielectricEnvironment_lett();
  virtual DielCubeRep get_cuberep(const CubeLatSpec& clsp) = 0;
  virtual float epsext_value() const = 0;
  static int instances() {return instanceCount;}
  void operator delete(void *p, size_t s);
  //!nowrap!-
private:
  static int instanceCount;
  short referenceCount;
};


inline DielCubeRep
DielectricEnvironment::get_cuberep(const CubeLatSpec& cls)
{return rep->get_cuberep(cls);}

inline float
DielectricEnvironment::epsext_value() const
{return rep->epsext_value();}

#if SWIGPP_LITERAL_INCLUDE

%addmethods DielectricEnvironment_lett {
//
// Return a Numeric_3D_Array_FLOAT of eps values
// This is used by all the DielectricEnvironment_lett derived classes,
// and works via an implicit downcasted self->get_cuberep() call.
//

  Numeric_3D_Array_FLOAT * get_cuberep(const CubeLatSpec& cls)
  {
    int ngr = cls.get_grid_dim();
    int nsq = ngr * ngr;
    int ncube = nsq * ngr;
    Numeric_3D_Array_FLOAT *array_data = new Numeric_3D_Array_FLOAT [ncube];
    if (array_data != 0) {
      const DielCubeRep dcr = self->get_cuberep(cls);
      int h = 0;
      for (int i = 0; i < ngr; ++i) {
        for (int j = 0; j < ngr; ++j) {
          for (int k = 0; k < ngr; ++k) {
            int fortind = i + j*ngr + k*nsq;
            array_data[fortind] = dcr[h];
            ++h;
          }
        }
      }
    }
    return array_data;
  }
};

#endif // SWIGPP_LITERAL_INCLUDE


#endif

// DielectricEnvironment.h ends here
