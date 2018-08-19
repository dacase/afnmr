/* Base class for Elecly. Env. and envelope for  subclasses -*- C++ -*-
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

$Id: ElectrolyteEnvironment.h,v 2.12 2004/12/06 17:58:15 bashford Exp $
*/
#ifndef _ElectrolyteEnvironment_h
#define _ElectrolyteEnvironment_h 1

#include <iostream>
#include "MEAD/Coord.h"

class CubeLatSpec;

class ElyCubeRep {
  // Why make the ctor private and require all these frieds?  FIXME?
friend class ElectrolyteEnvironment;
friend class UniformElectrolyte;
friend class ElectrolyteByAtoms;
friend class ElySphere;
public:
  inline int is_guoy(int h) {return isgarr[h];}
  ~ElyCubeRep();
private:
  ElyCubeRep(const CubeLatSpec& s);
  int *isgarr;
};

class ElstatPot;

class ElectrolyteEnvironment_lett;

//DON'T !wrap!
class ElectrolyteEnvironment {
  //!nowrap!+
friend class ElstatPot;
  //!nowrap!-
public:
  //!nowrap!+
  ElectrolyteEnvironment ();
  ElectrolyteEnvironment (const ElectrolyteEnvironment&);
  ElectrolyteEnvironment& operator= (const ElectrolyteEnvironment&);
  //!nowrap!-
  ElectrolyteEnvironment (ElectrolyteEnvironment_lett*);
  //!nowrap!+
  virtual float ionic_strength() const ;
  virtual ElyCubeRep* get_cuberep(const CubeLatSpec& cls) const ;
  virtual ~ElectrolyteEnvironment();
  static int instances() {return instanceCount;}

  //!nowrap!-

private:
  ElectrolyteEnvironment_lett *rep;
  short referenceCount;
  static int instanceCount;
};

//!wrap!
class ElectrolyteEnvironment_lett {
  //!nowrap!+
friend class ElectrolyteEnvironment;
  //!nowrap!-
public:
  //!nowrap!+
  ElectrolyteEnvironment_lett();
  ElectrolyteEnvironment_lett(const ElectrolyteEnvironment_lett&);
  ElectrolyteEnvironment_lett& operator= (const ElectrolyteEnvironment_lett&);
  virtual ~ElectrolyteEnvironment_lett();
  virtual float ionic_strength() const = 0;
  virtual ElyCubeRep* get_cuberep(const CubeLatSpec& cls) const = 0;
  static int instances() {return instanceCount;}
  void operator delete(void *p, size_t s);
  //!nowrap!-
private:
  static int instanceCount;
  short referenceCount;
};

inline float
ElectrolyteEnvironment::ionic_strength() const
{return rep->ionic_strength();}

inline ElyCubeRep*
ElectrolyteEnvironment::get_cuberep(const CubeLatSpec& cls) const
{return rep->get_cuberep(cls);}


#if SWIGPP_LITERAL_INCLUDE

%addmethods ElectrolyteEnvironment_lett {
//
// Return a Numeric_3D_Array_INT of AccTag's (interior or exterior regions)
// This is used by all the ElectrolyteEnvironment_lett derived classes,
// and works via an implicit downcasted self->get_cuberep() call.

  Numeric_3D_Array_INT * get_cuberep(const CubeLatSpec& cls)
  {
    int ngr = cls.get_grid_dim();
    int nsq = ngr * ngr;
    int ncube = nsq * ngr;
    Numeric_3D_Array_INT *array_data = new Numeric_3D_Array_INT [ncube];
    if (array_data != 0) {
      ElyCubeRep *ecr = self->get_cuberep(cls);
      int h = 0;
      for (int i = 0; i < ngr; ++i) {
        for (int j = 0; j < ngr; ++j) {
          for (int k = 0; k < ngr; ++k) {
            int fortind = i + j*ngr + k*nsq;
            array_data[fortind] = ecr->is_guoy(h);
            ++h;
          }
        }
      }
      delete ecr;
    }
    return array_data;
  }
};

#endif // SWIGPP_LITERAL_INCLUDE


#endif

// ElectrolyteEnvironment.h ends here
