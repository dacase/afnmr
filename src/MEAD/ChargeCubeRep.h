/* Charge on a cubic lattice suitable for Finite Diff calcs. -*- C++ -*-

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

$Id: ChargeCubeRep.h,v 2.5 2007/05/28 01:26:42 bashford Exp $
*/
#ifndef _ChargeCubeRep_h
#define _ChargeCubeRep_h 1

#include <iostream>
#include "MEAD/CubeLatSpec.h"
#include "MEAD/globals.h"

class ChargeCubeRep {
public:
  virtual ~ChargeCubeRep () {};
  virtual void add(int h, float ch) = 0;
  virtual int firstindex() = 0;
  virtual int nextindex() = 0;
  virtual int valid(int) = 0;
  virtual float get(int) = 0;
  virtual int is_sparse() {return is_sprs;}
  virtual int num_charged_points() = 0;
  virtual CubeLatSpec get_cubelatspec() {return cls;}
protected:
  int is_sprs;
  CubeLatSpec cls;
};

class SparseChargeCubeRep : public ChargeCubeRep {
public:
  virtual ~SparseChargeCubeRep ();
  SparseChargeCubeRep(CubeLatSpec, size_t num_charged_points);
  virtual void add(int h, float ch);
  virtual int firstindex();
  virtual int nextindex();
  virtual int valid(int);
  virtual float get(int);
  virtual int num_charged_points();
private:
  struct ChargedPoint {
    int h;
    float ch;
  };
  size_t lat_cube;
  ChargedPoint *chrp;
  ChargedPoint *last_added;
  ChargedPoint *last_indexed;
  ChargedPoint *highest_allowed;
  int iscomplete;
};


inline void SparseChargeCubeRep::add(int h, float ch)
{
  if (last_added >= highest_allowed)
    ::error("ERROR: Higher than highest!");
  else {
    ++last_added;
    last_added->h = h;
    last_added->ch = ch;
  }
}

inline int SparseChargeCubeRep::firstindex()
{
  if (last_added < chrp) // case of an empty one
    return lat_cube;
  else {
    last_indexed = chrp;
    return last_indexed->h;
  }
}

inline int SparseChargeCubeRep::nextindex()
{
  ++last_indexed;
  if (last_indexed > last_added)
    return lat_cube;
  else
    return last_indexed->h;
}

inline int SparseChargeCubeRep::valid(int h)
{
  if (last_added >= chrp && last_indexed <= last_added && h==last_indexed->h)
    return 1;
  else
    return 0;
}

inline float SparseChargeCubeRep::get(int h)
{
  if (valid(h))
    return last_indexed->ch;
  else {
    ::error("ERROR: attempted to get illegally!");
    return  0.0;
  }
}

inline int SparseChargeCubeRep::num_charged_points()
{
  return (last_added - chrp) + 1;
}


#endif

// ChargeCubeRep.h ends here
