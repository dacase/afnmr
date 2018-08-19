// This is -*- C++ -*-
/* Specification of the sequence of grids to use for finite diff. proc.

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

$Id: FinDiffMethod.h,v 2.8 2004/12/06 17:58:26 bashford Exp $
*/
#ifndef _FinDiffMethod_h
#define _FinDiffMethod_h 1

#include "MEAD/Coord.h"
#include "MEAD/CenteringStyle.h"
#include <string>
using std::string;

class CubeLatSpec;
struct FDMLink;

class FinDiffMethodRep {
friend class FinDiffMethod;
public:
  FinDiffMethodRep();
  void read (const string& filename);
  void add_level(int ngrd, float spc, Coord cntr);
  //!nowrap!
  void add_level(int ngrd, float spc, CenteringStyle);
  void resolve (Coord geom_cent, Coord center_of_intr);
  CubeLatSpec * get_coarsest();
  CubeLatSpec * get_finer();
  ostream& print(ostream&) const;
private:
  FDMLink *coarsest, *finest, *ptr;
  short resolved;
  short referenceCount;
};

//!wrap!
class FinDiffMethod {
public:
  FinDiffMethod();
  FinDiffMethod(const FinDiffMethod& f);
  //!nowrap!
  FinDiffMethod & operator= (const FinDiffMethod& f);
  ~FinDiffMethod();
  void read (const string& filename);
  void add_level (int ngrd, float spc, Coord cntr);
  //!nowrap!
  void add_level(int ngrd, float spc, CenteringStyle);
  void resolve (Coord geom_cent, Coord center_of_intr);
  //!nowrap!+
  CubeLatSpec * get_coarsest();
  CubeLatSpec * get_finer();
  ostream& print(ostream&) const;
  //!nowrap!-
  bool is_resolved() const { return rep->resolved != 0; }
private:
  FinDiffMethodRep * rep;
};

ostream & operator<< (ostream& ost, const FinDiffMethod& fdm);

#if SWIGPP_LITERAL_INCLUDE

%addmethods FinDiffMethod {
//
// Rewrap this (overloaded) method with CenteringStyle as an int
  void add_level(int ngrd, float spc, int sty);
//
// Python doesn't want you to redefine the print method, so call it write
  void write() { self->print(cout); }
//
// Return a list of CubeLatSpec levels in fdm
  list_CubeLatSpec * levels()
  {
    list_CubeLatSpec *lcls = new list_CubeLatSpec(0);
    CubeLatSpec *cls = 0;
    if (cls = self->get_coarsest()) {
      lcls->push_back(*cls);
      while (cls = self->get_finer()) {
        lcls->push_back(*cls);
      }
    }
    return lcls;
  }

};

%wrapper %{
//
// The wrapped (overloaded) add_level method for CenteringStyle as int
// Note: There is no range checking for the enum assignment, so the
//       result could be undefined if an out-of-range int is used.
//
void FinDiffMethod_add_level(FinDiffMethod* self, int ngrd, float spc, int sty)
{
  CenteringStyle style = CenteringStyle(sty);
  self->add_level(ngrd, spc, style);
}

%}

#endif // SWIGPP_LITERAL_INCLUDE

#endif

// FinDiffMethod.h ends here
