/* -*- C++ -*-
Dielectric environment due to a membrane plus some Atoms,
eg. a membrane protein.

    Copyright (c) 1995 by Donald Bashford

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

$Id: DielMembAtoms.h,v 2.7 2001/06/26 18:32:21 bergsma Exp $
*/

#ifndef DielMembAtoms_h
#define DielMembAtoms_h

#include "MEAD/DielectricSlab.h"
#include "MEAD/SolvAccVol.h"
#include "MEAD/AtomSet.h"

//!wrap!
class TwoValueDielMembAtoms : public DielectricEnvironment_lett {
public:
  TwoValueDielMembAtoms(const AtomSet& ats, float epsin, float zl, float zu, Coord c, float holrad);
  TwoValueDielMembAtoms(const SolvAccVol& sav, float epsin, float zl, float zu, Coord c, float holrad);
  ~TwoValueDielMembAtoms() {}
  //!nowrap!+
  DielCubeRep get_cuberep(const CubeLatSpec& cls);
  DielectricSlab* get_slab() {return membranept;}  // FIXME?  Dangerous?
  //!nowrap!-
  float epsext_value() const {return epsext;}
private:
  float epsin, epsext;
  float epsmemb;
  DielectricEnvironment membrane;
  DielectricSlab* membranept;
  Coord hole_center;
  float hole_radius;
// FIXME!  If this were a ref or ptr, not need include SolvAccVol.h
  SolvAccVol solv_acc;
};

// Given membrane boundaries and two AtomSets, region within atset1
// has ein1, region within atset2 (but outside of atset1) has ein2
// region outside of both atsets, but inside the membrane boundaries
// also has ein2, region outside of both atsets and outside membran
// boundaries gets epsext.
//!wrap!
class ThreeValueDielMembAtomsAtoms : public DielectricEnvironment_lett {
public:
  ThreeValueDielMembAtomsAtoms(const AtomSet& a1, float solrad1, const AtomSet& a2, float solrad2, float zl, float zu, Coord c, float holrad, float ein1, float ein2, float epsext);
  ~ThreeValueDielMembAtomsAtoms() {}
  //!nowrap!+
  DielCubeRep get_cuberep(const CubeLatSpec& cls);
  DielectricSlab* get_slab() {return membranept;}  // FIXME?  Dangerous?
  //!nowrap!-
  float epsext_value() const {return epsext;}
private:
  DielectricEnvironment membrane;
  DielectricSlab* membranept;
  Coord hole_center;
  float hole_radius;
  float epsin1, epsin2, epsext;
  SolvAccVol solv_acc1, solv_acc2;
};
 
#endif

// DielMembAtoms.h ends here
