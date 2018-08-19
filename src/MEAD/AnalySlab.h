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

$Id: AnalySlab.h,v 2.7 2004/12/06 17:57:51 bashford Exp $
*/
#ifndef _AnalySlab_h
#define _AnalySlab_h 1

#include "MEAD/ElstatPot.h"

#define ANALYSLAB_MAXTERM 20

class DielectricSlab;

//!wrap!
class AnalySlab : public AnalyticEP {
public:
  //!nowrap!
  AnalySlab (DielectricSlab*, ChargeDist_lett*, ElectrolyteEnvironment_lett*,
	int maxterm=ANALYSLAB_MAXTERM);  // No. of terms in image expression.
  ~AnalySlab() {}
  virtual void solve();
  virtual float value(Coord) const ;
  virtual Coord field(Coord x) const ;
  virtual Coord displacement(Coord x) const ;
  float get_zlower() const;
  float get_zupper() const;
  int get_maxterm() const;
//  virtual float operator * (const ChargeDist& c) const
//    {return c.analySlabMultiply(*this);}
private: 
  ChargeDist_lett* rho;
  DielectricSlab* eps;
  ElectrolyteEnvironment_lett* electrolyte;
  float epsslab, epsext;
  float zupper, zlower;

  float zmiddle;
  float twob;
  float b;
  float qinfac;
  float qoutfac;
  int maxterm;

  int solved;
};

inline float AnalySlab::get_zlower() const {return zlower;}
inline float AnalySlab::get_zupper() const {return zupper;}
inline int AnalySlab::get_maxterm() const {return maxterm;}

#if SWIGPP_LITERAL_INCLUDE

%addmethods AnalySlab {
//
// swigpp.el strips the default value of the maxterm argument in the
// above constructor. So, wrap this constructor with the default
// argument both present and removed.
  AnalySlab (DielectricSlab*, ChargeDist_lett*, ElectrolyteEnvironment_lett*, int);
  AnalySlab (DielectricSlab*, ChargeDist_lett*, ElectrolyteEnvironment_lett*);
//
// Multiply with a ChargeDist_lett derived class.
// (The __rmul__ version is defined in the ChargeDist_lett derived classes
//
  float operator* (const AnalySlab& asl, const ChargeDist_lett& cdl);
//
// Scale operations
//
  ElstatPotCombination operator* (AnalySlab& asl, float scale);
  ElstatPotCombination operator* (float scale, AnalySlab& asl);
  ElstatPotCombination operator/ (AnalySlab& asl, float scale);
};

%wrapper %{

BEGIN_CPLUSPLUS_SECTION

// Provide the wrapped constructors.
AnalySlab * new_AnalySlab (DielectricSlab* ds, ChargeDist_lett* cd, ElectrolyteEnvironment_lett* ee, int maxterm){
  return new AnalySlab(ds, cd, ee, maxterm);
}
AnalySlab * new_AnalySlab (DielectricSlab* ds, ChargeDist_lett* cd, ElectrolyteEnvironment_lett* ee){
  return new AnalySlab(ds, cd, ee, ANALYSLAB_MAXTERM);
}

END_CPLUSPLUS_SECTION

%}

#endif // SWIGPP_LITERAL_INCLUDE

#endif

// AnalySlab.h ends here
