/* Electrostatic potential for a dielectric sphere -*- C++ -*-

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

$Id: AnalySphere.h,v 2.8 2007/05/28 01:26:41 bashford Exp $
*/
#ifndef AnalySphere_h
#define AnalySphere_h

#include "MEAD/ElstatPot.h"
#include "MEAD/ElySphere.h"
#include "MEAD/Polynomial.h"
#include <vector>

#define ANALYSPHERE_MAXTERM 50

class DielectricSphere;

//!wrap!
class AnalySphere : public AnalyticEP {
public:
  //!nowrap!
  AnalySphere(DielectricSphere*, ChargeDist_lett*, ElySphere*,
	int maxterm=ANALYSPHERE_MAXTERM);  // No. of polynomial terms
  ~AnalySphere() {}
  virtual void solve();
  virtual float value(Coord c) const ;
  virtual Coord field(Coord x) const ;
  virtual Coord displacement(Coord x) const ;
  float get_rad_diel() const;
  float get_rad_ely() const;
  int get_maxterm() const;
private:
  DielectricSphere* eps;
  ChargeDist_lett* rho;
  ElectrolyteEnvironment_lett* electrolyte;
  float epsin, epsext;
  float rad_diel;
  float rad_ely;
  Coord center;
  float ionic_strength;
  float kappa;
  bool solved;
  bool central_charge_only_nosalt;
  float cencharge;
  float monopole_reacfield;
  vector<Polynomial> legendre;
  //vector<Polynomial> kirkwood;
  Polynomial* kirkwood;
  vector<double> Ka;
  vector<double> Pa;
  //vector<double> D;
  static const double epsilon;
  int l;
  int size;
  double* A;
  double* B;
  double* C;
  double* D;
};

inline float AnalySphere::get_rad_diel() const {return rad_diel;}
inline float AnalySphere::get_rad_ely() const  {return rad_ely;}
inline int AnalySphere::get_maxterm() const {return l;}

#if SWIGPP_LITERAL_INCLUDE

%addmethods AnalySphere {
//
// swigpp.el strips the default value of the maxterm argument in the
// above constructor. So, wrap this constructor with the default
// argument both present and removed.
  AnalySphere(DielectricSphere*, ChargeDist_lett*, ElySphere*, int);
  AnalySphere(DielectricSphere*, ChargeDist_lett*, ElySphere*);
//
// Multiply with a ChargeDist_lett derived class.
// (The __rmul__ version is defined in the ChargeDist_lett derived classes
//
  float operator* (const AnalySphere& asp, const ChargeDist_lett& cdl);
//
// Scale operations
//
  ElstatPotCombination operator* (AnalySphere& asp, float scale);
  ElstatPotCombination operator* (float scale, AnalySphere& asp);
  ElstatPotCombination operator/ (AnalySphere& asp, float scale);
};

%wrapper %{

BEGIN_CPLUSPLUS_SECTION

// Provide the wrapped constructors.
AnalySphere * new_AnalySphere(DielectricSphere* ds, ChargeDist_lett* cd, ElySphere* esp, int l){
  return new AnalySphere(ds, cd, esp, l);
}
AnalySphere * new_AnalySphere(DielectricSphere* ds, ChargeDist_lett* cd, ElySphere* esp){
  return new AnalySphere(ds, cd, esp, ANALYSPHERE_MAXTERM);
}

END_CPLUSPLUS_SECTION

%}

#endif // SWIGPP_LITERAL_INCLUDE

#endif

// AnalySphere.h ends here
