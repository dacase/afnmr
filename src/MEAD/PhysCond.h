// This is -*- C++ -*-
#ifndef _PhysCond_h
#define _PhysCond_h 1

/*
$Id: PhysCond.h,v 1.11 2007/05/28 01:26:42 bashford Exp $
*/

/*
    This source code file is part of the MEAD (Macroscopic
    Electrostatics with Atomic Detail) package of objects and
    programs, copyright (C) 1990 by Donald Bashford of the Department
    Molecular Biology of The Scripps Research Institute.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 1, or (at your option)
    any later version.

    The GNU General Public License should be in a file called COPYING.GNU.
    Some comments about copying and distribution by D. Bashford should
    be in a file called COPYING.DB

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

    Donald Bashford can be contacted by electronic mail by the address

       bashford@scripps.edu

    or by paper mail at

       Department of Molecular Biology
       The Scripps Research Institute
       10666 North Torrey Pines Road
       La Jolla, California  92037
*/

/*
PhysCond is a purely static class containing parameters describing
the physical conditions (dielectric constants, ionic strength,
solvent probe radius, etc.), some units conversion factors, and
the necessary access functions to set and get these values.
The main program should have a section that reads command line flags
to set these values.
*/

#include "MEAD/globals.h"

//!wrap!
class PhysCond {
public:
  //!nowrap!+
  static void set_epsext(float x) {epsext = x;
				   kappasq = hueck/epsext;};
  static void set_solrad(float x) {solrad = x;};
  static void set_sterln(float x) {sterln = x;};
  static void set_ionicstr(float x) {ionicstr = x;
				     hueck=8*pi*conconv*ionicstr/kBolt/T;
				     kappasq = hueck/epsext;};
  static void set_T(float x) {T = x;
			      ln10kT = ln10 * kBolt * T;
			      hueck=8*pi*conconv*ionicstr/kBolt/T;
			      kappasq = hueck/epsext;};
  static void set_kBolt(float x) {kBolt = x;
				  ln10kT = ln10 * kBolt * T;
				  hueck=8*pi*conconv*ionicstr/kBolt/T;
				  kappasq = hueck/epsext;};
  static void set_conconv(float x) {conconv = x;
				    hueck=8*pi*conconv*ionicstr/kBolt/T;
				    kappasq = hueck/epsext;};
  static void set_econv(float x) {econv = x;};

// FIXME!  Bohr radius and proton charge are only important because of
// AmstFitDens.  What other constants need to change when these change?
  static void set_bohr_radius(float x) {bohr_radius = x;}
  static void set_proton_charge(float x) {proton_charge = x;}

  static float get_epsext() {return epsext;};
  static float get_solrad() {return solrad;};
  static float get_sterln() {return sterln;};
  static float get_ionicstr() {return ionicstr;};
  static float get_T() {return T;};
  static double get_kBolt() {return kBolt;};
  static double get_conconv() {return conconv;};
  static double get_econv() {return econv;};
  static double get_hueck() {return hueck;};
  static double get_kappasq() {return kappasq;};
  static double get_ln10kT() {return ln10kT;};
  static double get_bohr_radius() {return bohr_radius;}
  static double get_proton_charge() {return proton_charge;}
  static void print();
  static void print(ostream& os);
  //!nowrap!-

private:
  static float epsext;
  static float solrad;
  static float sterln;
  static float ionicstr;
  static float T;
  static double kBolt;
  static double conconv;
  static double econv;
  static double hueck;
  static double kappasq;
  static double ln10kT;
  static double bohr_radius;
  static double proton_charge;
};

ostream& operator <<(ostream& os, PhysCond& p);

#if SWIGPP_LITERAL_INCLUDE

%addmethods PhysCond {
//
// Kind of a trick for Python to create an "instance" of _p_PhysCond
// A dummy constructor just returns a non-zero address to fool swig.
// It will never be used, I promise!
// Some bit of python code creates the instance (see below).
  PhysCond() {return (PhysCond *) 0x01;}
//
// Wrap all the private members as attributes.
// I think it provides a nicer interface than the get/set functions.
// These ones can be set:
  float epsext;
  float solrad;
  float sterln;
  float ionicstr;
  float T;
  float kBolt;
  float conconv;
  float econv;
  float bohr_radius;
  float proton_charge;
// These are read-only:
  float hueck;
  float kappasq;
  float ln10kT;
// And, of course, the print method can't be named print or Python gets
// confused, so name it write.
  void write() {PhysCond::print();}
};

// Good place to wrap up the global Blab level
%inline %{
class Blab {
public:
  Blab() : level(0) {}
  int get_level() {return level;}
  void set_level(int l) {level=l;}
private:
  int level;
};
%}

%addmethods Blab {
  int level;
};

%wrapper %{

// Get/Set functions for Blab level
void Blab_level_set(Blab *self, int l){
  switch (l) {
    case 0:
      blab1pt = &cnull;
      blab2pt = &cnull;
      blab3pt = &cnull;
      self->set_level(l);
      break;
    case 1:
      blab1pt = &cout;
      blab2pt = &cnull;
      blab3pt = &cnull;
      self->set_level(l);
      break;
    case 2:
      blab1pt = &cout;
      blab2pt = &cout;
      blab3pt = &cnull;
      self->set_level(l);
      break;
    case 3:
      blab1pt = &cout;
      blab2pt = &cout;
      blab3pt = &cout;
      self->set_level(l);
      break;
    default:
      _SWIG_exception(SWIG_ValueError, "The blab level must be between 0 and 3");
      break;
    }
}
int Blab_level_get(Blab *self) {return self->get_level();}

// Provide the get/set functions for the aboved wrapped attributes
void PhysCond_epsext_set(PhysCond *self, float x) {PhysCond::set_epsext(x);}
void PhysCond_solrad_set(PhysCond *self, float x) {PhysCond::set_solrad(x);}
void PhysCond_sterln_set(PhysCond *self, float x) {PhysCond::set_sterln(x);}
void PhysCond_ionicstr_set(PhysCond *self, float x) {PhysCond::set_ionicstr(x);}
void PhysCond_T_set(PhysCond *self, float x) {PhysCond::set_T(x);}
void PhysCond_kBolt_set(PhysCond *self, float x) {PhysCond::set_kBolt(x);}
void PhysCond_conconv_set(PhysCond *self, float x) {PhysCond::set_conconv(x);}
void PhysCond_econv_set(PhysCond *self, float x) {PhysCond::set_econv(x);}
void PhysCond_bohr_radius_set(PhysCond *self, float x) {PhysCond::set_bohr_radius(x);}
void PhysCond_proton_charge_set(PhysCond *self, float x) {PhysCond::set_proton_charge(x);}
void PhysCond_hueck_set(PhysCond *self, float x) {
  _SWIG_exception(SWIG_ValueError, "The hueck attribute is read-only");
}
void PhysCond_kappasq_set(PhysCond *self, float x) {
  _SWIG_exception(SWIG_ValueError, "The kappasq attribute is read-only");
}
void PhysCond_ln10kT_set(PhysCond *self, float x) {
  _SWIG_exception(SWIG_ValueError, "The ln10kT attribute is read-only");
}

float PhysCond_epsext_get(PhysCond *self) {return PhysCond::get_epsext();}
float PhysCond_solrad_get(PhysCond *self) {return PhysCond::get_solrad();}
float PhysCond_sterln_get(PhysCond *self) {return PhysCond::get_sterln();}
float PhysCond_ionicstr_get(PhysCond *self) {return PhysCond::get_ionicstr();}
float PhysCond_T_get(PhysCond *self) {return PhysCond::get_T();}
float PhysCond_kBolt_get(PhysCond *self) {return PhysCond::get_kBolt();}
float PhysCond_conconv_get(PhysCond *self) {return PhysCond::get_conconv();}
float PhysCond_econv_get(PhysCond *self) {return PhysCond::get_econv();}
float PhysCond_bohr_radius_get(PhysCond *self) {return PhysCond::get_bohr_radius();}
float PhysCond_proton_charge_get(PhysCond *self) {return PhysCond::get_proton_charge();}
float PhysCond_hueck_get(PhysCond *self) {return PhysCond::get_hueck();}
float PhysCond_kappasq_get(PhysCond *self) {return PhysCond::get_kappasq();}
float PhysCond_ln10kT_get(PhysCond *self) {return PhysCond::get_ln10kT();}

%}

// Include the python code to create the PhysCond instance
%pragma(python) include="PhysCondInstance.py"

#endif // SWIGPP_LITERAL_INCLUDE

#endif

// PhysCond.h ends here
