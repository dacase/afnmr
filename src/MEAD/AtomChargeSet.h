// This is -*- C++ -*-
/* ChargDensity subclass for sets of Atoms

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

$Id: AtomChargeSet.h,v 2.16 2007/05/28 01:26:41 bashford Exp $
*/
#ifndef _AtomChargeSet_h
#define _AtomChargeSet_h 1

#include "MEAD/ChargeDist.h"
#include "MEAD/AtomSet.h"
#include "MEAD/ManyPointCharge.h"
#include "MEAD/OnePointCharge.h"

class OnePointCharge;
class ManyPointCharge;

inline PointCharge
convert_to_PointCharge (const std::pair<const AtomID,Atom>& p)
{
  return PointCharge(p.second.coord, p.second.charge);
}

//!wrap!
class AtomChargeSet : public ChargeDist_lett, public AtomSet {
public:
  AtomChargeSet();
  AtomChargeSet(const AtomSet& atl);
  AtomChargeSet(const AtomChargeSet&);
  ~AtomChargeSet() {}
  //!nowrap!+
  AtomChargeSet& operator=(const AtomChargeSet&);
  virtual ChargeCubeRep* get_cuberep(const CubeLatSpec&,
				     bool warn_outside=false);
  //!nowrap!-
  virtual float total_charge() const ;
  virtual int has_charges() const ;
  virtual float vacuum_coulomb(const Coord& c) const;
  virtual size_t number_points() const { return size(); }
  //!nowrap!+
  virtual ChargeDist::const_iterator pc_begin() const
    {return make_PointCharge_const_iterator(begin());}
  virtual ChargeDist::const_iterator pc_end() const
    {return make_PointCharge_const_iterator(end());}

  //    {return STL_iter_to_CDL_iter_ptr<pair_to_point>(end());}

  // Functions inherited from AtomSet that update might the charges
  // must be redefined so as not to screw up the charge statistics FIXME!

  virtual bool dispatch_addChargeDist (const ChargeDist_lett& c,
				       ChargeDist &result) const ;
  ChargeDist addChargeDist (const ManyPointCharge&) const ;
  ChargeDist addChargeDist (const OnePointCharge&) const ;
  ChargeDist addChargeDist (const AtomChargeSet&) const ;
  //!nowrap!-


// Some resonable substitutes for these are needed.  FIXME!
//  virtual float debyeMultiply(const Debye& phi) const;
//  virtual float analySphereMultiply(const AnalySphere& phi) const;
//  virtual float analySlabMultiply(const AnalySlab& phi) const;
//  virtual float finDiffElstatPotMultiply(const FinDiffElstatPot& phi) const;
//  virtual ChargeDist operator + (const ChargeDist& c) const
//    {return c.atomChargeSetAdd(*this);}
//  virtual ChargeDist operator - (const ChargeDist& c) const
//    {return c.atomChargeSetSubtract(*this);}
//  virtual ChargeDist atomChargeSetSubtract(const AtomChargeSet &acs) const;
// #ifndef NO_AMST
//  virtual ChargeDist amstFitDensAdd(const AmstFitDens &acs) const;
//  virtual ChargeDist amstFitDensSubtract(const AmstFitDens &acs) const;
// #endif


private:
  void update_charge_stats() const ;
  class charge_accumulator {
    float charge;
  public:
    charge_accumulator() : charge(0.0) {}
    void operator() (const value_type& p) {charge += p.second.charge;}
    float get_sum() const {return charge;}
    void reset() {charge=0;}
  };
  struct is_charged_t {
    bool operator() (const value_type& p) {return p.second.charge != 0.0;}
  } is_charged;
};

// Added to enable the Python interface
// bergsma 5/21/01
AtomChargeSet * multiplyChargeDist (const AtomChargeSet& acs, float a);
AtomChargeSet * divideChargeDist (const AtomChargeSet& acs, float a);
AtomChargeSet * addChargeDist (const AtomChargeSet& acs1, const AtomChargeSet& acs2);

#if SWIGPP_LITERAL_INCLUDE

%addmethods AtomChargeSet {
//
// The operation for multiplying by an ElstatPot_lett derived class.
// The general __rmul__ version is defined in the ElstatPot_lett class.
   float operator* (const AtomChargeSet& acs, const ElstatPot_lett& ep);
//
// The operation for multiplying by an ElstatPotCombination.
// The __ruml__ version is defined in ElstatPotCombination.
   float operator* (const AtomChargeSet& acs, const ElstatPotCombination& epc);
//
// These unusual operators return a pointer to a newly constructed
// object in part because we don't want an expensive copy operation
// everytime __op__ is invoked, but also because they take advantage of
// the addChargeDist routines that are normally used for constructing
// new letter objects for the ChargeDist envelope.
//
// The ways of scaling AtomChargeSets
  AtomChargeSet* __mul__ (float a);
  AtomChargeSet* __rmul__ (float a);
  AtomChargeSet* __div__ (float a);
// The ways of adding AtomChargeSet with all other ChargeDist's.
// No __radd__ is needed since Python sees all these as InstanceType's
  AtomChargeSet* __add__ (const AtomChargeSet& acs);
  ManyPointCharge* __add__ (const ManyPointCharge& mpc);
  ManyPointCharge* __add__ (const OnePointCharge& opc);
//
// Return a list of Pointcharge items
  list_PointCharge * pointcharges()
  {
    list_PointCharge *lpc = new list_PointCharge(0);
    for (AtomSet::const_iterator ind = self->begin(); ind != self->end(); ++ind)
    {
      const Atom& a = ind->second;
      lpc->push_back(PointCharge(a.coord, a.charge));
    }
    return lpc;
  }
};

%wrapper %{

BEGIN_CPLUSPLUS_SECTION

// Define our special __mul__ functions here
AtomChargeSet * AtomChargeSet___mul__ (const AtomChargeSet *acs, float a)
{ return ::multiplyChargeDist(*acs, a); }
AtomChargeSet * AtomChargeSet___rmul__ (const AtomChargeSet *acs, float a)
{ return ::multiplyChargeDist(*acs, a); }
AtomChargeSet * AtomChargeSet___div__ (const AtomChargeSet *acs, float a)
{ return ::divideChargeDist(*acs, a); }

// Define our special __add__ functions here
AtomChargeSet * AtomChargeSet___add__ (const AtomChargeSet *acs1, const AtomChargeSet& acs2)
{ return ::addChargeDist(*acs1, acs2); }
ManyPointCharge * AtomChargeSet___add__ (const AtomChargeSet *acs, const ManyPointCharge& mpc)
{ return ::addChargeDist(*acs, mpc); }
ManyPointCharge * AtomChargeSet___add__ (const AtomChargeSet *acs, const OnePointCharge& opc)
{ return ::addChargeDist(*acs, opc); }

END_CPLUSPLUS_SECTION

%}

#endif // SWIGPP_LITERAL_INCLUDE

#endif

// AtomChargeSet.h ends here
