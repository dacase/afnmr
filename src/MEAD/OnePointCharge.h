// -*- C++ -*-
#ifndef _OnePointCharge_h
#define _OnePointCharge_h 1
#include "MEAD/ChargeDist.h"

//!wrap!
class OnePointCharge : public ChargeDist_lett {
public:
  OnePointCharge();
  OnePointCharge(float ch, Coord crd);

  ~OnePointCharge() {}
  virtual float vacuum_coulomb(const Coord& c) const ;
  //!nowrap!
  virtual ChargeCubeRep* get_cuberep(const CubeLatSpec&,
				     bool warn_outside=false);
  virtual float total_charge() const { return _pc.charge; }
  virtual int has_charges() const { return ! (_pc.charge==0.0); }
  virtual size_t number_points() const { return 1; }

  float get_charge() const {return _pc.charge;}
  Coord get_coord() const {return _pc.coord;}

  //!nowrap!+
  operator PointCharge () const {return _pc;}

  // Iterator for a container of one, and only one PointCharge
  class const_iterator {
  private:
    const OnePointCharge* _opcptr;
    int _ctr; // zero is "pointing to data", otherwise "off the end"
  public:
    const_iterator(const OnePointCharge* opcptr, int ctr)
      : _opcptr(opcptr), _ctr(ctr) {}
    const PointCharge& operator*()
      {return _opcptr->_pc;} // Who needs range checking!
    const_iterator& operator++() { ++_ctr; return *this;}
    bool operator!=(const const_iterator& x) const
      { // Are same containers pointed to and are _ctr values equivalent
	return (_opcptr != x._opcptr)
	  || ( (_ctr==0 ? 0 : 1) != (x._ctr==0 ? 0 : 1));
      }
  };
  friend class OnePointCharge::const_iterator;
  const_iterator begin() const {return const_iterator(this,0);}
  const_iterator end() const {return const_iterator(this,1);}
  ChargeDist::const_iterator pc_begin() const
      {return make_PointCharge_const_iterator(begin());}
  ChargeDist::const_iterator pc_end() const
      {return make_PointCharge_const_iterator(end());}

  virtual bool dispatch_addChargeDist (const ChargeDist_lett& c,
				       ChargeDist &result) const ;
  ChargeDist addChargeDist (const OnePointCharge& c) const ;
  //!nowrap!-


private:
  PointCharge _pc;
};

inline OnePointCharge::OnePointCharge() : ChargeDist_lett(), _pc() {}
inline OnePointCharge::OnePointCharge(float ch, Coord crd) : ChargeDist_lett(), _pc(ch,crd) {}

// Added to enable the Python interface
// bergsma 5/25/01
OnePointCharge * multiplyChargeDist (const OnePointCharge& opc, float a);
OnePointCharge * divideChargeDist (const OnePointCharge& opc, float a);

#if SWIGPP_LITERAL_INCLUDE

%addmethods OnePointCharge {
//
// These class methods need to be added because their definitions are in the
// class body and need to be explicitly defined for SWIG
//
//  float total_charge() { return PointCharge (*self).charge; }
//  int has_charges() { return ! PointCharge (*self).charge==0.0; }
//  int number_points() { return 1; }

//  float get_charge() { return PointCharge (*self).charge; }
//  Coord get_coord() { return PointCharge (*self).coord; }
//
// Return a copy of it's PointCharge
  PointCharge pointcharge() { return PointCharge (*self); };
//
// The operation for multiplying by an ElstatPot_lett derived class.
// The general __rmul__ version is defined in the ElstatPot_lett class.
  float operator* (const OnePointCharge& opc, const ElstatPot_lett& ep);
//
// The operation for multiplying by an ElstatPotCombination derived class.
// The __ruml__ version is in ElstatPotCombination
  float operator* (const OnePointCharge& opc, const ElstatPotCombination& epc);
//
// The ways of scaling a OnePointCharge
//
  OnePointCharge* __mul__ (float a);
  OnePointCharge* __rmul__ (float a);
  OnePointCharge* __div__ (float a);
//
// The ways of adding OnePointCharge with all other ChargeDist's.
// No __radd__ is needed since Python sees all these as InstanceType's
//
  ManyPointCharge* __add__ (const AtomChargeSet& acs);
  ManyPointCharge* __add__ (const ManyPointCharge& mpc);
// Have to pick one or the other, I suppose ManyPointCharge is best
//  OnePointCharge* __add__ (const OnePointCharge& opc);
  ManyPointCharge* __add__ (const OnePointCharge& opc);
 
};

%wrapper %{

BEGIN_CPLUSPLUS_SECTION

// Define our special __mul__ functions here
OnePointCharge* OnePointCharge___mul__ (const OnePointCharge *opc, float a)
{ return ::multiplyChargeDist(*opc, a); }
OnePointCharge* OnePointCharge___rmul__ (const OnePointCharge *opc, float a)
{ return ::multiplyChargeDist(*opc, a); }
OnePointCharge* OnePointCharge___div__ (const OnePointCharge *opc, float a)
{ return ::divideChargeDist(*opc, a); }

// Define our special __add__ functions here
ManyPointCharge * OnePointCharge___add__ (const OnePointCharge *opc, const AtomChargeSet& acs)
{ return ::addChargeDist(acs, *opc); }
ManyPointCharge * OnePointCharge___add__ (const OnePointCharge *opc, const ManyPointCharge& mpc)
{ return ::addChargeDist(mpc, *opc); }
//
// This is a bit unusual, since it can return either a One or Many.
// We'll always make it a Many!
ManyPointCharge * OnePointCharge___add__ (const OnePointCharge *opc1, const OnePointCharge& opc2)
{
  ChargeDist_lett *cdp = ::addChargeDist(*opc1, opc2);
  int npoints = cdp->number_points();
  if (npoints > 1) {
    // We know this has to be a ManyPointCharge!
    return static_cast<ManyPointCharge *> (cdp);
  }
  else {
    // We know this has to be a OnePointCharge!
    // But for some reason this cast does not always work.
    // OnePointCharge *opc = static_cast<OnePointCharge *> (cdp);
    // float chg = opc->get_charge();
    float chg = opc1->get_charge() + opc2.get_charge();
    Coord crd = opc1->get_coord();
    ManyPointCharge *mpc = new ManyPointCharge(1, &chg, &crd);
    delete cdp;
    return mpc;
  }
}

END_CPLUSPLUS_SECTION

%}

#endif // SWIGPP_LITERAL_INCLUDE

#endif

// OnePointCharge.h ends here
