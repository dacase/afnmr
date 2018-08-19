// -*- C++ -*-
#ifndef _ManyPointCharge_h
#define _ManyPointCharge_h 1
#include "MEAD/ChargeDist.h"
#include <list>
using std::list;

#include "MEAD/OnePointCharge.h"

struct point_to_point {
  inline PointCharge operator() (const PointCharge& p) {return p;}
};

// Necessary for swigpp.el to process the base class arguments properly
// and for swig typemap conversions
typedef std::list<PointCharge> list_PointCharge;

//!wrap!
class ManyPointCharge : public ChargeDist_lett, public list_PointCharge {
public:
  ManyPointCharge();
  ManyPointCharge(const list_PointCharge& lpc);
  ManyPointCharge(int num_chgs, float *ch, Coord *crd);
  ~ManyPointCharge() {}
  //!nowrap!
  virtual ChargeCubeRep* get_cuberep(const CubeLatSpec&,
				     bool warn_outside=false);
  virtual float total_charge() const ;
  virtual int has_charges() const ;
  virtual float vacuum_coulomb(const Coord& c) const ;
  virtual size_t number_points() const { return size(); }
  // publicize the const versions of begin and end
  //!nowrap!+
  ChargeDist::const_iterator pc_begin() const
    {return make_PointCharge_const_iterator(list<PointCharge>::begin());}
  ChargeDist::const_iterator pc_end() const
    {return make_PointCharge_const_iterator(list<PointCharge>::end());}


  virtual bool dispatch_addChargeDist (const ChargeDist_lett& c,
				       ChargeDist &result) const ;
  ChargeDist addChargeDist (const ManyPointCharge&) const ;
  ChargeDist addChargeDist (const OnePointCharge&) const ;
  //!nowrap!-

private:
  int has_chgs;
  float tot_chg;

};

// Templated addChargeDist's for adding any two ChargeDist types
// that return a new ManyPointCharge
// template<class T> ManyPointCharge * addChargeDist (const T& cd, const OnePointCharge& opc);
// template<class T> ManyPointCharge * addChargeDist (const T& cd, const ManyPointCharge& mpc);
//
// Unfortunately, GNU does not instantiate template functions automatically,
// given only the declarations, and the export command is not working, so...
// Need to define these template functions here in this .h file as
// opposed to the .cc file.
// bergsma - 5/24/01
//
template<class T> ManyPointCharge * addChargeDist (const T& cd, const OnePointCharge& opc)
{
  Coord c_coord = opc.get_coord();
  list<PointCharge> newch;
  bool not_found = true;
  for(typename T::const_iterator i=cd.begin(); i!=cd.end(); ++i) {
  // First try to add the charge to a PointCharge with matching coord.
    const PointCharge& pc = convert_to_PointCharge(*i);
    if (not_found && c_coord == pc.coord) {
      newch.push_back(PointCharge(c_coord, pc.charge + opc.get_charge()));
      not_found = false;
    }
    else
      newch.push_back(pc);
  }

  if (not_found) // No matching coord, so it's a new PointCharge
    newch.push_back(opc);

  return new ManyPointCharge(newch);
}

template<class T> ManyPointCharge * addChargeDist (const T& cd, const ManyPointCharge& mpc)
{
  list<PointCharge> newch(mpc.begin(), mpc.end());
  for (typename T::const_iterator j = cd.begin(); j!=cd.end(); ++j) {
    const PointCharge& pc = convert_to_PointCharge(*j);
    // First try to add the charge to a PointCharge with matching coord.
    std::list<PointCharge>::iterator i=newch.begin();
    for ( ; i != newch.end(); ++i) {
      if (pc.coord == i->coord) {
        i->charge += pc.charge;
        break;
      }
    }
    if (i == newch.end()) // No matching coord, so it's a new PointCharge
      newch.push_back(pc);
  }
  return new ManyPointCharge(newch);
}

// Added to enable the Python interface
// bergsma 5/22/01
ManyPointCharge * multiplyChargeDist (const ManyPointCharge& mpc, float a);
ManyPointCharge * divideChargeDist (const ManyPointCharge& mpc, float a);

#if SWIGPP_LITERAL_INCLUDE

%addmethods ManyPointCharge {
//
// The operation for multiplying by an ElstatPot_lett derived class.
// The general __rmul__ version is defined in the ElstatPot_lett class.
   float operator* (const ManyPointCharge& mpc, const ElstatPot_lett& ep);
//
// The operation for multiplying by an ElstatPotCombination.
// The __rmul__ version is defined in ElstatPotCombination
   float operator* (const ManyPointCharge& mpc, const ElstatPotCombination& epc);
//
// These unusual operators return a pointer to a newly constructed
// object in part because we don't want an expensive copy operation
// everytime __op__ is invoked, but also because they take advantage of
// the addChargeDist routines that are normally used for constructing
// new letter objects for the ChargeDist envelope.
//
// The ways of scaling ManyPointCharges
//
  ManyPointCharge* __mul__ (float a);
  ManyPointCharge* __rmul__ (float a);
  ManyPointCharge* __div__ (float a);
//
// The ways of adding a ManyPointCharge with all other ChargeDist's.
// No __radd__ is needed since Python sees all these as InstanceType's
//
  ManyPointCharge* __add__ (const AtomChargeSet& acs);
  ManyPointCharge* __add__ (const ManyPointCharge& mpc);
  ManyPointCharge* __add__ (const OnePointCharge& opc);
//
// Return a list of Pointcharge items
//
  list_PointCharge * pointcharges()
  {
    return new list_PointCharge(self->begin(), self->end());
  }
//
// Methods to emulate a Python list
//
// Append an item
  void append(const PointCharge& pc) { list_append(self, pc); }
// Length of list
  int __len__ () { return list___len__(self); }
// Add all items from mpc to end of list
  void extend (const ManyPointCharge& mpc) { list_extend(self, mpc); }
// Number of occurrences of item in list
  int count (const PointCharge& pc) { return list_count(self, pc); }
// Index of first occurrence of item
  int index (const PointCharge& pc) { return list_index(self, pc); }
// Insert item at index.
  void insert (int index, const PointCharge& pc) { list_insert(self, index, pc); }
// Remove first occurrence of item
  void remove (const PointCharge& pc) { list_remove(self, pc); }
// Get item at index
  PointCharge __getitem__ (int index) { return list___getitem__(self, index); }
// Set item at index to value
  void __setitem__ (int index, const PointCharge& value) { list___setitem__(self, index, value); }
// Delete item at index
  void __delitem__ (int index) { list___delitem__(self, index); }

};

%wrapper %{

BEGIN_CPLUSPLUS_SECTION

// Define our special __mul__ functions here
ManyPointCharge* ManyPointCharge___mul__ (const ManyPointCharge *mpc, float a)
{ return ::multiplyChargeDist(*mpc, a); }
ManyPointCharge* ManyPointCharge___rmul__ (const ManyPointCharge *mpc, float a)
{ return ::multiplyChargeDist(*mpc, a); }
ManyPointCharge* ManyPointCharge___div__ (const ManyPointCharge *mpc, float a)
{ return ::divideChargeDist(*mpc, a); }

// Define our special __add__ functions here
ManyPointCharge * ManyPointCharge___add__ (const ManyPointCharge *mpc, const AtomChargeSet& acs)
{ return ::addChargeDist(acs, *mpc); }
ManyPointCharge * ManyPointCharge___add__ (const ManyPointCharge *mpc1, const ManyPointCharge& mpc2)
{ return ::addChargeDist(*mpc1, mpc2); }
ManyPointCharge * ManyPointCharge___add__ (const ManyPointCharge *mpc, const OnePointCharge& opc)
{ return ::addChargeDist(*mpc, opc); }

END_CPLUSPLUS_SECTION

%}

#endif // SWIGPP_LITERAL_INCLUDE

#endif

// ManyPointCharge.h ends here
