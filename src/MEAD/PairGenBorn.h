// -*- C++ -*-
#ifndef _PairGenBorn_h
#define _PairGenBorn_h 1

#include "MEAD/ElstatPot.h"
#include "MEAD/Coord.h"
#include "MEAD/AtomID.h"

#include <map>
using std::map;

class PairGenBorn : public AnalyticEP {
public:
  PairGenBorn(class TwoValueDielectricByAtoms*,
	      class AtomChargeSet*, class ElectrolyteEnvironment_lett*);
  virtual void solve() ;
  virtual float value(Coord) const ;
  virtual float value(const class AtomID&) const ;
  virtual Coord field(Coord x) const ;
  virtual Coord displacement(Coord x) const ;
  float solvation_energy();
private:
  class TwoValueDielectricByAtoms* tvdba;
  class AtomChargeSet* acs;
  class ElectrolyteEnvironment_lett* ue;
  map<AtomID, float> potmap;
};

float operator* (const AtomChargeSet& c, const PairGenBorn& e);
inline float operator* (const PairGenBorn& e, const AtomChargeSet& c)
{ return operator*(c,e); }

#endif

// PairGenBorn.h ends here
