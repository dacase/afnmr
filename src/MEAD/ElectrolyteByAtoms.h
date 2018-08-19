// -*- C++ -*-
#ifndef ElectrolyteByAtoms_h
#define ElectrolyteByAtoms_h 1

#include "MEAD/ElectrolyteEnvironment.h"
#include "MEAD/AtomSet.h"
#include "MEAD/PhysCond.h"

struct IonAccVol;
class SolvAccVol;

//!wrap!
class ElectrolyteByAtoms : public ElectrolyteEnvironment_lett {
public:
  //!nowrap!
  ElectrolyteByAtoms() : ionic_str(0.0) {}
  ElectrolyteByAtoms(const AtomSet& ats,
			     float ionic_strnth = PhysCond::get_ionicstr(),
			     float exclus_radius = PhysCond::get_sterln());
  ElectrolyteByAtoms(const AtomSet& ats, float ionic_strnth, const SolvAccVol& sav);
  ~ElectrolyteByAtoms() {}
  virtual float ionic_strength() const {return ionic_str;}
  //!nowrap!
  virtual ElyCubeRep* get_cuberep(const CubeLatSpec& cls) const ;

private:
  SolvAccVol* sav;
  float ionic_str;
  IonAccVol* ion_acc_vol;
};

#if SWIGPP_LITERAL_INCLUDE

%addmethods ElectrolyteByAtoms {
//
// swigpp.el strips the default values of the float arguments in the
// above constructor. So, wrap this constructor with each default
// argument removed and define these special constructors below
  ElectrolyteByAtoms(const AtomSet& ats);
  ElectrolyteByAtoms(const AtomSet& ats, float ionic_strnth);

};

%wrapper %{

BEGIN_CPLUSPLUS_SECTION

// Provide the wrapped constructors with the default arguments supplied
ElectrolyteByAtoms * new_ElectrolyteByAtoms(const AtomSet& ats)
{
  return new ElectrolyteByAtoms(ats, PhysCond::get_ionicstr(), PhysCond::get_sterln());
}
ElectrolyteByAtoms * new_ElectrolyteByAtoms(const AtomSet& ats, float ionic_strnth)
{
  return new ElectrolyteByAtoms(ats, ionic_strnth, PhysCond::get_sterln());
}

END_CPLUSPLUS_SECTION

%}

#endif // SWIGPP_LITERAL_INCLUDE

#endif

// ElectrolyteByAtoms.h ends here
