// -*- C++ -*-
#ifndef DielByAtoms_h
#define DielByAtoms_h 1

#include "MEAD/DielectricEnvironment.h"
#include "MEAD/SolvAccVol.h"
#include "MEAD/AtomSet.h"


//!wrap!
class TwoValueDielectricByAtoms : public DielectricEnvironment_lett {
public:
//  TwoValueDielectricByAtoms (const AtomSet&);
  TwoValueDielectricByAtoms (const AtomSet&, float ein);
  TwoValueDielectricByAtoms (const AtomSet&, float ein, float eext, float rad);
  TwoValueDielectricByAtoms (const SolvAccVol&, float ein);
  TwoValueDielectricByAtoms (const SolvAccVol&, float ein, float eext);
  virtual ~TwoValueDielectricByAtoms();
  float epsext_value() const {return epsext;}
  float epsin_value() const {return epsin;}
  //!nowrap!+
  virtual DielCubeRep get_cuberep(const CubeLatSpec& clsp);
  const AtomSet& get_atomset() const; // For GB/DL's sake
  //!nowrap!-
private:
  float epsin, epsext;
  SolvAccVol solv_acc;
  AtomSet _ats;   // For GB/DL's sake
  bool _ats_defined;
};

//!wrap!
class ThreeValueDielectricByAtoms : public DielectricEnvironment_lett {
public:
  ThreeValueDielectricByAtoms(const AtomSet& a1, float ein1, const AtomSet& a2, float ein2);
  ThreeValueDielectricByAtoms(const AtomSet& a1, float ein1, float solrad1, const AtomSet& a2, float ein2, float solrad2, float epsext);
  ThreeValueDielectricByAtoms(const SolvAccVol& a1, float ein1, const SolvAccVol& a2, float ein2);
  ThreeValueDielectricByAtoms(const SolvAccVol& a1, float ein1, const SolvAccVol& a2, float ein2, float epsext);
  virtual ~ThreeValueDielectricByAtoms();
  float epsext_value() const {return epsext;}
  //!nowrap!
  virtual DielCubeRep get_cuberep(const CubeLatSpec& clsp);
private:
  float epsin1, epsin2, epsext;
  SolvAccVol solv_acc1, solv_acc2;
};

#endif

// DielByAtoms.h ends here
