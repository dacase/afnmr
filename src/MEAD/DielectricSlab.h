// -*- C++ -*-
#ifndef DielectricSlab_h
#define DielectricSlab_h 1

#include "MEAD/DielectricEnvironment.h"

//!wrap!
class DielectricSlab : public DielectricEnvironment_lett {
public:
  DielectricSlab(float eslab, float eext, float zl, float zu);
  ~DielectricSlab() {}
  //!nowrap!
  virtual DielCubeRep get_cuberep(const CubeLatSpec& clsp);
  float epsslab_value() const {return epsslab;}
  float epsext_value() const {return epsext;}
  float zupper_value() const {return zupper;}
  float zlower_value() const {return zlower;}
private:
  float epsslab, epsext;
  float zupper, zlower;
};

#endif

// DielectricSlab.h ends here
