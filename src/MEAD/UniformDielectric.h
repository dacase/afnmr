// -*- C++ -*-
#ifndef UniformDielectric_h
#define UniformDielectric_h 1

#include "MEAD/DielectricEnvironment.h"

//!wrap!
class UniformDielectric : public DielectricEnvironment_lett {
public:
  UniformDielectric(float e);
  ~UniformDielectric() {}
  //!nowrap!
  virtual DielCubeRep get_cuberep(const CubeLatSpec& cls);
  virtual float epsext_value() const {return eps;}
protected:
  float eps;
};

#endif

// UniformDielectric.h ends here
