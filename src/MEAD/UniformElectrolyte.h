// -*- C++ -*-
#ifndef UniformElectrolyte_h
#define UniformElectrolyte_h 1

#include "MEAD/ElectrolyteEnvironment.h"

//!wrap!
class UniformElectrolyte : public ElectrolyteEnvironment_lett {
public:
  //!nowrap!
  UniformElectrolyte() : ionic_str(0.0) {}
  UniformElectrolyte(float ionstr);
  ~UniformElectrolyte() {}
  virtual float ionic_strength() const {return ionic_str;}
  //!nowrap!
  virtual ElyCubeRep* get_cuberep(const CubeLatSpec& cls) const ;


private:
  float ionic_str;
};

#endif

// UniformElectrolyte.h ends here
