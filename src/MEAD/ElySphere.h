// -*- C++ -*-
#ifndef ElySphere_h
#define ElySphere_h 1

#include "MEAD/ElectrolyteEnvironment.h"

//!wrap!
class ElySphere : public ElectrolyteEnvironment_lett {
public:
  //!nowrap!
  ElySphere() : ionic_str(0.0) {}
  ElySphere(float ionstr, Coord center, float r);
  ~ElySphere() {}
  virtual float ionic_strength() const {return ionic_str;}
  //!nowrap!
  virtual ElyCubeRep* get_cuberep(const CubeLatSpec& cls) const ;
  float get_radius() const {return radius;}
  Coord get_center() const {return center;}

private:
  Coord center;
  float radius;
  float ionic_str;
};

#endif

// ElySphere.h ends here
