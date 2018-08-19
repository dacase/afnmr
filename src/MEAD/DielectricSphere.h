// -*- C++ -*-
#ifndef DielectricSphere_h
#define DielectricSphere_h 1

#include "MEAD/DielectricEnvironment.h"

//!wrap!
class DielectricSphere : public DielectricEnvironment_lett {
public:
  DielectricSphere (float ein, float eext, float rad, Coord cntr);
  ~DielectricSphere () {}
  //!nowrap!
  virtual DielCubeRep get_cuberep(const CubeLatSpec& clsp);
  virtual float epsext_value() const {return epsext;}
  float radius_value() const {return radius;}
  float epsin_value() const {return epsin;}
  Coord get_center() const {return center;}
private:
  float epsin, epsext;
  float radius;
  Coord center;
};

#endif

// DielectricSphere.h ends here
