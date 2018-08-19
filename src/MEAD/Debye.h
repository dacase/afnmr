// -*- C++ -*-
#ifndef Debye_h
#define Debye_h

#include "MEAD/ElstatPot.h"

class UniformDielectric;

//!wrap!
class Debye : public AnalyticEP {
public:
  Debye(UniformDielectric*, ChargeDist_lett*, ElectrolyteEnvironment_lett*);
  ~Debye() {}
  virtual void solve() {} // A NO-OP
  virtual float value(Coord c) const ;
  virtual Coord field(Coord x) const ;
  virtual Coord displacement(Coord x) const ;
  float get_kappa() const;
//  virtual float operator * (const ChargeDist& c) const ;
private:
  float dielectric_constant;
  float ionic_strength;
  float kappa;
  UniformDielectric* eps;
  ChargeDist_lett* rho;
  ElectrolyteEnvironment_lett* electrolyte;
};

inline float Debye::get_kappa() const {return kappa;}

#if SWIGPP_LITERAL_INCLUDE

%addmethods Debye {
//
// Multiply with a ChargeDist_lett derived class.
// (The __rmul__ version is defined in the ChargeDist_lett derived classes
//
  float operator* (const Debye& dep, const ChargeDist_lett& cdl);
//
// Scale operations
//
  ElstatPotCombination operator* (Debye& dep, float scale);
  ElstatPotCombination operator* (float scale, Debye& dep);
  ElstatPotCombination operator/ (Debye& dep, float scale);
};

#endif // SWIGPP_LITERAL_INCLUDE


#endif

// Debye.h ends here
