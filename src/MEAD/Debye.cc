#include "MEAD/Debye.h"
#include "MEAD/UniformDielectric.h"
#include "MEAD/ElectrolyteEnvironment.h"
#include "MEAD/ChargeDist.h"
#include "MEAD/PhysCond.h"

#include <math.h>


Debye::Debye(UniformDielectric* e, ChargeDist_lett* r,
	     ElectrolyteEnvironment_lett* ely)
: AnalyticEP(e, r, ely), eps(e), rho(r), electrolyte(ely)
{
  blab3 << "Entering Debye ctor." << endl;
  dielectric_constant = e->epsext_value();
  ionic_strength = ely->ionic_strength();
  if (ionic_strength) {
    float kappasq = float(8.0 * pi * PhysCond::get_conconv() * ionic_strength
      / (PhysCond::get_kBolt() * PhysCond::get_T() * PhysCond::get_epsext()));
    kappa = sqrt(kappasq);
  }
  else
    kappa = 0.0;
}


float Debye::value(Coord c) const
{
  if (kappa) {
    float pot = 0;
    for (ChargeDist::const_iterator b = rho->pc_begin();
	 b!=rho->pc_end(); ++b) {
      const PointCharge p = *b;
      const Coord d = c - p.coord;
      const float r = sqrt(d*d);
      if (r != 0.0) // skip singularities FIXME?
	pot += p.charge * exp(-kappa*r) / r / dielectric_constant;
    }
    return pot;
  }
  else
    return rho->vacuum_coulomb(c) / dielectric_constant;
}

Coord Debye::field(Coord c) const
{
  ::error("SORRY, Debye::field not implemented yet.  Exiting\n");
  return Coord(0,0,0);
}

Coord Debye::displacement(Coord c) const
{
  ::error("SORRY, Debye::displacement not implemented yet.  Exiting\n");
  return Coord(0,0,0);
}

/*
float Debye::operator* (const ChargeDist& c) const
{
  return 0.0;
}
*/

// Debye.cc ends here
