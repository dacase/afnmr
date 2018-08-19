#include "MEAD/AnalySlab.h"
#include "MEAD/DielectricSlab.h"
#include "MEAD/ElectrolyteEnvironment.h"
#include "MEAD/ChargeDist.h"
#include "MEAD/PhysCond.h"
#include "MEAD/globals.h"
#include <math.h>


AnalySlab::AnalySlab (DielectricSlab* e, ChargeDist_lett *r,
			  ElectrolyteEnvironment_lett* ely, int mt)
: AnalyticEP(e, r, ely), eps(e), rho(r), electrolyte(ely), maxterm(mt)
{
  if (ely->ionic_strength() != 0) {
    error ("AnalySlab constructor: Slab with non-zero ionic strength\n",
	   "not implemented.  Sorry.\n");
  }
  solved = 0;
  epsslab = eps->epsslab_value();
  epsext = eps->epsext_value();
  zupper = eps->zupper_value();
  zlower = eps->zlower_value();

  zmiddle = (zlower + zupper)/2.0F;
  twob = zupper - zlower;
  b = twob/2.0F;
  qinfac = -(epsext-epsslab)/(epsext+epsslab);
  qoutfac = 2*epsext/(epsext+epsslab);

}

void AnalySlab::solve()
{
  if (solved) {
    std::cerr << "WARNING: AnalySlabe::solve called a second time"
	      << std::endl;
    return;
  }
  solved = 1;
}


float AnalySlab::value(Coord c) const
{
  if (!solved) {
    error ("ERROR: AnalySlab::value:\n",
	   "Setup not done yet, must call solve first\n");
  }

  float z = c.z - zmiddle;
  float phi_tot = 0.0;
  for (ChargeDist::const_iterator iso = rho->pc_begin();
       iso!=rho->pc_end(); ++iso) {
    PointCharge pc = *iso;
    pc.coord.z -= zmiddle;
    float dx = c.x - pc.coord.x;
    float dy = c.y - pc.coord.y;
    float coord_rhosq = dx*dx + dy*dy;
    if (pc.coord.z < b && pc.coord.z > -b) { // Source inside slab
      // VALID IF SOURCE POINT INSIDE
      float phi = 0.0;
      if ((z<0.0?-z:z) < b) {
	float zplus = pc.coord.z;
	float zminus = pc.coord.z;
	float zsq = z-zplus;
	zsq *= zsq;
	float rsq = zsq + coord_rhosq;
	float qin = pc.charge;
	// The Coulomb term
	if (rsq != 0) {   // skip singularities FIXME?
	  phi +=  qin / sqrt (rsq) / epsslab;
	}
	for (int n=1; n<=maxterm; ++n) {
	  float tmp_minus = zminus;
	  zminus = -twob - zplus;
	  zplus = twob - tmp_minus;
	  qin *= qinfac;
	  float zsq = z-zminus;
	  zsq *= zsq;
	  float rsqminus =  zsq + coord_rhosq;
	  zsq = z-zplus;
	  zsq *= zsq;
	  float rsqplus =  zsq + coord_rhosq;
	  phi += qin/sqrt(rsqminus) / epsslab + qin/sqrt(rsqplus) / epsslab;
	}
      }
      else if (z <= -b) {
	phi = 0.0;
	float zplus = pc.coord.z;
	float zminus = pc.coord.z;
	float qin = pc.charge;
	for (int n=1; n<=maxterm; ++n) {
	  float qout = qoutfac * qin;
	  float zsq = z-zplus;
	  zsq *= zsq;
	  float rsqplus =  zsq + coord_rhosq;
	  phi += qout/sqrt(rsqplus) / epsext;
	  float tmp_minus = zminus;
	  zminus = -twob - zplus;
	  zplus = twob - tmp_minus;
	  qin *= qinfac;
	}
      }
      else if (z >= b) {
	phi = 0.0;
	float zplus = pc.coord.z;
	float zminus = pc.coord.z;
	float qin = pc.charge;
	for (int n=1; n<=maxterm; ++n) {
	  float qout = qoutfac * qin;
	  float zsq = z-zminus;
	  zsq *= zsq;
	  float rsqminus =  zsq + coord_rhosq;
	  phi += qout/sqrt(rsqminus) / epsext;
	  float tmp_minus = zminus;
	  zminus = -twob - zplus;
	  zplus = twob - tmp_minus;
	  qin *= qinfac;
	}
      }
      phi_tot += phi;
    }
    else { // Source above or below slab
      float zsource, zfpt;
      if (pc.coord.z <= -b) { //Source below so do a sign swap
	zsource = -pc.coord.z;
	zfpt = -z;
      }
      else { // Source above
	zsource = pc.coord.z;
	zfpt = z;
      }
      float phi = 0.0;
      if ((zfpt<0.0 ? -zfpt : zfpt) < b) {
	float zplus = zsource;
	float zminus = zsource;
	float qin = 2.0F*epsslab / (epsslab + epsext) * pc.charge;
	float dz = zfpt-zplus;
	phi = qin / sqrt(dz*dz + coord_rhosq) /epsslab;
	int nterm = maxterm % 2 ? maxterm+1 : maxterm;
	for (int n=2; n<nterm; ++n) {
	  float tmp_zminus = zminus;
	  zminus = -twob - zplus;
	  zplus = twob - tmp_zminus;
	  qin *= qinfac;
	  float dz;
	  if (n%2) // Only odd terms contribute plus terms
	    dz = zfpt-zplus;
	  else // Only even terms contribute minus terms.
	    dz = zfpt-zminus;
	  phi += qin / sqrt(dz*dz + coord_rhosq) /epsslab;
	}
      }
      else if (zfpt >= b) {
	// zero term.
	float qout = pc.charge;
	float dz = zfpt-zsource;
	if (float rsq = dz*dz + coord_rhosq) // skip singularities FIXME?
	  phi += qout / sqrt(rsq) / epsext;
	// -1 term.
	float zminus = twob - zsource;
	float zplus = zminus;
	float qin = 2.0F*epsslab / (epsslab + epsext) * pc.charge;
	qout *= -qinfac;
	dz = zfpt-zminus;
	phi += qout / sqrt(dz*dz + coord_rhosq) / epsext;
	zminus = -twob - zplus;
	zplus = zsource;
	// Rest of (negative) terms.  Only odd ones contribute.
	// The zplus and zminus come from the -(|n|-1) term.
	for (int n=2; n<=maxterm; ++n) {
	  qout = qoutfac * qin;
	  qin *= qinfac;
	  if (n%2) {
	    float dz = zfpt-zminus;
	    phi += qout / sqrt(dz*dz + coord_rhosq) / epsext;
	  }
	  float tmp_minus = zminus;
	  zminus = -twob - zplus;
	  zplus = twob - tmp_minus;
	}
      }
      else if (zfpt <= -b) {
	phi = 0;

	// -1 term.
	float qin = 2.0F*epsslab / (epsslab + epsext) * pc.charge;
	float zminus = zsource;
	float zplus = zsource;
	// Rest of (negative) terms.  Only even ones contribute.
	// The zplus and zminus come from the -(|n|-1) term.
	for (int n=2; n<=maxterm; ++n) {
	  float qout = qoutfac * qin;
	  qin *= qinfac;
	  if ((n%2) == 0) {
	    float dz = zfpt-zplus;
	    phi += qout / sqrt(dz*dz + coord_rhosq) / epsext;
	  }
	  float tmp_minus = zminus;
	  zminus = -twob - zplus;
	  zplus = twob - tmp_minus;
	}
      }
      phi_tot += phi;
    }
  }
  return phi_tot;
}

Coord AnalySlab::field(Coord c) const
{
  ::error("SORRY, AnalySlab::field not implemented yet.  Exiting\n");
  return Coord(0,0,0);
}

Coord AnalySlab::displacement(Coord c) const
{
  ::error("SORRY, AnalySlab::displacement not implemented yet.  Exiting\n");
  return Coord(0,0,0);
}

// AnalySlab.cc ends here
