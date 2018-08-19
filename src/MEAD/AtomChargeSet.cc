/* ChargDensity subclass for sets of Atoms

    Copyright (c) 1993--1995 by Donald Bashford

    This source code file is part of the MEAD (Macroscopic
    Electrostatics with Atomic Detail) package of objects and
    programs.  MEAD is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 1, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; see the file COPYING.  If not, write to
    the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA
    02139, USA.

    Donald Bashford can be contacted by electronic mail by the address,
    bashford@scripps.edu, or by paper mail at Department of Molecular
    Biology, The Scripps Research Institute, 10666 North Torrey Pines
    Road, La Jolla, California 92037.

$Id: AtomChargeSet.cc,v 2.16 2007/05/28 01:26:41 bashford Exp $
*/

#include "MEAD/AtomChargeSet.h"

#include "MEAD/OnePointCharge.h"
#include "MEAD/ManyPointCharge.h"
#include "MEAD/globals.h"
#include "MEAD/CubeLatSpec.h"
#include "MEAD/ChargeCubeRep.h"
#include "MEAD/Bigmem.h"
#include <string.h>
#include <math.h>
#include <algorithm>


AtomChargeSet::AtomChargeSet()
 : ChargeDist_lett(), AtomSet()
{
}

AtomChargeSet::AtomChargeSet(const AtomSet& atl)
 : ChargeDist_lett(), AtomSet(atl)
{
}

AtomChargeSet::AtomChargeSet(const AtomChargeSet& acs)
 : ChargeDist_lett(acs), AtomSet(acs)
{
}

AtomChargeSet&
AtomChargeSet::operator=(const AtomChargeSet& acs)
{
  ChargeDist_lett::operator=(acs);
  AtomSet::operator=(acs);
  return *this;
}


float AtomChargeSet::total_charge() const
{
  charge_accumulator acc;
  acc = std::for_each(begin(), end(), acc);
  return acc.get_sum();
}

int AtomChargeSet::has_charges() const
{
  for (const_iterator i = begin(); i != end(); ++i)
    if (i->second.charge != 0) return 1;
  return 0;
}

float AtomChargeSet::vacuum_coulomb(const Coord& c) const
{
  float coul = 0;
  typedef AtomSet::const_iterator Iter;
  for (Iter i = begin(); i != end(); ++i) {
    const Atom& thisatom = i->second;
    if (!thisatom.charge) continue; // Go to next atom if this is uncharged.
    const Coord r = thisatom.coord - c;
    const float dist = sqrt(r*r);
    if (dist != 0.0)  // skip singularities FIXME?
      coul += thisatom.charge/dist;
  }
  return coul;
}

ChargeCubeRep* AtomChargeSet::get_cuberep(const CubeLatSpec& cls,
					  bool warn_outside)
{
  blab2 << "Entered AtomChargeSet::get_cuberep" << endl;
  int grid_dim = cls.get_grid_dim();
  int nsq = grid_dim * grid_dim;
  int ncube = nsq * grid_dim;
  // Allocate a temp array for charges of grid points;
  float *chtmp = (float *) big_rigid_malloc ((sizeof (float)) * ncube);
  memset ((char *) chtmp, 0, ncube * sizeof (float)); // zero the array
  float grlen = (float) (grid_dim - 1);
  Coord topcorner(grlen, grlen, grlen);
  Coord gridcenter = topcorner / 2.0;
  Coord bottomcorner(0.0, 0.0, 0.0);
  int num_charges_out = 0;
  typedef AtomSet::const_iterator Iter;
  for (Iter itr = begin(); itr != end(); ++itr) {
    const Atom& thisatom = itr->second;
    float chartmp = thisatom.charge;
    if (!chartmp) continue; // Go to next atom if this is uncharged.
    Coord thiscoord = thisatom.coord;
    // r will be the atom position scaled and offset to grid units.
    Coord r = (thiscoord - cls.get_center()) / cls.get_spacing()
      + gridcenter;
//    if (r < bottomcorner || r > topcorner) {
//      cout << "WARNING PBEquations::chmake: charge outside of box\n";
//      continue;
//    }
    int i = (int) r.x;
    int j = (int) r.y;
    int k = (int) r.z;
    int ip = i+1;
    int jp = j+1;
    int kp = k+1;
    if (i<0 || ip>=grid_dim || j<0 || jp>=grid_dim || k<0 || kp>=grid_dim) {
      ++num_charges_out;
      continue;
    }
    float dx = r.x - (float) i;
    float dy = r.y - (float) j;
    float dz = r.z - (float) k;
    chtmp[i*nsq+j*grid_dim+k] += chartmp * (1.0F-dx)*(1.0F-dy)*(1.0F-dz);
    chtmp[i*nsq+j*grid_dim+kp] += chartmp * (1.0F-dx)*(1.0F-dy)*(dz);
    chtmp[i*nsq+jp*grid_dim+k] += chartmp * (1.0F-dx)*(dy)*(1.0F-dz);
    chtmp[i*nsq+jp*grid_dim+kp] += chartmp * (1.0F-dx)*(dy)*(dz);
    chtmp[ip*nsq+j*grid_dim+k] += chartmp * (dx)*(1.0F-dy)*(1.0F-dz);
    chtmp[ip*nsq+j*grid_dim+kp] += chartmp * (dx)*(1.0F-dy)*(dz);
    chtmp[ip*nsq+jp*grid_dim+k] += chartmp * (dx)*(dy)*(1.0F-dz);
    chtmp[ip*nsq+jp*grid_dim+kp] += chartmp * (dx)*(dy)*(dz);
  }
  if (num_charges_out && warn_outside)
    cerr << "WARNING: AtomChargeSet::get_cuberep: " << num_charges_out
      << "charges fall outside of the lattice,\n" << cls << endl;

  // chtmp done.  Get ready to make the actual ChargeCubeRep.

  // Note that charges on outer boundary points will be dropped.
  // This is as it should be.

  // FIXME?  But this costs a bit to do the skipping.  Is it worth it?
  // Perhaps it is only really important to skip the -x, x sides of the box,
  // which would be cheap.

  // First count the number of charged grid points.
  int chcount=0;
  int nm1 = grid_dim - 1;
  int i;
  for (i = 1; i<nm1; ++i) {
    int iterm = i*nsq;
    for (int j=1; j<nm1; ++j) {
      int jterm = j*grid_dim;
      for (int k=1; k<nm1; ++k) {
	int h = iterm + jterm + k;
	if (chtmp[h] != 0.0)
	  ++chcount;
      }
    }
  }

  // Create and fill the (Sparse) ChargeCubeRep
  ChargeCubeRep *ccrp = new SparseChargeCubeRep(cls, chcount);
  for (i = 1; i<nm1; ++i) {
    int iterm = i*nsq;
    for (int j=1; j<nm1; ++j) {
      int jterm = j*grid_dim;
      for (int k=1; k<nm1; ++k) {
	int h = iterm + jterm + k;
	if (chtmp[h] != 0.0)
	  ccrp->add(h, chtmp[h]);
      }
    }
  }

  big_rigid_free(chtmp);
  blab2 << "Exiting AtomChargeSet::get_cuberep" << endl;
  return ccrp;
}


/* Some resonable substitutes for these are needed.  FIXME!

float AtomChargeSet::debyeMultiply(const Debye& phi) const
{
  blab3 << "Multiplying{" << "}\ntimes\n{"
    << "}" << endl;
  return 1.0;
}

float AtomChargeSet::analySphereMultiply(const AnalySphere& phi) const
{
  blab3 << "Multiplying{" << "}\ntimes\n{"
    << "}" << endl;
  return 1.0;
}

float AtomChargeSet::analySlabMultiply(const AnalySlab& phi) const
{
  blab3 << "Multiplying{"  << "}\ntimes\n{"
    << "}" << endl;
  return 1.0;
}

float
AtomChargeSet::finDiffElstatPotMultiply(const FinDiffElstatPot& phi) const
{
  blab3 << "Multiplying{" << "}\ntimes\n{"
    << "}" << endl;
  float product = 0;
  for (Pix ind = first(); ind; next(ind)) {
    Atom a = contents(ind);
    float pot = phi.value(a.coord);
    product += pot * a.charge;
  }
  return product;
}

*/

// These now call the addChargeDist template functions.
// bergsma - 5/24/01

ChargeDist AtomChargeSet::addChargeDist(const OnePointCharge &opc) const
{
  blab3 << "AtomCharge::addChargeDist(const OnePointCharge&) called"
    << endl;

  return ChargeDist(::addChargeDist (*this, opc));
}

ChargeDist AtomChargeSet::addChargeDist(const ManyPointCharge &mpc) const
{
  blab3 << "AtomCharge::addChargeDist(const ManyPointCharge&) called"
    << endl;

  return ChargeDist(::addChargeDist (*this, mpc));
}


AtomChargeSet * addChargeDist(const AtomChargeSet &acs1, const AtomChargeSet &acs2)
{
  AtomSet ats = acs1;
  for (AtomChargeSet::const_iterator j = acs2.begin(); j!=acs2.end(); ++j) {
    const Atom& acs2_ar = j->second;
    AtomID acs2_key = j->first;
    if (ats.contains(acs2_key)) {
      Atom& atat = ats[acs2_key];
      if (acs2_ar.coord == atat.coord) {
	// The same atom in the same place, so add the charges
	atat.charge += acs2_ar.charge;
      }
      else { // Same atom, but in a different place
	::error("AtomChargeSet::addChargeDist: ",
		"adding sets containing the same atom with different coords ",
		"is forbidden!");
	// FIXME.  Perhaps the real solution here is to back out and create
	// the result as a ManyPointCharge.  Use exceptions to implement?
      }
    }
    else { // Atom not in current set, so arrange to add it
      ats.insert(acs2_ar);
    }
  }
  return new AtomChargeSet(ats);
}

ChargeDist AtomChargeSet::addChargeDist(const AtomChargeSet &acs) const
{
  blab3 << "AtomChargeSet::addChargeDist(const AtomChargeSet&) called"
	<< endl;

  return ChargeDist(::addChargeDist(*this, acs));
}

// These could be used to construct ChargeDist envelopes as is done for
// the addChargeDist routines above.
//
AtomChargeSet * multiplyChargeDist(const AtomChargeSet& acs, float a)
{
  list<Atom> lat;
  for (AtomChargeSet::const_iterator i = acs.begin(); i != acs.end(); ++i) {
    Atom at = i->second;
    at.charge *= a;
    lat.push_back(at);
  }
  AtomSet ats(lat);
  return new AtomChargeSet(ats);
}

AtomChargeSet * divideChargeDist(const AtomChargeSet& acs, float a)
{
  if (a == 0) error("ERROR: AtomChargeSet::divideChargeDist: divide by zero");
  list<Atom> lat;
  for (AtomChargeSet::const_iterator i = acs.begin(); i != acs.end(); ++i) {
    Atom at = i->second;
    at.charge /= a;
    lat.push_back(at);
  }
  AtomSet ats(lat);
  return new AtomChargeSet(ats);
}

#ifdef EXPLICIT_INSTANTIATION
// Explicit instantiation of templates for addChargeDist
template ManyPointCharge * addChargeDist(const AtomChargeSet&, const OnePointCharge&);
template ManyPointCharge * addChargeDist(const AtomChargeSet&, const ManyPointCharge&);
// Explicit instantiation of templates to be used by dispatch function:
template bool try_addChargeDist(const OnePointCharge&, const AtomChargeSet&,
				const ChargeDist_lett&, ChargeDist&);
template bool try_addChargeDist(const ManyPointCharge&, const AtomChargeSet&,
				const ChargeDist_lett&, ChargeDist&);
template bool try_addChargeDist(const AtomChargeSet&, const AtomChargeSet&,
				const ChargeDist_lett&, ChargeDist&);
#endif

bool AtomChargeSet::dispatch_addChargeDist (const ChargeDist_lett& c,
					    ChargeDist &result) const
{
  blab3 << "ManyPointCharge::dispatch_addChargeDist called" << endl;
  if (try_addChargeDist (OnePointCharge(), *this, c, result) ) ;
  else if (try_addChargeDist (ManyPointCharge(), *this, c, result) ) ;
  else if (try_addChargeDist (AtomChargeSet(), *this, c, result) ) ;
  else
    return false;
  return true;
}



/*
ChargeDist AtomChargeSet::atomChargeSetSubtract(const AtomChargeSet &acs) const
{
  AtomSet at(*(AtomSet *)this);
  for (Pix ind = acs.first(); ind; acs.next(ind)) {
    Atom atm = acs.contents(ind);
    atm.charge *= -1; // since this is subtraction.
    AtomID key = acs.key(ind);
    if (at.contains(key)) {
      Atom old_atm = at[key];
      if (old_atm.coord == atm.coord) {
	atm.charge += old_atm.charge;
	at[key] = atm;
      }
      else {
	cerr << "WARNING: AtomChargeSet::atomChargeSetSubtract:\n"
	  << "Coords do not match" << endl;
      }
    }
    else {
      at[key] = atm;
    }
  }
  ChargeDist cd(at);
  float sum = total_chrg - acs.total_charge();
  if (cd.total_charge() != sum) {
    cerr << "WARNING: AtomChargeSet::atomChargeSetSubtract: "
	 << "Charge totals non-additive (" << sum << " vs. "
	 << cd.total_charge() << ")" << endl;
  }
  return  cd;
}

#ifndef NO_AMST
ChargeDist AtomChargeSet::amstFitDensAdd(const AmstFitDens &acs) const
{::error ("ERROR: AtomChargeSet::amstFitDensAdd not implemented");
 return *(ChargeDist*)this;}

ChargeDist AtomChargeSet::amstFitDensSubtract(const AmstFitDens &acs) const
{::error ("ERROR: AtomChargeSet::amstFitDensSubtract not implemented");
 return *(ChargeDist*)this;}
#endif
*/

// AtomChargeSet.cc ends here
