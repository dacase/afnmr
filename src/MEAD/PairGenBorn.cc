#include "MEAD/PairGenBorn.h"
#include "MEAD/globals.h"
#include "MEAD/AtomChargeSet.h"
#include "MEAD/ElectrolyteEnvironment.h"
#include "MEAD/PhysCond.h"
#include "MEAD/DielByAtoms.h"
#include <iostream>
#include <math.h>


PairGenBorn::PairGenBorn
(TwoValueDielectricByAtoms* tvdbai, AtomChargeSet* acsi,
 ElectrolyteEnvironment_lett* uei)
  : AnalyticEP(tvdbai, acsi, uei), tvdba(tvdbai), acs(acsi), ue(uei)
{
}

#include <vector>
using std::vector;
#include <map>
using std::map;

void PairGenBorn::solve()
{
  blab2 << "PairGenBorn::solve entered" << endl;
  if (tvdba->epsin_value() != 1.0)
    ::error("ERROR PairGenBorn::solve: values of interior diel. const.\n",
	    "other than 1.0 not allowed in this implementation.  Sorry.\n");

  /* define a set of scaling factors for different atom types */
  map<char,float> scalefact_map;
  scalefact_map['N'] = 0.79F;
  scalefact_map['H'] = 0.85F;
  scalefact_map['C'] = 0.72F;
  scalefact_map['O'] = 0.85F;
  scalefact_map['S'] = 0.96F;
  scalefact_map['P'] = 0.86F;
  scalefact_map['F'] = 0.88F;

  const float BOFFSET = 0.09F; // rborn offset. A.O. 02.27.98

  /* Calculate the solvation energy using GB method */

  const AtomSet& dats = tvdba->get_atomset();

  vector<float> rrescaled;
  rrescaled.reserve(dats.size());
  map<AtomID,float> reff;


  /*    get the "effective" Born radii via the approximate pairwise method
	Use Eqs 9-11 of Hawkins, Cramer, Truhlar, J. Phys. Chem. 100:19824
	(1996).  Below, lij and uij and the inverses of the L and U
	parameters defined in this reference.                             */


  /* get physical conditions  */
  const float diel_ext = tvdba->epsext_value();
  float kappa = 0;
  const float ionic_strength = ue->ionic_strength();
  if (ionic_strength) {
    float kappasq =
      8.0 * pi * PhysCond::get_conconv() * ionic_strength
      / (PhysCond::get_kBolt() * PhysCond::get_T() * diel_ext);
    kappa = sqrt(kappasq);
  }
  else
    kappa = 0.0;



  // Loop to get the (empirically) rescaled radii for later use in
  // calculating effective radii.  Also initialize potmap to zeros
  for (AtomSet::const_iterator ind = dats.begin(); ind!=dats.end(); ++ind) {
    if(ind->second.rad < 0.01) { // Check if any rborn == 0.0
	cerr << "WARNING: suspiciously small born radius  " <<
	ind->second.rad << "  of atom " << ind->second.atname
	     << "  in residue " << ind->second.resname << endl;
    }
    std::map<char,float>::const_iterator isf
      = scalefact_map.find(ind->second.atname[0]);
    if (isf != scalefact_map.end())
      rrescaled.push_back(ind->second.rad * isf->second);
    else {
      cerr << "Atom type " << ind->second.atname[0]
	   << " unknown.  Treating as carbon" << endl;
      // next line A.O. 02.27.98
      rrescaled.push_back((ind->second.rad - BOFFSET) * scalefact_map['C']);
    }
    potmap.insert( pair<AtomID,float> (ind->first, 0.0) );
  }

	/* Compute the GB effective radii */

  for (AtomSet::const_iterator ati = dats.begin(); ati!=dats.end(); ++ati) {
    const Coord ci = ati->second.coord;
    const float radi = ati->second.rad - BOFFSET; // A.O. 02.27.98
    float sumi = 1.0F/radi;

    { // Both AtomSet dats, and the rrescaled vector are iterated over
      std::vector<float>::const_iterator rsj= rrescaled.begin();
      for (AtomSet::const_iterator atj = dats.begin();
	   atj!=dats.end(); ++atj, ++rsj) {
	if( ati==atj ) continue;
	const Coord disp = atj->second.coord - ci;
	const float dij = sqrt( disp*disp );

	const float sj =  *rsj;

	if (radi < dij+sj) { // rescaled j not enclosed by i
	  // 1/lij = i radius if no overlap,
	  //          otherwise dist from i to near surf of j
	  const float lij = 1.0F / ( radi > dij-sj ? radi : dij-sj );
	  // 1/uij = dist to far surface of j
	  const float uij = 1.0F / (dij+sj);
	  const float diff_lusq = uij*uij - lij*lij;
	  const float temp3 = lij - uij + 0.25F*dij*diff_lusq
	    + (0.5F/dij)*log(uij/lij)
	    + (0.25F*sj*sj/dij)*(-diff_lusq);
	  sumi -= 0.5F*temp3;
	}
      }
    }
    const float min_rad_allowed = radi + BOFFSET; // A.O. 02.27.98
    const float calc_rad = 1.0F/sumi;
    reff[ati->first] = calc_rad > min_rad_allowed ? calc_rad : min_rad_allowed;
  }

  // Compute the GB potentials
  for (AtomChargeSet::const_iterator chi = acs->begin();
       chi != acs->end(); ++chi) {
    Coord ci = chi->second.coord;
    const float qi = chi->second.charge;
    AtomSet::const_iterator fnd = dats.find(chi->first);
    if (fnd == dats.end() || fnd->second.coord != ci)
      ::error("ERROR: PairGenBorn::solve: charge carrying atom has no\n",
	      "equivalent among the dielectric-defining atoms\n");
    const float reffi = reff[chi->first];
    { // The AtomSet dats, the reff map, and the potmap are all iterated over
      std::map<AtomID,float>::const_iterator rej= reff.begin();
      std::map<AtomID,float>::iterator potj= potmap.begin();
      for (AtomSet::const_iterator atj = dats.begin();
	   atj!=dats.end(); ++atj, ++rej, ++potj) {
	const Coord disp = ci - atj->second.coord;
	const float r2 = disp*disp;
	const float rb2 = reffi * rej->second;
	const float fgb = sqrt(r2 + rb2*exp(-r2/(4.0F*rb2))); // Still notation
	const float dielfac = 1.0F - exp(-(kappa)*fgb*0.73F)/(diel_ext);
	// The "generalized Born" part of V, self-energy term included
	potj->second += -dielfac * qi / fgb;
	// Add back the vacuum potential (no self-energy part)
	if(fnd != atj) { // This much checked probably not needed (FIXME?)
	  if (r2 == 0.0)
	    cerr << "WARNING PairGenBorn::solve: "
		 << "atom pair with zero distance\n"
		 << "will exclude from the vacuum part" << endl;
	  else
	    potj->second += qi/sqrt(r2);
	}
      }
    }
  }
  blab2 << "PairGenBorn::solve finised" << endl;
}




float PairGenBorn::value(Coord c) const
{
  const AtomSet& atset = tvdba->get_atomset();
  for (AtomSet::const_iterator i = atset.begin(); i!=atset.end(); ++i) {
    if (i->second.coord == c) {
      std::map<AtomID,float>::const_iterator k = potmap.find(i->first);
      return k->second;
    }
  }
  cerr << "WARNING PairGenBorn::value(Coord): "
       << " No atom with coord = " << c << " found " << endl;
  return 0;
}

float PairGenBorn::value(const AtomID& aid) const
{
  std::map<AtomID,float>::const_iterator k = potmap.find(aid);
  if (k == potmap.end()) ::error("PairGenBorn::value(AtomID)",
				 "atom not found");
  return k->second;
}

Coord PairGenBorn::field(Coord x) const
{
  ::error("PairGenBorn::field not implemented");
  return Coord();
}

Coord PairGenBorn::displacement(Coord x) const
{
  ::error("PairGenBorn::field not implemented");
  return Coord();
}

float PairGenBorn::solvation_energy()
{
  ::error("PairGenBorn::field not implemented");
  return 0.0;
}


float operator* (const AtomChargeSet& c, const PairGenBorn& e)
{
  float prod = 0;
  for(AtomChargeSet::const_iterator i = c.begin(); i != c.end(); ++i)
    prod += i->second.charge * e.value(i->first);
  return prod;
}

// PairGenBorn.cc ends here
