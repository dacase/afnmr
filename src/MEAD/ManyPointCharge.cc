#include "MEAD/ManyPointCharge.h"
#include "MEAD/OnePointCharge.h"
#include "MEAD/globals.h"

#include <math.h>


ManyPointCharge::ManyPointCharge()
: tot_chg(0), has_chgs(0)
{}

ManyPointCharge::ManyPointCharge(int num_chgs, float ch[], Coord crd[])
{
  tot_chg = 0;
  has_chgs = 0;
  for (int i=0; i<num_chgs; ++i) {
    push_back(PointCharge(crd[i], ch[i]));
    tot_chg += ch[i];
    has_chgs = has_chgs || (ch[i]!=0.0) ? 1 : 0;
  }
}

ManyPointCharge::ManyPointCharge(const list<PointCharge>& lpc)
  : list<PointCharge>(lpc)
{
  tot_chg = 0;
  has_chgs = 0;
  for (const_iterator i= begin(); i!=end(); ++i) {
    float ch = i->charge;
    tot_chg += ch;
    has_chgs = has_chgs || (ch!=0.0) ? 1 : 0;
  }
}


ChargeCubeRep* ManyPointCharge::get_cuberep(const CubeLatSpec&,
					    bool warn_outside)
{cerr << "ManyPointCharge::get_cuberep not implemented" << endl; return 0;}

/*
 * I changed these to always recompute since charges can
 * been updated through STL list operations, or from Python.
 * Was there some reason to be otherwise? FIXME!
 * bergsma 6/25/01
 */
float ManyPointCharge::total_charge() const
{
  float total_chg = 0.0;
  for (const_iterator i= begin(); i!=end(); ++i) {
    total_chg += i->charge;
  }
  return total_chg;
}
//{ return tot_chg; }

int ManyPointCharge::has_charges() const
{
  for (const_iterator i= begin(); i!=end(); ++i) {
    if (i->charge != 0.0) return 1;
  }
  return 0;
}
//{ return has_chgs; }


float ManyPointCharge::vacuum_coulomb(const Coord& c) const
{
  float pot=0.0;
  for (const_iterator i=begin(); i!=end(); ++i) {
    const PointCharge& pc = *i;
    const Coord r = c-pc.coord;
    const float rsq = r*r;
    // skip singularities FIXME?
    if (!rsq) continue;
    pot += pc.charge / sqrt(r*r);
  }
  return pot;
}

// These now call the addChargeDist template functions
// bergsma - 5/24/01

ChargeDist ManyPointCharge::addChargeDist (const OnePointCharge& opc) const
{
  blab3 << "ManyPointCharge::addChargeDist(const OnePointCharge&) called"
    << endl;

  return ChargeDist(::addChargeDist(*this, opc));
}

ChargeDist ManyPointCharge::addChargeDist (const ManyPointCharge& mpc) const
{
  blab3 << "ManyPointCharge::addChargeDist(const ManyPointCharge&) called"
    << endl;

  return ChargeDist(::addChargeDist(*this, mpc));
}

// These could be used to construct ChargeDist envelopes as is done for
// the addChargeDist routines above.
//
ManyPointCharge * multiplyChargeDist(const ManyPointCharge& mpc, float a)
{
  list<PointCharge> newch(mpc.begin(), mpc.end());
  std::list<PointCharge>::iterator i=newch.begin();
  for ( ; i != newch.end(); ++i) {
    i->charge *= a;
  }
  return new ManyPointCharge(newch); 
}

ManyPointCharge * divideChargeDist(const ManyPointCharge& mpc, float a)
{
  if (a == 0) error("ERROR: ManyPointCharge::divideChargeDist: divide by zero");
  list<PointCharge> newch(mpc.begin(), mpc.end());
  std::list<PointCharge>::iterator i=newch.begin();
  for ( ; i != newch.end(); ++i) {
    i->charge /= a;
  }
  return new ManyPointCharge(newch); 
}

#ifdef EXPLICIT_INSTANTIATION
// Explicit instantiation of templates for addChargeDist
template ManyPointCharge * addChargeDist(const ManyPointCharge&, const OnePointCharge&);
template ManyPointCharge * addChargeDist(const ManyPointCharge&, const ManyPointCharge&);
// Explicit instantiation of templates to be used by dispatch function:
template bool try_addChargeDist(const OnePointCharge&, const ManyPointCharge&,
				const ChargeDist_lett&, ChargeDist&);
template bool try_addChargeDist(const ManyPointCharge&, const ManyPointCharge&,
				const ChargeDist_lett&, ChargeDist&);
#endif

bool ManyPointCharge::dispatch_addChargeDist (const ChargeDist_lett& c,
					      ChargeDist &result) const
{
  blab3 << "ManyPointCharge::dispatch_addChargeDist called" << endl;
  if (try_addChargeDist (OnePointCharge(), *this, c, result) ) ;
  else if (try_addChargeDist (ManyPointCharge(), *this, c, result) ) ;
  else
    return false;
  return true;
}

// ManyPointCharge.cc ends here
