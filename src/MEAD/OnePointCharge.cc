#include "MEAD/OnePointCharge.h"
#include "MEAD/ManyPointCharge.h"
#include "MEAD/ChargeCubeRep.h"
#include "MEAD/globals.h"
#include <list>
#include <math.h>

ChargeCubeRep*
OnePointCharge::get_cuberep(const CubeLatSpec& cls, bool warn_outside)
{
  blab2 << "Entered OnePointCharge::get_cuberep" << endl;
  int grid_dim = cls.get_grid_dim();
  int nsq = grid_dim * grid_dim;
  float grlen = (float) (grid_dim - 1);

  Coord topcorner(grlen, grlen, grlen);
  Coord gridcenter = topcorner / 2.0;
  Coord r = (_pc.coord - cls.get_center()) / cls.get_spacing() + gridcenter;

  int i = (int) r.x;
  int j = (int) r.y;
  int k = (int) r.z;
  int ip = i+1;
  int jp = j+1;
  int kp = k+1;

  if (i<0 || ip>=grid_dim || j<0 || jp>=grid_dim || k<0 || kp>=grid_dim) {
    cerr << "WARNING: OnePointCharge::get_cubelat: charge outside lattice"
	 << endl;
    return new SparseChargeCubeRep(cls, 0);
  }
  float dx = r.x - (float) i;
  float dy = r.y - (float) j;
  float dz = r.z - (float) k;

  list< std::pair<int,float> > prelim;

  if (i > 0) {
    if (j > 0) {
      if (k > 0)
	prelim.push_back(std::pair<int,float>(i*nsq+j*grid_dim+k,
					 _pc.charge * (1.0F-dx)*(1.0F-dy)*(1.0F-dz)));
      if (kp < grid_dim-1)
	prelim.push_back(std::pair<int,float>(i*nsq+j*grid_dim+kp,
					 _pc.charge * (1.0F-dx)*(1.0F-dy)*(dz)));
    }
    if (jp < grid_dim-1) {
      if (k > 0)
	prelim.push_back(std::pair<int,float>(i*nsq+jp*grid_dim+k,
					 _pc.charge * (1.0F-dx)*(dy)*(1.0F-dz)));
      if (kp < grid_dim-1)
	prelim.push_back(std::pair<int,float>(i*nsq+jp*grid_dim+kp,
					 _pc.charge * (1.0F-dx)*(dy)*(dz)));
    }
  }
  if (ip < grid_dim-1) {
    if (j > 0) {
      if (k > 0)
	prelim.push_back(std::pair<int,float>(ip*nsq+j*grid_dim+k,
					 _pc.charge * (dx)*(1.0F-dy)*(1.0F-dz)));
      if (kp < grid_dim-1)
	prelim.push_back(std::pair<int,float>(ip*nsq+j*grid_dim+kp,
					 _pc.charge * (dx)*(1.0F-dy)*(dz)));
    }
    if (jp < grid_dim-1) {
      if (k > 0)
	prelim.push_back(std::pair<int,float>(ip*nsq+jp*grid_dim+k,
					 _pc.charge * (dx)*(dy)*(1.0F-dz)));
      if (kp < grid_dim-1)
	prelim.push_back(std::pair<int,float>(ip*nsq+jp*grid_dim+kp,
					 _pc.charge * (dx)*(dy)*(dz)));
    }
  }

  ChargeCubeRep *ccrp = new SparseChargeCubeRep(cls, prelim.size());
  for (std::list< std::pair<int,float> >::const_iterator itr = prelim.begin();
       itr != prelim.end(); ++itr) {
    ccrp->add((*itr).first, (*itr).second);
  }


  blab2 << "Exiting OnePointCharge::get_cuberep" << endl;
  return ccrp;
}


float OnePointCharge::vacuum_coulomb(const Coord& c) const
{
  const Coord r = c-_pc.coord;
  const float rsq = r*r;
  // skip singularities FIXME?
  const float pot = rsq ? _pc.charge / sqrt(rsq) : 0.0F;
  return pot;
}

ChargeDist_lett * addChargeDist (const OnePointCharge& opc1, const OnePointCharge& opc2)
{
  ChargeDist_lett* cdp = 0;
  PointCharge pc1 = opc1;
  PointCharge pc2 = opc2;
  if (pc2.coord == pc1.coord)
    cdp = new OnePointCharge(pc1.charge + pc2.charge, pc1.coord);
  else { // Yuck!  FIXME!
    float ch[2];
    ch[0] = pc1.charge; ch[1] = pc2.charge;
    Coord crd[2];
    crd[0] = pc1.coord; crd[1] = pc2.coord;
    cdp = new ManyPointCharge(2, ch, crd);
  }
  if (!cdp)
    error("ERROR: OnePointCharge::addChargeDist: new fails");
  return (cdp);
}

ChargeDist OnePointCharge::addChargeDist (const OnePointCharge& c) const
{
  blab3 << "OnePointCharge::addChargeDist (const OnePointCharge&) called"
    << endl;

  return ChargeDist(::addChargeDist(*this, c));
}

// These could be used to construct ChargeDist envelopes as is done for
// the addChargeDist routines above.
//
OnePointCharge * multiplyChargeDist(const OnePointCharge& opc, float a)
{
  PointCharge pc = opc;
  pc.charge *= a;
  return new OnePointCharge(pc.charge, pc.coord);
}

OnePointCharge * divideChargeDist(const OnePointCharge& opc, float a)
{
  if (a == 0) error("ERROR: OnePointCharge::divideChargeDist: divide by zero");
  PointCharge pc = opc;
  pc.charge /= a;
  return new OnePointCharge(pc.charge, pc.coord);
}

#ifdef EXPLICIT_INSTANTIATION
// Explicit instantiation of template to be used by dispatch function:
template bool try_addChargeDist(const OnePointCharge&, const OnePointCharge&,
			   const ChargeDist_lett&, ChargeDist&);
#endif

bool OnePointCharge::dispatch_addChargeDist (const ChargeDist_lett& c,
				       ChargeDist &result) const
{
//  blab2 << "OnePointCharge::dispatch_addChargeDist called" << endl;
  if (try_addChargeDist (OnePointCharge(), *this, c, result) ) ;
  else
    return false;
  return true;
}

// OnePointCharge.cc ends here
