#include "MEAD/DielMembAtoms.h"
#include "MEAD/globals.h"
#include "MEAD/PhysCond.h"
#include <iostream>

TwoValueDielMembAtoms::TwoValueDielMembAtoms(const AtomSet& a, float ein,
					     float zl, float zu,
					     Coord c, float holrad)
: solv_acc(a), epsin(ein), epsext(PhysCond::get_epsext()),
  epsmemb(epsin), hole_center(c), hole_radius(holrad)
{
  membranept = new DielectricSlab(epsin, PhysCond::get_epsext(), zl, zu);
  membrane = DielectricEnvironment(membranept);
}

TwoValueDielMembAtoms::TwoValueDielMembAtoms(const SolvAccVol& sav, float ein,
					     float zl, float zu,
					     Coord c, float holrad)
: solv_acc(sav), epsin(ein), epsext(PhysCond::get_epsext()),
  epsmemb(epsin), hole_center(c), hole_radius(holrad)
{
  membranept = new DielectricSlab(epsin, PhysCond::get_epsext(), zl, zu);
  membrane = DielectricEnvironment(membranept);
}


DielCubeRep TwoValueDielMembAtoms::get_cuberep(const CubeLatSpec& cls)
{
  blab3 << "TwoValueDielMembAtoms::get_cuberep called" << endl;

  if (solv_acc.check_is_calculated()) {
    blab2 << "TwoValueDielMembAtoms::get_cuberep using existing "
      << "analytical representation to make a cubic representation" << endl;
  }
  else {
    blab2 << "TwoValueDielMembAtoms::get_cuberep Analytical "
	  << "representation, so call SolvAccVol::anal_calc to calculate it "
	  << endl;
    solv_acc.anal_calc();
  }

  float spacing = cls.get_spacing();
  Coord grid_center_in_space = cls.get_center();
  int grid_dim = cls.get_grid_dim();
  int nsq = grid_dim*grid_dim;
  int ncube = grid_dim*grid_dim*grid_dim;
  float grlen = (float) (grid_dim - 1);
  float halfgrlen = grlen/2;
  Coord grid_center_in_grid (halfgrlen, halfgrlen, halfgrlen);
  float zlbou = (membranept->zlower_value() - grid_center_in_space.z)/spacing
    + grid_center_in_grid.z;
  float zubou = (membranept->zupper_value() - grid_center_in_space.z)/spacing
    + grid_center_in_grid.z;
  Coord holecen_in_grid = (hole_center - grid_center_in_space)/spacing
    + grid_center_in_grid;
  float holerad_in_grid = hole_radius/spacing;
  float holesq = holerad_in_grid*holerad_in_grid;


  AccTag* acc_array = new AccTag[ncube];
  solv_acc.calc_cuberep(cls, acc_array);
  blab3 << "Back in TwoValueDielMembAtoms::get_cuberep.\n"
    << "Using acc_array with modifications to suit membrane..." << endl;

  DielCubeRep dcr(cls);

  for (int k=0; k<grid_dim; ++k) {
    float fk = (float) k;
    if (fk > zlbou && fk < zubou) {
      for (int i = 0; i<grid_dim; ++i) {
	int iterm = i*nsq;
	float xr = ((float) i) - holecen_in_grid.x;
	float xrsq = xr*xr;
	for (int j=0; j<grid_dim; ++j) {
	  int h = iterm + j*grid_dim + k;
// FIXME?  This scheme is not logical for the case epsin==epsext
// but in the TwoValue case it should not matter.  What about ThreeValue?
	  if (acc_array[h] == exterior) {
	    float yr = ((float) j) - holecen_in_grid.y;
	    float rsq = xrsq + yr*yr;
	    if (rsq > holesq)
	      dcr[h] = epsmemb;
	    else
	      dcr[h] = epsext;
	  }
	  else if (acc_array[h] == interior)
	    dcr[h] = epsin;
	  else {
	    cerr << "acc_array[" << h << "] has unexpected value, "
	      << acc_array[h] << endl;
	    ::error("ERROR in ");
	  }
	}
      }
    }
    else { // We are not in the membrane bounds
      for (int i = 0; i<grid_dim; ++i) {
	int iterm = i*nsq;
	for (int j=0; j<grid_dim; ++j) {
	  int h = iterm + j*grid_dim + k;
	  if (acc_array[h] == interior)
	    dcr[h] = epsin;
	  else if (acc_array[h] == exterior)
	    dcr[h] = epsext;
	  else {
	    cerr << "acc_array[" << h << "] has unexpected value, "
	      << acc_array[h] << endl;
	    ::error("ERROR in ");
	  }
	}
      }
    }
  }

  delete [] acc_array;
  dcr.declare_defined();
  blab3 << "Exit from TwoValueDielMembAtoms::get_cuberep"
    << endl;
  return dcr;
}

ThreeValueDielMembAtomsAtoms::ThreeValueDielMembAtomsAtoms
(const AtomSet& a1, float solrad1,
 const AtomSet& a2, float solrad2,
 float zl, float zu, Coord c, float holrad,
 float ein1, float ein2, float eext)
  : solv_acc1(a1, solrad1), solv_acc2(a2, solrad2),
    hole_center(c), hole_radius(holrad),
    epsin1(ein1), epsin2(ein2), epsext(eext)
{
  membranept = new DielectricSlab(epsin2, epsext, zl, zu);
  membrane = DielectricEnvironment(membranept);
}


DielCubeRep 
ThreeValueDielMembAtomsAtoms::get_cuberep(const CubeLatSpec& cls)
{
  blab3 << "ThreeValueDielMembAtomsAtoms::get_cuberep called" << endl;

  if (solv_acc1.check_is_calculated()) {
    blab2 << "Using existing analytical representation to make a cubic\n"
	  << "representation of solv_acc1" << endl;
  }
  else {
    blab2 << "No analytical representation for solvacc1, so call anal_calc\n"
	  << "to calculate it " << endl;
    solv_acc1.anal_calc();
  }
  if (solv_acc1.check_is_calculated()) {
    blab2 << "Using existing analytical representation to make a cubic\n"
	  << "representation of solv_acc2" << endl;
  }
  else {
    blab2 << "No analytical representation for solvacc2, so call anal_calc\n"
	  << "to calculate it " << endl;
    solv_acc2.anal_calc();
  }

  float spacing = cls.get_spacing();
  Coord grid_center_in_space = cls.get_center();
  int grid_dim = cls.get_grid_dim();
  int nsq = grid_dim*grid_dim;
  int ncube = grid_dim*grid_dim*grid_dim;
  float grlen = (float) (grid_dim - 1);
  float halfgrlen = grlen/2;
  Coord grid_center_in_grid (halfgrlen, halfgrlen, halfgrlen);
  float zlbou = (membranept->zlower_value() - grid_center_in_space.z)/spacing
    + grid_center_in_grid.z;
  float zubou = (membranept->zupper_value() - grid_center_in_space.z)/spacing
    + grid_center_in_grid.z;
  Coord holecen_in_grid = (hole_center - grid_center_in_space)/spacing
    + grid_center_in_grid;
  float holerad_in_grid = hole_radius/spacing;
  float holesq = holerad_in_grid*holerad_in_grid;


  AccTag* acc_array1 = new AccTag[ncube];
  AccTag* acc_array2 = new AccTag[ncube];

  solv_acc1.calc_cuberep(cls, acc_array1);
  solv_acc2.calc_cuberep(cls, acc_array2);
  blab3 << "Back in ThreeValueDielMembAtomsAtoms::get_cuberep.\n"
	<< "Using acc_array{1,2} with modifications to suit membrane..." << endl;

  DielCubeRep dcr(cls);

  for (int k=0; k<grid_dim; ++k) {
    float fk = (float) k;
    if (fk > zlbou && fk < zubou) {
      for (int i = 0; i<grid_dim; ++i) {
	int iterm = i*nsq;
	float xr = ((float) i) - holecen_in_grid.x;
	float xrsq = xr*xr;
	for (int j=0; j<grid_dim; ++j) {
	  int h = iterm + j*grid_dim + k;

	  if (acc_array1[h] == interior)
	    dcr[h] = epsin1;
	  else if (acc_array2[h] == interior)
	    dcr[h] = epsin2;
	  else {
	    // since we are within the membrane z-bounds,
	    // we can only be epsext if we are in the hole.
	    float yr = ((float) j) - holecen_in_grid.y;
	    float rsq = xrsq + yr*yr;
	    if (rsq > holesq)
	      dcr[h] = epsin2;
	    else
	      dcr[h] = epsext;
	  }
	}
      }
    }
    else { // We are not in the membrane bounds
      for (int i = 0; i<grid_dim; ++i) {
	int iterm = i*nsq;
	for (int j=0; j<grid_dim; ++j) {
	  int h = iterm + j*grid_dim + k;

	  if (acc_array1[h] == interior)
	    dcr[h] = epsin1;
	  else if (acc_array2[h] == interior)
	    dcr[h] = epsin2;
	  else 
	    dcr[h] = epsext;
	}
      }
    }
  }

  delete [] acc_array1;
  delete [] acc_array2;
  dcr.declare_defined();
  blab3 << "Exit from ThreeValueDielMembAtomsAtoms::get_cuberep"
    << endl;
  return dcr;
}

// DielMembAtoms.cc ends here
