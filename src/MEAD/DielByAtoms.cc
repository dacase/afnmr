#include "MEAD/DielByAtoms.h"
#include "MEAD/AtomSet.h"
#include "MEAD/SolvAccVol.h"
#include "MEAD/PhysCond.h"


TwoValueDielectricByAtoms::TwoValueDielectricByAtoms (const AtomSet & a,
						      float ein)
  : DielectricEnvironment_lett(), _ats(a), _ats_defined(true),
    solv_acc(a), epsin(ein), epsext(PhysCond::get_epsext())

{
}

TwoValueDielectricByAtoms::TwoValueDielectricByAtoms (const AtomSet & a,
						      float ein, float eext,
						      float rad)
 : DielectricEnvironment_lett(), _ats(a), _ats_defined(true),
    solv_acc(a, rad), epsin(ein), epsext(eext)
{
}

TwoValueDielectricByAtoms::TwoValueDielectricByAtoms (const SolvAccVol & a,
						      float ein)
  : DielectricEnvironment_lett(), _ats(), _ats_defined(false),
    solv_acc(a), epsin(ein), epsext(PhysCond::get_epsext())

{
}

TwoValueDielectricByAtoms::TwoValueDielectricByAtoms (const SolvAccVol & a,
						      float ein, float eext)
 : DielectricEnvironment_lett(), _ats(), _ats_defined(false),
    solv_acc(a), epsin(ein), epsext(eext)
{
}

// For GB/DL's sake
const AtomSet& TwoValueDielectricByAtoms::get_atomset() const
{
  if (!_ats_defined)
    ::error("TwoValueDielectricByAtoms::get_atomset: atom set not defined");
  return _ats;
}

DielCubeRep
TwoValueDielectricByAtoms::get_cuberep (const CubeLatSpec& cls)
{
  blab3 << "TwoValueDielectricByAtoms::get_cuberep entered"
    << endl;

  if (solv_acc.check_is_calculated()) {
    blab2 << "TwoValueDielectricByAtoms::get_cuberep using existing "
      << "analytical representation to make a cubic representation" << endl;
  }
  else {
    solv_acc.anal_calc();
  }
  DielCubeRep dcr(cls);
  int grid_dim = cls.get_grid_dim();
  int ncube = grid_dim*grid_dim*grid_dim;
  AccTag* acc_array = new AccTag[ncube];
  solv_acc.calc_cuberep(cls, acc_array);
  for (int h=0; h<ncube; ++h) {
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
  delete [] acc_array;
  dcr.declare_defined();
  blab3 << "exit from TwoValueDielectricByAtoms::get_cuberep\n" << endl;
  return dcr;

}


TwoValueDielectricByAtoms::~TwoValueDielectricByAtoms()
{
}


ThreeValueDielectricByAtoms::ThreeValueDielectricByAtoms
(const AtomSet& a1, float ein1,
 const AtomSet& a2, float ein2)
 : DielectricEnvironment_lett(),
   solv_acc1(a1), solv_acc2(a2)
{
  epsin1 = ein1;
  epsin2 = ein2;
  epsext = PhysCond::get_epsext();
}

ThreeValueDielectricByAtoms::ThreeValueDielectricByAtoms
(const AtomSet& a1, float ein1, float solrad1,
 const AtomSet& a2, float ein2, float solrad2,
 float eext)
 : DielectricEnvironment_lett(), solv_acc1(a1, solrad1), solv_acc2(a2, solrad2)
{
  epsin1 = ein1;
  epsin2 = ein2;
  epsext = eext;
}

ThreeValueDielectricByAtoms::ThreeValueDielectricByAtoms
(const SolvAccVol& a1, float ein1,
 const SolvAccVol& a2, float ein2)
 : DielectricEnvironment_lett(),
   solv_acc1(a1), solv_acc2(a2)
{
  epsin1 = ein1;
  epsin2 = ein2;
  epsext = PhysCond::get_epsext();
}

ThreeValueDielectricByAtoms::ThreeValueDielectricByAtoms
(const SolvAccVol& a1, float ein1,
 const SolvAccVol& a2, float ein2,
 float eext)
 : DielectricEnvironment_lett(), solv_acc1(a1), solv_acc2(a2)
{
  epsin1 = ein1;
  epsin2 = ein2;
  epsext = eext;
}


DielCubeRep
ThreeValueDielectricByAtoms::get_cuberep (const CubeLatSpec& cls)
{
  blab3 << "ThreeValueDielectricByAtoms::get_cuberep entered" << endl;

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
  DielCubeRep dcr(cls);
  int grid_dim = cls.get_grid_dim();
  int ncube = grid_dim*grid_dim*grid_dim;
  AccTag* acc_array1 = new AccTag[ncube];
  AccTag* acc_array2 = new AccTag[ncube];
  solv_acc1.calc_cuberep(cls, acc_array1);
  solv_acc2.calc_cuberep(cls, acc_array2);
  for (int h=0; h<ncube; ++h) {
    if (acc_array1[h] == interior)
      dcr[h] = epsin1;
    else if (acc_array2[h] == interior)
      dcr[h] = epsin2;
    else if (acc_array1[h] == exterior && acc_array2[h] == exterior)
      dcr[h] = epsext;
    else {
      cerr << "ThreeValueDielectric::get_cuberep: "
	<< "acc_arrays have unexpected value at h=" << h << endl;
      ::error("ERROR");
    }
  }
  delete [] acc_array1;
  delete [] acc_array2;
  dcr.declare_defined();
  blab3 << "exit from TwoValueDielectricByAtoms::get_cuberep\n" << endl;
  return dcr;

}

ThreeValueDielectricByAtoms::~ThreeValueDielectricByAtoms()
{
}

// DielByAtoms.cc ends here
