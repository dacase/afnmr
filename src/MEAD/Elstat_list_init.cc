#include "MEAD/ElstatPot.h"
#include "MEAD/ElstatMaker.h"
#include "MEAD/FDElstatMaker.h"
#include "MEAD/AnalyMaker.h"
#include "MEAD/DielectricEnvironment.h"
#include "MEAD/ElectrolyteEnvironment.h"
#include "MEAD/ChargeDist.h"

ElstatMaker* ElstatMaker::list = 0;
FDElstatMaker* FDElstatMaker::list = 0;
AnalyMaker* AnalyMaker::list = 0;

#include "MEAD/UniformDielectric.h"
#include "MEAD/Debye.h"

#include "MEAD/DielectricSphere.h"
#include "MEAD/AnalySphere.h"

#include "MEAD/DielectricSlab.h"
#include "MEAD/AnalySlab.h"

#include "MEAD/DielByAtoms.h"
#include "MEAD/FinDiffElstatPot.h"
#include "MEAD/ElectrolyteByAtoms.h"

#include "MEAD/DielMembAtoms.h"


class TwoValDielMembAtoms_Maker : public DerivedAnalyMaker<AnalySlab,
TwoValueDielMembAtoms, ChargeDist_lett, ElectrolyteEnvironment_lett> {
public:
  AnalySlab* make_it(TwoValueDielMembAtoms* der_dept,
		     ChargeDist_lett* der_cdpt,
		     ElectrolyteEnvironment_lett* der_eept) const
      {
	return new AnalySlab(der_dept->get_slab(), der_cdpt, der_eept);
      }
};

class ThreeValueDielMembAtomsAtoms_Maker : public DerivedAnalyMaker<AnalySlab,
ThreeValueDielMembAtomsAtoms, ChargeDist_lett, ElectrolyteEnvironment_lett> {
public:
  AnalySlab* make_it(ThreeValueDielMembAtomsAtoms* der_dept,
		     ChargeDist_lett* der_cdpt,
		     ElectrolyteEnvironment_lett* der_eept) const
      {
	return new AnalySlab(der_dept->get_slab(), der_cdpt, der_eept);
      }
};

class DefaultAnalyMaker : public DerivedAnalyMaker<Debye,
DielectricEnvironment_lett, ChargeDist_lett, ElectrolyteEnvironment_lett> {
public:
  Debye* make_it(DielectricEnvironment_lett* der_dept,
		 ChargeDist_lett* der_cdpt,
		 ElectrolyteEnvironment_lett* der_eept) const
    {
      return new Debye(new UniformDielectric(der_dept->epsext_value()),
		       der_cdpt, der_eept);
    }
};


// By defualt, initialization of the Maker lists uses global
// objects with ctors, but if DO_INITS_IN_FUNC is defined, it
// becomes necessary to call an initialization function before
// using MEAD (ugh!)

#ifndef DO_INITS_IN_FUNC

// For a Uniform Dielectric and non-specialized Charge And Electrolyte,
// the ElstatPot_lett will be a Debye and the analytic approx will be a Debye

DerivedElstatMaker<Debye,
   UniformDielectric, ChargeDist_lett, ElectrolyteEnvironment_lett>
debye_maker_ini;

// For a Dielectric Sphere and non-specialized Charge And Electrolyte,
// the ElstatPot_lett will be an AnalySphere.

/*DerivedElstatMaker<AnalySphere,
   DielectricSphere, ChargeDist_lett, ElectrolyteEnvironment_lett>
   AnalySphere_maker_ini; */
//DL
DerivedElstatMaker<AnalySphere,
   DielectricSphere, ChargeDist_lett, ElySphere>
   AnalySphere_maker_ini;

// For a Dielectric Slab and non-specialized Charge And Electrolyte,
// the ElstatPot_lett will be an AnalySlab.

DerivedElstatMaker<AnalySlab,
   DielectricSlab, ChargeDist_lett, ElectrolyteEnvironment_lett>
AnalySlab_maker_ini;

// For a TwoValueDielectricByAtoms and non-specialized Charge And Electrolyte
// the ElstatPot_lett will be a FinDiffElstatPot.

DerivedElstatMaker<FinDiffElstatPot,
  TwoValueDielectricByAtoms, ChargeDist_lett, ElectrolyteEnvironment_lett>
NoFDM_FinDiff_maker_ini;

// This default applies to anything that includes a FinDiffMaker.


DerivedFDElstatMaker<FinDiffElstatPot,
   DielectricEnvironment_lett, ChargeDist_lett, ElectrolyteEnvironment_lett>
FinDiff_maker_ini;

// ------------- List of AnalyMakers -------------

TwoValDielMembAtoms_Maker TwoValMemb_maker_ini;

ThreeValueDielMembAtomsAtoms_Maker ThreeValueDielMembAtomsAtoms_maker_ini;

// For cases not otherwise specified the analytic approx will be a Debye.
// BUG FIXME!  This assumes last one defined ends up as the default;
// but if definitions are spread over several files, no guarantee who is last.

DefaultAnalyMaker default_approx_maker_ini;


#else

// By defualt this branch of the DO_INITS_IN_FUNC is inactive

void do_elstat_inits()
{
  cout << "do_elstat_inits starting" << endl;

  // This is ugly!  It does a bunch of news, and leaves the pointers
  // dangling!   The reason is that we want the ctors invoked, since
  // they do the work of building the global list, but we don't want
  // the objects to be destroyed when this functions exits.

  new DerivedElstatMaker<Debye,
    UniformDielectric, ChargeDist_lett, ElectrolyteEnvironment_lett>();

// For a Dielectric Sphere and non-specialized Charge And Electrolyte,
// the ElstatPot_lett will be an AnalySphere.

/*DerivedElstatMaker<AnalySphere,
   DielectricSphere, ChargeDist_lett, ElectrolyteEnvironment_lett>
   AnalySphere_maker_ini; */
//DL
    new DerivedElstatMaker<AnalySphere,
    DielectricSphere, ChargeDist_lett, ElySphere>();

  // For a Dielectric Slab and non-specialized Charge And Electrolyte,
  // the ElstatPot_lett will be an AnalySlab.


    new DerivedElstatMaker<AnalySlab,
    DielectricSlab, ChargeDist_lett, ElectrolyteEnvironment_lett>();

// For a TwoValueDielectricByAtoms and non-specialized Charge And Electrolyte
// the ElstatPot_lett will be a FinDiffElstatPot.


    new DerivedElstatMaker<FinDiffElstatPot,
    TwoValueDielectricByAtoms, ChargeDist_lett,
    ElectrolyteEnvironment_lett>();

  // This default applies to anything that includes a FinDiffMaker.



    new DerivedFDElstatMaker<FinDiffElstatPot,
    DielectricEnvironment_lett, ChargeDist_lett,
    ElectrolyteEnvironment_lett>();

  // ------------- List of AnalyMakers -------------

  new TwoValDielMembAtoms_Maker();
  new ThreeValueDielMembAtomsAtoms_Maker();

  // For cases not otherwise specified the analytic approx will be a Debye.
  // BUG FIXME!  This assumes last one defined ends up as the default;
  // but if definitions are spread over several files, no guarantee who is last.

  new DefaultAnalyMaker();

  cout << "do_elstat_inits finishing" << endl;
}

#endif

// Elstat_list_init.cc ends here
