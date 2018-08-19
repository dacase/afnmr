#include "MEAD/globals.h"
#include "MEAD/FDElstatMaker.h"
#include "MEAD/ElstatPot.h"
#include "MEAD/DielectricEnvironment.h"
#include "MEAD/ElectrolyteEnvironment.h"
#include "MEAD/ChargeDist.h"

#include <iostream>


ElstatPot_lett* FDElstatMaker::maker(FinDiffMethod fdm,
				     DielectricEnvironment_lett* dept,
				     ChargeDist_lett* cdpt,
				     ElectrolyteEnvironment_lett* eept)
{
  DCEsignature argsig(dept, cdpt, eept);
  blab3 << "Entered FDElstatMaker::maker with argsig:" << endl;
  if (blab3pt == &cout) argsig.tell();
  FDElstatMaker* p=0;
  for (p=list; p; p=p->next) {
    if (p->operator==(argsig)) return p->derived_maker(fdm, dept, cdpt, eept);
  }

  argsig.charge_type = &typeid(ChargeDist_lett);
  for (p=list; p; p=p->next) {
    if (p->operator==(argsig)) return p->derived_maker(fdm, dept, cdpt, eept);
  }
  argsig.diel_type = &typeid(DielectricEnvironment_lett);
  for (p=list; p; p=p->next) {
    if (p->operator==(argsig)) return p->derived_maker(fdm, dept, cdpt, eept);
  }
  argsig.ely_type = &typeid(ElectrolyteEnvironment_lett);
  for (p=list; p; p=p->next) {
    if (p->operator==(argsig)) return p->derived_maker(fdm, dept, cdpt, eept);
  }

  cerr << "ERROR: FDElstatMaker::maker failed.  Returning null pointer."
    << endl;
  return 0;
}

// FDElstatMaker.cc ends here
