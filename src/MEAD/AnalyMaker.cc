#include "MEAD/AnalyMaker.h"
#include "MEAD/ElstatPot.h"
#include "MEAD/DielectricEnvironment.h"
#include "MEAD/ElectrolyteEnvironment.h"
#include "MEAD/ChargeDist.h"

#include <iostream>

AnalyticEP* AnalyMaker::maker(DielectricEnvironment_lett* dept,
				   ChargeDist_lett* cdpt,
				   ElectrolyteEnvironment_lett* eept)
{
  DCEsignature argsig(dept, cdpt, eept);
  AnalyMaker* p=0;
  for (p=list; p; p=p->next) {
    if (p->operator==(argsig)) return p->derived_maker(dept, cdpt, eept);
  }
  argsig.ely_type = &typeid(ElectrolyteEnvironment_lett);
  for (p=list; p; p=p->next) {
    if (p->operator==(argsig)) return p->derived_maker(dept, cdpt, eept);
  }
  argsig.charge_type = &typeid(ChargeDist_lett);
  for (p=list; p; p=p->next) {
    if (p->operator==(argsig)) return p->derived_maker(dept, cdpt, eept);
  }
  argsig.diel_type = &typeid(DielectricEnvironment_lett);
  for (p=list; p; p=p->next) {
    if (p->operator==(argsig)) return p->derived_maker(dept, cdpt, eept);
  }

  std::cerr << "ERROR: AnalyMaker::maker failed.  Returning null pointer."
	    << std::endl;
  return 0;
}

// AnalyMaker.cc ends here
