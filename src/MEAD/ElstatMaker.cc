#include "MEAD/ElstatMaker.h"
#include "MEAD/ElstatPot.h"
#include "MEAD/DielectricEnvironment.h"
#include "MEAD/ElectrolyteEnvironment.h"
#include "MEAD/ChargeDist.h"

#include <iostream>

DCEsignature::DCEsignature(const DielectricEnvironment_lett* dpt,
			   const ChargeDist_lett* cpt,
			   const ElectrolyteEnvironment_lett* ept)
{
  diel_type = &typeid(*dpt);
  charge_type = &typeid(*cpt);
  ely_type = &typeid(*ept);
}

DCEsignature::DCEsignature(const type_info& d,
			   const type_info& c,
			   const type_info& e)
{
  diel_type = &d;
  charge_type = &c;
  ely_type = &e;
}


bool DCEsignature::operator==(const DCEsignature& a) const
{
  return (*diel_type == *a.diel_type && *charge_type == *a.charge_type
	  && *ely_type == *a.ely_type);
}

void DCEsignature::tell() const
{
  cout << "My DielectricEnvironment is a " << diel_type->name()
    << ",\nmy ElectrolyteEnvironment is a " << ely_type->name()
      << "\nand my ChargeDist is a " << charge_type->name()
	<< endl;
}



ElstatPot_lett* ElstatMaker::maker(DielectricEnvironment_lett* dept,
				   ChargeDist_lett* cdpt,
				   ElectrolyteEnvironment_lett* eept)
{
  DCEsignature argsig(dept, cdpt, eept);
  const DCEsignature argsig_copy = argsig;
  ElstatMaker* p=0;
  for (p=list; p; p=p->next) {
    if (p->operator==(argsig)) return p->derived_maker(dept, cdpt, eept);
  }
  // BEGIN change one argument
  argsig.ely_type = &typeid(ElectrolyteEnvironment_lett);
  for (p=list; p; p=p->next) {
    if (p->operator==(argsig)) return p->derived_maker(dept, cdpt, eept);
  }
  argsig = argsig_copy;
  argsig.charge_type = &typeid(ChargeDist_lett);
  for (p=list; p; p=p->next) {
    if (p->operator==(argsig)) return p->derived_maker(dept, cdpt, eept);
  }
  // END change one argument
  //
  // BEGIN change two arguments
  argsig = argsig_copy;
  argsig.charge_type = &typeid(ChargeDist_lett);
  argsig.ely_type = &typeid(ElectrolyteEnvironment_lett);
  for (p=list; p; p=p->next) {
    if (p->operator==(argsig)) return p->derived_maker(dept, cdpt, eept);
  }
  argsig = argsig_copy;
  argsig.diel_type = &typeid(DielectricEnvironment_lett);
  argsig.ely_type = &typeid(ElectrolyteEnvironment_lett);
  for (p=list; p; p=p->next) {
    if (p->operator==(argsig)) return p->derived_maker(dept, cdpt, eept);
  }
  argsig = argsig_copy;
  argsig.diel_type = &typeid(DielectricEnvironment_lett);
  argsig.charge_type = &typeid(ChargeDist_lett);
  for (p=list; p; p=p->next) {
    if (p->operator==(argsig)) return p->derived_maker(dept, cdpt, eept);
  }
  // END change two arguments

  cerr << "ERROR: ElstatMaker::maker failed.  Returning null pointer." << endl;
  return 0;
}

// ElstatMaker.cc ends here
