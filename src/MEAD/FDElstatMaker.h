// -*- C++ -*-
#ifndef FDElstatMaker_h
#define  FDElstatMaker_h 1

#include <typeinfo>
using std::type_info;
using std::bad_cast;
using std::bad_typeid;

#include <iostream>
#include "MEAD/globals.h"
#include "MEAD/FinDiffMethod.h"
#include "MEAD/ElstatMaker.h"

class FDElstatMaker {
public:
  FDElstatMaker(const type_info& d, const type_info& c, const type_info& e)
    : sig(d,c,e)
      {
	//	cout << "Entering DerivedFDElstatMaker ctor" << endl;
	next = list;
	list = this;
      }
  bool operator==(const DCEsignature& s) {return sig==s;}
  static ElstatPot_lett* maker(FinDiffMethod,
			       DielectricEnvironment_lett*,
			       ChargeDist_lett*,
			       ElectrolyteEnvironment_lett*);
private:
  virtual ElstatPot_lett* derived_maker(FinDiffMethod,
					DielectricEnvironment_lett*,
					ChargeDist_lett*,
					ElectrolyteEnvironment_lett*) const =0;
  DCEsignature sig;
  FDElstatMaker *next;
  static FDElstatMaker* list;
};


template<class Elstat_T, class Diel_T, class Charge_T, class Ely_T>
class DerivedFDElstatMaker : public FDElstatMaker {
public:
  DerivedFDElstatMaker()
    : FDElstatMaker(typeid(Diel_T), typeid(Charge_T), typeid(Ely_T)) {}
  ElstatPot_lett* derived_maker(FinDiffMethod fdm,
				DielectricEnvironment_lett* dept,
				ChargeDist_lett* cdpt,
				ElectrolyteEnvironment_lett* eept) const
    {
      Diel_T* der_dept = dynamic_cast<Diel_T*>(dept);
      Charge_T* der_cdpt = dynamic_cast<Charge_T*>(cdpt);
      Ely_T* der_eept = dynamic_cast<Ely_T*>(eept);
      Elstat_T* retval = 0;
      if (der_dept && der_cdpt && der_eept)
	retval = new Elstat_T(fdm, der_dept, der_cdpt, der_eept);
      else
	error ("ERROR: DerivedFDElstatMaker::derived__maker: casts failed");
      return retval;
    }
};


#endif

// FDElstatMaker.h ends here
