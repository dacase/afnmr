// -*- C++ -*-
#ifndef AnalyMaker_h
#define  AnalyMaker_h 1

#include <typeinfo>
using std::type_info;
using std::bad_cast;
using std::bad_typeid;

#include <iostream>

#include "MEAD/ElstatMaker.h"

class DielectricEnvironment_lett;
class ChargeDist_lett;
class ElectrolyteEnvironment_lett;
class AnalyticEP;

class AnalyMaker {
public:
  AnalyMaker(const type_info& d, const type_info& c, const type_info& e)
    : sig(d,c,e)
      {
	//	cout << "Entering DerivedAnalyMaker ctor" << endl;
	next = list;
	list = this;
      }
  bool operator==(const DCEsignature& s) {return sig==s;}
  static AnalyticEP* maker(DielectricEnvironment_lett*,
			       ChargeDist_lett*,
			       ElectrolyteEnvironment_lett*);
private:
  virtual AnalyticEP* derived_maker(DielectricEnvironment_lett*,
				   ChargeDist_lett*,
				   ElectrolyteEnvironment_lett*) const = 0;
  DCEsignature sig;
  AnalyMaker *next;
  static AnalyMaker* list;
};


// Unlike DerivedElstatMaker, this template makes an abstract class.
// The user must provide its make_it function, in the simplest case this is
//    {new Analy_T(der_dept, der_cdpt, der_eept);}
// but approximations will generally need to be built by hand.

template<class Analy_T, class Diel_T, class Charge_T, class Ely_T>
class DerivedAnalyMaker : public AnalyMaker {
public:
  DerivedAnalyMaker()
    : AnalyMaker(typeid(Diel_T), typeid(Charge_T), typeid(Ely_T)) {}
  AnalyticEP* derived_maker(DielectricEnvironment_lett* dept,
			   ChargeDist_lett* cdpt,
			   ElectrolyteEnvironment_lett* eept) const
    {
      Diel_T* der_dept = dynamic_cast<Diel_T*>(dept);
      Charge_T* der_cdpt = dynamic_cast<Charge_T*>(cdpt);
      Ely_T* der_eept = dynamic_cast<Ely_T*>(eept);
      Analy_T* retval = 0;
      if (der_dept && der_cdpt && der_eept)
	retval = make_it(der_dept, der_cdpt, der_eept);
      else
	error("ERROR: DerivedAnalyMaker::derived__maker: casts failed");
      return retval;
    }
  virtual Analy_T* make_it(Diel_T* der_dept, Charge_T* der_cdpt,
			   Ely_T* der_eept) const
#ifdef _UNICOS
			     {}  // Unicos CC bug prevents making it abstract.
#else
                              = 0;
#endif
};


#endif

// AnalyMaker.h ends here
