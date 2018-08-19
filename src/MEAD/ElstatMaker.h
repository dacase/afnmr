// -*- C++ -*-
#ifndef ElstatMaker_h
#define  ElstatMaker_h 1

#include <typeinfo>
using std::type_info;
using std::bad_cast;
using std::bad_typeid;

#include <iostream>
#include "MEAD/globals.h"

class DielectricEnvironment_lett;
class ChargeDist_lett;
class ElectrolyteEnvironment_lett;
class ElstatPot_lett;

class DCEsignature {
public:
  DCEsignature(const DielectricEnvironment_lett* dpt,
	       const ChargeDist_lett* cpt,
	       const ElectrolyteEnvironment_lett* ept);
  DCEsignature(const type_info& d, const type_info& c, const type_info& e);
  bool operator==(const DCEsignature& a) const;
  void tell() const ;
  const type_info* diel_type;
  const type_info* charge_type;
  const type_info* ely_type;
};

class ElstatMaker {
public:
  ElstatMaker(const type_info& d, const type_info& c, const type_info& e)
    : sig(d,c,e)
      {
	//	cout << "Entering DerivedElstatMaker ctor" << endl;
	next = list;
	list = this;
      }
  bool operator==(const DCEsignature& s) {return sig==s;}
  static ElstatPot_lett* maker(DielectricEnvironment_lett*,
			       ChargeDist_lett*,
			       ElectrolyteEnvironment_lett*);
private:
  virtual ElstatPot_lett* derived_maker(DielectricEnvironment_lett*,
				   ChargeDist_lett*,
				   ElectrolyteEnvironment_lett*) const = 0;
  DCEsignature sig;
  ElstatMaker *next;
  static ElstatMaker* list;
};


template<class Elstat_T, class Diel_T, class Charge_T, class Ely_T>
class DerivedElstatMaker : public ElstatMaker {
public:
  DerivedElstatMaker()
    : ElstatMaker(typeid(Diel_T), typeid(Charge_T), typeid(Ely_T)) {}
  ElstatPot_lett* derived_maker(DielectricEnvironment_lett* dept,
			   ChargeDist_lett* cdpt,
			   ElectrolyteEnvironment_lett* eept) const
    {
      Diel_T* der_dept = dynamic_cast<Diel_T*>(dept);
      Charge_T* der_cdpt = dynamic_cast<Charge_T*>(cdpt);
      Ely_T* der_eept = dynamic_cast<Ely_T*>(eept);
      Elstat_T* retval = 0;
      if (der_dept && der_cdpt && der_eept)
	retval = new Elstat_T(der_dept, der_cdpt, der_eept);
      else
	error ("ERROR: DerivedElstatMaker::derived__maker: casts failed");
      return retval;
    }
};


#endif

// ElstatMaker.h ends here
