// This is -*- C++ -*-
/*
$Id: AtomID.h,v 2.11 2004/12/06 17:57:57 bashford Exp $
*/
#ifndef _AtomID_h
#define _AtomID_h 1

#include <string>
using std::string;
#include <iostream>
#include "MEAD/Atom.h"

//!wrap!
class AtomID {
public:
  AtomID ();
  AtomID (int residue_number, const string& atom_name, const string& chainid = "");
  AtomID (const AtomID&);
  AtomID (const Atom&);
  ~AtomID () {}
  //!nowrap!
  AtomID& operator= (const AtomID&);
  bool operator== (const AtomID&) const;
  bool operator!= (const AtomID& a) const;
  bool operator<(const AtomID&) const;
  bool operator>(const AtomID&) const;
  //!nowrap!+
  ostream& print(ostream&);
  int get_resnum() const {return resnum;}
  string get_atname() const {return atname;}
  string get_chainid() const {return chainid;}
  //!nowrap!-
private:
  int resnum;
  string atname;
  string chainid;
};

#if SWIGPP_LITERAL_INCLUDE

%addmethods AtomID {
// Wrap these private members as read-only attributes for Python, but
// this is a bit of a hack (see below).
// It'd probably be better to just wrap the get_member() functions,
// but for Python's use it seems more natural to treat these simple types
// as attributes.
  int resnum;
  string atname;
  string chainid;

// swigpp.el strips the default value of the string argument in the
// above constructor. So, wrap this constructor with the default
// argument removed and...
  AtomID (int residue_number, const string& atom_name);
};

%wrapper %{
// ...provide the wrapped constructor with the default argument supplied
AtomID * new_AtomID(int residue_number, const string& atom_name){
  return new AtomID(residue_number, atom_name, "");
}
// Provide the read-only get/set functions for the aboved wrapped attributes
void AtomID_resnum_set(AtomID *self, int resnum){
   _SWIG_exception(SWIG_ValueError, "The resnum attribute is read-only");
}
void AtomID_atname_set(AtomID *self, string *atname){
   _SWIG_exception(SWIG_ValueError, "The atname attribute is read-only");
}
void AtomID_chainid_set(AtomID *self, string *chainid){
   _SWIG_exception(SWIG_ValueError, "The chainid attribute is read-only");
}
int AtomID_resnum_get(AtomID *self){return self->get_resnum();}
// These strings are returned by reference, so create global placeholders
// to stash the results to be returned (read-only).
// Note that this works fine from Python since assignment of
// simple objects like numbers and strings effectively creates a copy.
// However, using these global references certainly is not thread safe!
string _AtomID_atname_readonly;
string _AtomID_chainid_readonly;
string * AtomID_atname_get(AtomID *self){_AtomID_atname_readonly = self->get_atname(); return &_AtomID_atname_readonly;}
string * AtomID_chainid_get(AtomID *self){_AtomID_chainid_readonly = self->get_chainid(); return &_AtomID_chainid_readonly;}
%}

#endif // SWIGPP_LITERAL_INCLUDE

inline AtomID::AtomID (const AtomID& a)
  : resnum(a.resnum), atname(a.atname), chainid(a.chainid) {}

inline AtomID::AtomID (const Atom& a)
  : resnum(a.resnum), atname(a.atname), chainid(a.chainid) {}

inline bool AtomID::operator== (const AtomID& a) const
{return resnum == a.resnum && atname == a.atname && chainid == a.chainid;}

inline bool AtomID::operator!= (const AtomID& a) const
{return !(*this == a);}

inline AtomID& AtomID::operator= (const AtomID& a)
{resnum = a.resnum; atname = a.atname; chainid = a.chainid;
 return *this;}

inline bool AtomID::operator<(const AtomID& x) const
{
  if (chainid == x.chainid) {
    if (resnum == x.resnum)
      return atname < x.atname;
    else
      return resnum < x.resnum;
    }
  else
    return chainid < x.chainid;
}

// Added by Bergsma 5/2/01
inline bool AtomID::operator>(const AtomID& x) const
{
  if (chainid == x.chainid) {
    if (resnum == x.resnum)
      return atname > x.atname;
    else
      return resnum > x.resnum;
    }
  else
    return chainid > x.chainid;
}

inline ostream& operator<<(ostream& ost, AtomID aid)
{return aid.print(ost);}

#endif

// AtomID.h ends here
