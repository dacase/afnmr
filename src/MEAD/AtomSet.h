/* A set of atoms in an associative array keyed by AtomID -*- C++ -*-

    Copyright (c) 1993--1995 by Donald Bashford

    This source code file is part of the MEAD (Macroscopic
    Electrostatics with Atomic Detail) package of objects and
    programs.  MEAD is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 1, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; see the file COPYING.  If not, write to
    the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA
    02139, USA.

    Donald Bashford can be contacted by electronic mail by the address,
    bashford@scripps.edu, or by paper mail at Department of Molecular
    Biology, The Scripps Research Institute, 10666 North Torrey Pines
    Road, La Jolla, California 92037.

$Id: AtomSet.h,v 2.19 2006/05/23 21:23:11 bashford Exp $
*/
#ifndef _AtomSet_h
#define _AtomSet_h 1

#include <map>
#include <list>
#include <string>
using std::string;
#include "MEAD/Atom.h"
#include "MEAD/Coord.h"
#include "MEAD/AtomID.h"
#include "MEAD/globals.h"
#include <vector>
#include <sstream>

using std::vector;
using std::list;
using std::map;
using std::pair;
using std::ostringstream;

class DeltaAtomSet;

// Needed for swig typemap
typedef std::list<Atom> list_Atom;

//!wrap!
class AtomSet: public map<AtomID, Atom > {
public:
  AtomSet();
  AtomSet (const string& name);
  AtomSet (const list_Atom& lat);
  AtomSet (const AtomSet& a);
  ~AtomSet () {}
  //!nowrap!+
  AtomSet& operator=(const AtomSet& a);
  void build_from_vectors (vector<string> rnames, vector<string> anames,
                          vector<int> rnums, vector<Coord> coords,
                          vector<float> radii, vector<float> charges);
  //!nowrap!-
  void read();
  void read (const string& filename);
  Coord geom_cent() const;
  void set_coords_to(const class AtomSet&);
  //!nowrap!+
  void adjust_atoms_absolute(const DeltaAtomSet&); // See below
  string get_name() const {return name;}

  // some deviations from STL container functionality:

  // associative lookup, but don't create if missing
  const Atom& operator[] (const AtomID&) const;
  //!nowrap!-
  Atom& operator[] (const AtomID&);

  //!nowrap!+
  bool contains(const AtomID&) const; // Existence of element (not STL-like)
  pair<iterator,bool> insert (const Atom&); // easier-to-use insert
  //!nowrap!-
private:
  string name;
};

#if SWIGPP_LITERAL_INCLUDE

%addmethods AtomSet {
// Wrap this private member as a read-only attribute for Python, but
// this is a bit of a hack (see below).
// It'd probably be better to just wrap the get_name() function,
// but for Python it seems more natural to treat this simple type
// as an attribute.

  string name;
//
// Keep this for old times sake
  void build_from_vectors (vector_string& rnames, vector_string& anames,
                          vector_int& rnums, vector_Coord& coords,
                          vector_float& radii, vector_float& charges)
  {
    self->build_from_vectors (rnames, anames, rnums, coords, radii, charges);
  }
//
// Not really dict-like, but it's the only way to add new atoms
  void insert (const Atom& a) { self->insert(a); }
//
// Functions to emulate a Python dict
//
  void clear () { dict_clear(self); }
  AtomSet * copy() { return new AtomSet(*self); }
  int has_key (const AtomID& k) { return dict_has_key(self, k); }
  int __len__ () { return dict___len__(self); }
  void __delitem__ (const AtomID& k) { dict___delitem__(self, k); }
// Adds all elements from ats
  void update (const AtomSet& ats) { dict_update(self, ats); }
// This function will convert the map_AtomID_Atom
// to a list of (AtomID, Atom) pairs, which later is converted
// to a python list of tuples by a typemap
// Leaks memory - FIXME!
//  list_AtomID_Atom_pair * items() { return dict_items(self); }
// This function will return a list of AtomID keys from map_AtomID_Atom,
// which later is converted to a python list of keys by a typemap
  list_AtomID * keys() { return dict_keys(self); }
// This function will return a list of Atom values from map_AtomID_Atom,
// which later is converted to a python list of values by a typemap
  list_AtomPtr * values() { return dict_values(self); }

};

%wrapper %{
//
// Provide the read-only get/set functions for the aboved wrapped attribute
void AtomSet_name_set(AtomSet *self, string *name){
   _SWIG_exception(SWIG_ValueError, "The name attribute is read-only");
}
// This type is returned by reference, so create a global placeholder
// to stash the result to be returned (read-only)
// Note that this works OK from Python since assignment of
// simple objects like numbers and strings effectively creates a copy.
// However, using this global reference certainly is not thread safe!
string _AtomSet_name_readonly;
string * AtomSet_name_get(AtomSet *self){_AtomSet_name_readonly = self->get_name(); return &_AtomSet_name_readonly;}

%}

#endif // SWIGPP_LITERAL_INCLUDE

// associative lookup, but don't create if missing
inline const Atom& AtomSet::operator[] (const AtomID& k) const
{
  const_iterator i = find(k);
  if (i == end()) {
    ostringstream ost;
    ost << "AtomSet::operator[]: key, " << k << " not found.";
    ::error(ost.str());
  }
  return i->second;
}

inline Atom& AtomSet::operator[] (const AtomID& k)
{
  iterator i = find(k);
  if (i == end())
    ::error("AtomSet::operator[]: key not found");
  return i->second;
}

// Existence of element (not STL-like
inline bool AtomSet::contains(const AtomID& k) const
{
  return find(k) != end();
}

// easier to use insert (not quite STL-like)
inline pair<AtomSet::iterator,bool> AtomSet::insert (const Atom& a){
  AtomID k(a);
  if (contains(k)) {
    ::error("AtomSet::insert: attempt to insert a duplicate atom\n");
  }
  return map<AtomID,Atom>::insert(value_type(k, a));
}

// Describes a change in an AtomSet, typically involving just a few
// Atoms.  Atoms occuring in "from" are altered to be like their
// corresponding "to" entry if they also occur in "to" or deleted if
// they do not occur in "to".  Atoms occuring in "to" but not "from"
// will be added to the AtomSet if they are not already there, and
// altered if they are.  It is not an error to put into the "from",
// atoms that do not occur in the AtomSet; that is, specifying the
// deletion of non-existing Atoms is OK.  It is also not an error if
// the coords, charges, etc. of a "from" atom differ from the AtomSet,
// since the AtomSet atom will be deleted or overwritten anyway.

class DeltaAtomSet {
public:
  DeltaAtomSet(const AtomSet& from, const AtomSet& to);
  DeltaAtomSet(const list<Atom>& from, const list<Atom>& to);
  //  int differs();
  int coords_differ();
  //  int rads_differ();
private:
friend class AtomSet;
  AtomSet from;
  AtomSet to;
};

inline DeltaAtomSet::DeltaAtomSet(const AtomSet& f, const AtomSet& t)
 : from(f), to(t)
{}

inline
DeltaAtomSet::DeltaAtomSet (const list<Atom>& f, const list<Atom>& t)
 : from(f), to(t)
{}

#endif

// AtomSet.h ends here
