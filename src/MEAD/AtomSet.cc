/* An AtomSet is an associative array of Atoms with some I/O functions
    Copyright (c) 1993--1995 by Donald Bashford.

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


$Id: AtomSet.cc,v 1.30 2008/11/06 20:38:50 bashford Exp $
*/
#include "MEAD/AtomSet.h"
#include "MEAD/globals.h"
#include <fstream>

typedef std::list<Atom>::iterator alitr;
typedef std::list<Atom>::const_iterator calitr;

AtomSet::AtomSet(const AtomSet& a)
: map<AtomID,Atom>(a), name(a.name)
{
}


AtomSet::AtomSet()
: name("")
{
}
AtomSet::AtomSet (const string& nm)
: name(nm)
{

}

AtomSet::AtomSet(const list<Atom>& alin)
{
  for (calitr i=alin.begin(); i != alin.end(); ++i) {
    insert(*i);
  }
}

AtomSet& AtomSet::operator=(const AtomSet& a)
{
  map<AtomID,Atom>::operator=(a);
  name = a.name;
  return *this;
}


/* FIXME
   read is a temporary hack for testing purposes.
   Probably should redo this whole mess with flex.
*/
//#include <string.h>
#include <stdio.h>
#include <cstring>
using std::strncmp;
void AtomSet::read (const string& filename_string)
{
  blab2 << "AtomSet::read: starting\n";
  char linebuf[501];
  char wbuf1[100], wbuf2[100], wbuf3[100];
  int n, nfields;
  const char *filename = filename_string.c_str();

  ifstream atin (filename);
  if (!atin)
    ::error ("AtomSet::read: Cannot open input file, ", filename, "\n");

  clear();
  int atcount=0;
  nfields = 0;

  while (atin.get (linebuf, 500, '\n')) {
    // Read a line

    char junk;

    if (atin.get(junk) && junk != '\n')
      ::error ("readmol: The following line was too long:\n", linebuf, "\n");

    // Check for the word "ATOM" at start of line
    if (!(strncmp (linebuf, "ATOM", 4)==0
	  || strncmp (linebuf, "HETATM", 6)==0)) continue;

    // Read the info for a new atom from the line
    ++atcount;
    Atom at;
    // For backward compatibility, the optional chainid is read last.
    // Bergsma 5/2/01
    wbuf3[0] = '\0';
    n = sscanf (linebuf, "%*s %*d %s %s %d %f %f %f %f %f %s",
	    wbuf1, wbuf2, &at.resnum, &at.coord.x, &at.coord.y, &at.coord.z,
	    &at.charge, &at.rad, wbuf3);
    nfields = (n > nfields) ? n : nfields;
    at.atname = wbuf1;
    at.resname = wbuf2;
    at.chainid = wbuf3;

    AtomID key(at);
    if (contains(key)) {
      // Bergsma changed this to be a fatal error 5/2/01.
      ::error("AtomSet::read: attempt to insert a duplicate atom\n");
/*
 *    cerr << "WARNING: AtomSet::read: Input file, " << filename
 *	<< ",\ncontains duplicate entries for the atom, " << key
 *	  << ".\nLast one will override previous" << endl;
 *    (*this)[key] = at;
 */
    }
    else {
      pair<iterator,bool> p = insert(at);
      if (!p.second) ::error("AtomSet::read: unexpected insertion failure\n");
    }
  }
  if (atcount != size()) {
    cerr << "WARNING AtomSet::read: Map length = " << size()
	 << ", while atom count = " << atcount << endl;
  }
 if (nfields == 9)
   blab2 << "AtomSet::read: records contain chainid fields\n";

  blab1 << "AtomSet::read: finished with " << atcount << " atoms" << endl;
}

#include <vector>
using std::vector;

void AtomSet::build_from_vectors (vector<string> rnames, vector<string> anames,
                          vector<int> rnums, vector<Coord> coords,
                          vector<float> radii, vector<float> charges)
{
  int n = rnames.size();
  for (int i=0; i<n; ++i) {
    Atom at;
    at.atname = anames[i];
    at.resname = rnames[i];
    at.resnum = rnums[i];
    at.chainid = "";		// Don't worry about chainid here
    at.coord = coords[i];
    at.rad = radii[i];
    at.charge = charges[i];
    AtomID key(at);
    if (contains(key)) {
      // FIXME!  This should be a fatal error!
      cerr << "WARNING: AtomSet::build_from_arrays: "
	<< "arrays contain duplicate entries for the atom, " << key
	  << ".\nLast one will override previous" << endl;
      (*this)[key] = at;
    }
    else {
      pair<iterator,bool> p = insert(at);
      if (!p.second) ::error("AtomSet::build_from_arrays: unexpected insertion failure\n");
    }
  }
  if (n != size()) {
    cerr << "WARNING AtomSet::build_from_arrays: Map length = " << size()
	 << ", array size parameter = " << n << endl;
  }
}


void AtomSet::read ()
{
  string filename = name + DOT + "pqr";
  read(filename);
}

Coord AtomSet::geom_cent () const
{
  Coord geom_cent(0.0, 0.0, 0.0);
  for (const_iterator i = begin(); i!=end(); ++i) {
    geom_cent += i->second.coord;
  }
  geom_cent /= float(size());
  return geom_cent;
}

void AtomSet::adjust_atoms_absolute(const class DeltaAtomSet& delta)
{
  const_iterator b;
  // loop over atoms in the "from" set
  for (b = delta.from.begin(); b != delta.from.end(); ++b) {
    AtomID k = (*b).first;
    const_iterator ito = delta.to.find(k);
    iterator ithis = find(k);
    if (ito != delta.to.end()) { // if atom is in the "to" set...
      if (ithis != end()) {// ... it should also be in *this...
	// ... so alter it.
	(*ithis).second = (*ito).second;
      }
      else { // It's odd if it isn't in *this, but not a catastrophe.
	cerr << "WARNING AtomSet::adjust_atoms_absolute:\n   "
	     << "Atom in both \"from\" and \"to\" sets is not in *this set"
	     << "...\n   Adding it to this" << endl;
	insert(b->second);
	// FIXME! the above puts new atoms at the end of the list.
	// Would be better to but them with others of same residue.
      }

    }
    else if (ithis != end()) { // Atom not in "to" set but is in *this
      // so remove this atom
      erase(ithis);
    }
    else { // An atom in "from" is neither in "to" nor *this
      // senseless, but not a catastrophe.  Warn and do nothing.
	cerr << "WARNING AtomSet::adjust_atoms_absolute:\n   "
	     << "Atom in  \"from\", but neither \"to\" nor *this set"
	     << "...\n   Ignoring." << endl;
    }
  } // (end of for loop over "from" atoms)

  for (b = delta.to.begin(); b != delta.to.end(); ++b) {
    AtomID k = (*b).first;
    const_iterator ifrom = delta.from.find(k);
    if (ifrom == end()) {// if not already taken care of above,
      if (contains(k)) {// if we already have it, alter it
	(*this)[k] = b->second;
      }
      else {// otherwise add it.
	insert(b->second);
      }
    }
  }
}

void AtomSet::set_coords_to(const class AtomSet& a)
{
  for(const_iterator b = a.begin(); b != a.end(); ++b) {
    iterator ithis = find(b->first);
    if (ithis != end())
      ithis->second.coord = b->second.coord;
    else {
      cerr << "ERROR: AtomSet::set_coords_to the atom, " << (*b).first
	<< ",\nwhich was specified to have its coord altered,\n"
	  << "does not occur in the AtomSet, " << name << endl;
      ::error();
    }
  }
}

/*
int DeltaAtomSet::differs()
{
  for (Pix i = from.first(); i; from.next(i)) {
    AtomID k = from.key(i);
    Pix j = to.seek(k);
    if (j) {
      Atom atfrom = from.contents(i);
      Atom atto = to.contents(j);
// FIXME!  To much inner knowledge of Atom need here.
// Atom should have a more reasonable == oparator.
      if (!(atfrom == atto && atfrom.coord == atto.coord
	  && atfrom.charge == atto.charge && atfrom.rad == atto.rad))
	return 1;
    }
    else
      return 1;
  }
  for (Pix j = to.first(); j; to.next(j)) {
    AtomID k = to.key(j);
    if (!from.contains(k))
      return 1;
  }
  return 0;
}
*/

// FIXME?  This function might not be so efficient since it makes a full
// pass through each list.  Is there an STL aglorithm that could be used?
// Also, the scheme here is similar to that in rads_differ, etc, so there
// could be a single function that gets called with different predicates.

int DeltaAtomSet::coords_differ()
{
  for (AtomSet::const_iterator i = from.begin(); i!=from.end(); ++i) {
    const Atom& atfrom = i->second;
    AtomID key = i->first;
    if (!to.contains(key)) { // Not found, sets differ
      return 1;
    }
    else { // equiv atom is found, compare coords
      const Atom& atto = to[key];
      if (!(atfrom == atto && atfrom.coord == atto.coord))
	return 1;
    }
  }

  for (AtomSet::const_iterator j = to.begin(); j!=to.end(); ++j) {
    if (!from.contains(j->first)) { // Not found, sets differ
      return 1;
    }
  }
  return 0;
}

/*
int DeltaAtomSet::rads_differ()
{
  for (Pix i = from.first(); i; from.next(i)) {
    AtomID k = from.key(i);
    Pix j = to.seek(k);
    if (j) {
      Atom atfrom = from.contents(i);
      Atom atto = to.contents(j);
// FIXME!  To much inner knowledge of Atom need here.
// Atom should have a more reasonable == oparator.
      if (!(atfrom == atto && atfrom.rad == atto.rad))
	return 1;
    }
    else
      return 1;
  }
  for (Pix j = to.first(); j; to.next(j)) {
    AtomID k = to.key(j);
    if (!from.contains(k))
      return 1;
  }
  return 0;
}
*/

// AtomSet.cc ends here
