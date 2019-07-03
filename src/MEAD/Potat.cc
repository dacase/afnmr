/* Potentials at each atom with storage in files. -*- C++ -*-

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

$Id: Potat.cc,v 2.19 2013/10/29 21:08:15 bashford Exp $
*/
#include <string>
#include <iostream>
using std::ios;
#include "MEAD/Potat.h"
#include "MEAD/AtomSet.h"
#include "MEAD/AtomChargeSet.h"
#include "MEAD/ElstatPot.h"
#include "MEAD/globals.h"

typedef AtomSet::const_iterator ascitr;

Potat::Potat() {
  defined = 0;
}

Potat::Potat(const AtomSet& ats)
:_atset(ats)
{
  defined = 0;
}

void Potat::zero()
{
  for (ascitr i = _atset.begin(); i != _atset.end(); ++i)
    _map[i->first] = 0.0;
}

float Potat::operator*(const AtomChargeSet& acs) {
  if (!defined) {
    ::error("ERROR Potat, AtomChargeSet multiplication invoked for\n",
	    "a Potat with undefined potential values\n");
    return 0;
  }
  float prod = 0;
  for (ascitr b = acs.begin(); b!=acs.end(); ++b) {
    AtomID key = b->first;
    std::map<AtomID, Potval>::const_iterator i = _map.find(key);
    if (i == _map.end()) {
      cerr << "ERROR WARNING: Potat::operator*: The atom, " << key
	<< ",\nfrom the AtomChargeSet does not occur in this Potat" << endl;
    }
    else {
      prod += i->second * b->second.charge;
    }
  }
  return prod;
}

InPotat::InPotat(const AtomSet &as) : Potat(as) {}

int InPotat::read(const string& filename)
{
  blab3 << "InPotat::read has no new implementation, call read_oldstyle"
    << endl;
  return read_oldstyle(filename);
}

int InPotat::read_oldstyle(const string& filename_string)
{
  const char *filename = filename_string.c_str();
  blab3 << "InPotat::read_oldstyle entered" << endl;
  if (_atset.empty()) {
    cerr << "WARNING: InPotat::read_oldstyle: This potat has no atoms defined "
      << ",\nso no potat reading possible" << endl;
    return 0;
  }
  if (filename_string == "") {
    cerr << "WARNING InPotat::read_oldstyle: null filename, "
      << "no potat file read.\n";
    return 0;
  }
  blab2 << "Trying to read a potat from file, " << filename << endl;
  ifstream potatfile(filename, ios::in | ios::binary);
  if (!potatfile.good())
    return 0;
  int nat;
  potatfile.read((char *) &nat, sizeof (int));
  if (!potatfile) {
    cerr << "WARNING: InPotat::read_oldstyle:\n  Some kind of error occurred "
      << "while reading potat info from " << filename << endl;
    return 0;
  }
  if (nat != _atset.size()) {
    cerr << "WARNING: InPotat::read_oldstyle:\n"
      << "The number of atoms read from the file, " << filename	<< ", is "
	<< nat << ";\nbut the number expected is " << _atset.size()
	  << ".  No potat info read"  << endl;
    return 0;
  }
  float *bigarr = new float[nat*6];
  potatfile.read((char *) bigarr, (sizeof (float)) * nat * 6);
  if (!potatfile) {
    cerr << "WARNING: InPotat::read_oldstyle: Some kind of error occurred "
      << "while reading potat info from " << filename << endl;
    delete[] bigarr;
    return 0;
  }
  if (potatfile.get() != EOF) {
    cerr << "WARNING: OutPotat::read_oldstyle:\n"
      << "Some extra data at end of file, " << filename << " ignored." << endl;
  }
  int n = 0;
  for (ascitr b = _atset.begin(); b!=_atset.end(); ++b) {
    bool foundit = false;
    Atom a = b->second;
    if (n>=nat) {
      ::error("WARNING OutPotat::write_oldstyle: too many atoms in atset\n");
      delete[] bigarr;
      return 0;
    }
    float *file_atom_entry = bigarr + n*6;
    Coord filecoor(file_atom_entry[0], file_atom_entry[1], file_atom_entry[2]);
    Coord diff = filecoor - a.coord;
    float sumsq = diff*diff;
    float raddiff = file_atom_entry[3] - a.rad;
    sumsq += raddiff*raddiff;
    if (sumsq > 1.0e-4) { // No easy match, so scan list for it
      for (ascitr c = _atset.begin(); c!=_atset.end(); ++c) {
	a = c->second;
	Coord diff = a.coord - filecoor;
	float sumsq = diff*diff;
	float raddiff = file_atom_entry[3] - a.rad;
	sumsq += raddiff*raddiff;
	if (sumsq < 1.0e-4) {
	  foundit = true;
	  break;
	}
      }
    }
    else {
      foundit = true;
    }
    if (foundit) {
      AtomID key(a);
      if (_map.find(key) == _map.end()) { // Not yet in map. Good.
	_map[key] = file_atom_entry[5];
      }
      else {
	cerr << "WARNING InPotat::read_oldstyle:\n"
	     << "There seems to be either a duplicate atom for "
	     << key << ",\n or superposed atoms in the input." << endl;
	delete[] bigarr;
	return 0;
      }
    }
    else { // Didn't find it
      cerr << "WARNING InPotat::read_oldstyle:\n"
	   <<"   Discrepencies between atom info in file, " << filename
	   << "\n   and atom info expeced." << endl;
      delete[] bigarr;
      return 0;
    }
    ++n;
  }
  delete[] bigarr;
  blab2 << "Successfully read a potat from file, " << filename << endl;
  defined = 1;
  return 1;
}

OutPotat::OutPotat(const AtomSet& as, const ElstatPot& phi) : Potat(as)
{
  for (ascitr b = _atset.begin(); b!=_atset.end(); ++b) {
    const Atom& a = b->second;
    AtomID key = b->first;
    _map[key] = phi.value(a.coord);
  }
  defined = 1;
}

int OutPotat::write(const string& filename) const
{
  blab3 << "OutPotat::write has no new implementation, call write_oldstyle"
    << endl;
  return write_oldstyle(filename);
}

int OutPotat::write_oldstyle(const string& filename) const
{
  blab3 << "OutPotat::write_oldstyle entered" << endl;
  if (filename == "") {
    cerr << "WARNING OutPotat::write_oldstyle: null filename, "
      << "no potat file written\n";
    return 0;
  }
  blab2 << "Trying to write a potat to file, " << filename << endl;
  size_t nat = _atset.size();
  float *bigarr = new float[nat*6];
  size_t n=0;
  for (ascitr b = _atset.begin(); b!=_atset.end(); ++b) {
    const Atom& a = b->second;
    AtomID key = b->first;
    std::map<AtomID, Potval>::const_iterator ipot = _map.find(key);
    if (ipot == _map.end()) {
      cerr << "ERROR OutPotat::write_oldstyle: The atom, " << key
	<< ",\nis supposed to be in this potat, but isn't" << endl;
      delete[] bigarr;
      return 1;
    }
    if (n>=nat)
      error ("OutPotat::write_oldstyle: array overflow (n>=nat)\n");
    bigarr[n*6+0] = a.coord.x;
    bigarr[n*6+1] = a.coord.y;
    bigarr[n*6+2] = a.coord.z;
    bigarr[n*6+3] = a.rad;
    bigarr[n*6+4] = 0; // FIXME! The generating charge (if any) should go here
    bigarr[n*6+5] = ipot->second;
    ++n;
  }
  ofstream potatfile;
  safeopen(potatfile, filename);
  potatfile.write((char *) &nat, sizeof (int));
  potatfile.write((char *) bigarr, (sizeof (float)) * nat * 6);
  delete[] bigarr;
  if (potatfile.good()) {
    blab2 << "Successfully wrote a potat to file, " << filename << endl;
    return 1;
  }
  else {
    cerr << "WARNING: OutPotat::write_oldstyle:\n"
      << "Some kind of problem writing to " << filename << endl;
    return 0;
  }
}

// Potat.cc ends here
