/* Definition of some global constants, variables and functions
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

$Id: globals.cc,v 1.28 2007/11/26 22:38:35 bashford Exp $
*/



#include "MEAD/globals.h"
#include "MEAD/PhysCond.h"
#include <iostream>
#include <fstream>

#ifdef USE_EXCEPTIONS
#include "MEAD/MEADexcept.h"
#endif

// Strictly for debugging:
#ifdef POPULATION_COUNT
#include "MEAD/Vertex.h"
#include "MEAD/VertElem.h"
 int Vertex::population = 0;
 int Vertex::max_population = 0;
 int Vertex::num_created = 0;
 int Vertex::num_deleted = 0;
 int VertElem::population = 0;
 int VertElem::max_population = 0;
 int VertElem::num_created = 0;
 int VertElem::num_deleted = 0;
#endif

ofstream cnull ("/dev/null");
ostream *blab1pt = &cnull;
ostream *blab2pt = &cnull;
ostream *blab3pt = &cnull;




const float pi = 3.1415927F;
const float ln10 = 2.3025851F;
const double Boltzmann_constant = 1.380658e-23;  // J/K
// UNUSED const double Gas_constant = 8.314510;  // J/mol/K
const double Avogadro_Number = 6.0221367e+23; // particles/mol
const double elementary_charge = 1.60217733e-19; // Coulombs
const double light_speed = 299792458.0;  // m/s
const double inv_4pi_vacuum_permittivity = light_speed*light_speed*1.0e-7;
                                        // F/m
// The standard energy unit internally to the program is
// charge units squared divided by length, or elem_ch**2/Angst
const double energy_unit = elementary_charge* elementary_charge / 1.0e-10
                          * inv_4pi_vacuum_permittivity;

std::string DOT = ".";

float PhysCond::epsext = 80.0;
float PhysCond::solrad = 1.4F;
float PhysCond::sterln = 2.0F;
float PhysCond::ionicstr = 0.0F;
float PhysCond::T = 300.0;
double PhysCond::kBolt = Boltzmann_constant/energy_unit;
                                   // (proton ch sq)/(Angst * Kelvin)
                                   // or is 5.98339e-6 a better value??
double PhysCond::conconv =  Avogadro_Number * 1.0e-27;
                                     // ((particles/cubic Ang)/(moles/liter))
double PhysCond::econv = energy_unit * Avogadro_Number / 4186.8;
                             // ((kcal/mole)/(proton ch sq/Angst))
                             // Need more accurate value !!!
double PhysCond::hueck = 8 * pi * PhysCond::conconv * PhysCond::ionicstr
                      / PhysCond::kBolt / PhysCond::T;
double PhysCond::kappasq = PhysCond::hueck / PhysCond::epsext;
double PhysCond::ln10kT = ln10 * PhysCond::kBolt * PhysCond::T;
double PhysCond::bohr_radius = 0.529177F;
double PhysCond::proton_charge = 1.0F;

void error (const string& s1, const string& s2, const string& s3)
{
  cerr << s1 << s2 << s3 << "\n";
#ifdef USE_EXCEPTIONS
  throw MEADexcept(s1, s2, s3);
#else
  exit(1);
#endif
}

/* Create and open a file for output, first making sure it doesn't already
exist.
*/
void safeopen(ofstream& f, std::string filename_string)
{
  const char *filename = filename_string.c_str();
  {
    // Try to open it for reading as a portable existence test.
    ifstream file_existence(filename);
    if (file_existence.good()) {
      ::error("ERROR: safeopen: File, ", filename,
	      " already exists\n");
    }
    // This end of block will close/destroy file_existence
  }

  f.open(filename, std::ios::out);
  if (!f) {
    ::error("ERROR: safeopen: File, ", filename,
	    " could not be opened for writing\n");
  }
}

// globals.cc ends here
