/* Hold data on the physical conditions, e.g. temperature.
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

$Id: PhysCond.cc,v 1.10 2004/12/06 17:58:31 bashford Exp $
*/


/*
PhysCond is a purely static class containing parameters describing
the physical conditions (dielectric constants, ionic strength,
solvent probe radius, etc.), some units conversion factors, and
the necessary access functions to set and get these values.
The main program should have a section that reads command line flags
to set these values.
*/

#include "MEAD/PhysCond.h"
#include <iostream>
#include <math.h>


void PhysCond::print()
{
  cout << "Exterior dielectric constant, espext = " << epsext << "\n";
  cout << "Solvent probe radius, solrad = " << solrad << "\n";
  cout << "Ion exclusion layer thickness, sterln = " << sterln << "\n";
  cout << "Temperature, T = " << T << "\n";
  cout << "Ionic strength, ionicstr = " << ionicstr << "\n";
  cout << "Hueckel factor (epsext*kappasq) = " << hueck << "\n";
  if (kappasq)  cout << "Debye length 1/kappa = " << 1/sqrt(kappasq) << "\n";
  cout << "ln(10) * kBolt * T = " << ln10kT << "\n";
  cout << "kBolt = " << kBolt << endl;
  cout << "conconv = " << conconv << endl;
  cout << "econv = " << econv << endl;
  cout << "Bohr radius = " << bohr_radius << endl;
  cout << "Proton Charge = " << proton_charge << endl;
}

void PhysCond::print(ostream& os)
{
  os << "Exterior dielectric constant, espext = " << epsext << "\n";
  os << "Solvent probe radius, solrad = " << solrad << "\n";
  os << "Ion exclusion layer thickness, sterln = " << sterln << "\n";
  os << "Temperature, T = " << T << "\n";
  os << "Ionic strength, ionicstr = " << ionicstr << "\n";
  os << "Hueckel factor (epsext*kappasq) = " << hueck << "\n";
  if (kappasq)  os << "Debye length 1/kappa = " << 1/sqrt(kappasq) << "\n";
  os << "ln(10) * kBolt * T = " << ln10kT << "\n";
  os << "kBolt = " << kBolt << endl;
  os << "conconv = " << conconv << endl;
  os << "econv = " << econv << endl;
}

ostream & operator <<(ostream& os, PhysCond& p)
{
  p.print(os);
  return os;
}

// PhysCond.cc ends here
