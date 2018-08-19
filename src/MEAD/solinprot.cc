/* Calculate "solvation energy" of a group in protein/water system.

    Copyright (c) 1993--1995 by Donald Bashford

    This source code file is part of a Scripps in-house version of the
    MEAD (Macroscopic Electrostatics with Atomic Detail) package of
    objects and programs.  It should not be redistributed beyond
    Scripps without the permission of Donald Bashford.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

    Donald Bashford can be contacted by electronic mail by the address,
    bashford@scripps.edu, or by paper mail at Department of Molecular
    Biology, The Scripps Research Institute, 10666 North Torrey Pines
    Road, La Jolla, California 92037.

$Id: solinprot.cc,v 2.14 2013/10/30 21:30:41 bashford Exp $
*/
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <sstream>
using std::istringstream;
#include <fstream>
using std::ifstream;
#include <string>
using std::string;

#include "MEAD/PhysCond.h"
#include "MEAD/AtomSet.h"
#include "MEAD/ChargeDist.h"
#include "MEAD/AtomChargeSet.h"
#include "MEAD/DielectricEnvironment.h"
#include "MEAD/DielByAtoms.h"
#include "MEAD/DielMembAtoms.h"

#include "MEAD/ElectrolyteEnvironment.h"
#include "MEAD/ElectrolyteByAtoms.h"
#include "MEAD/FinDiffMethod.h"
#include "MEAD/ElstatPot.h"
#include "MEAD/FDGridLevel.h"

// For DEBUGGING
// extern void big_rigid_stats();
// extern void big_rigid_delete_all();
// extern "C" {
//   void malloc_stats();
// };

int main(int argc, char* argv[])
{
  string molname;
  string protname;
  float epsvac = 1.0;
  float epsin1 = 0.0;  // These zeros are error states.  A specification of
  float epsin2 = 0.0;  // epsin1 and epsin2 on the command line is required.
  float doReactionField = 0;
  float doProteinField = 0;

  int ismemb = 0;
  float ztop=0, zbot=0;
  float holerad = 0, holex = 0, holey = 0;

  // Parse the command line
  for (int iarg=1; iarg<argc; ++iarg) {
    if (iarg+2==argc) {  // These is the last one, must be molname protname
      molname = argv[iarg];
      ++iarg;
      protname = argv[iarg];
    }
    else if ((string) argv[iarg] == (string) "-epsin1") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -epsin1");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -epsin1 option");
      epsin1 = f;
    }
    else if ((string) argv[iarg] == (string) "-epsin2") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -epsin2");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -epsin2 option");
      epsin2 = f;
    }
    else if ((string) argv[iarg] == (string) "-epsext") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -epsext");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -epsext option");
      PhysCond::set_epsext(f);
    }
    else if ((string) argv[iarg] == (string) "-epsvac") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -epsvac");
      istringstream ist(argv[++iarg]);
      ist >> epsvac;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -epsvac option");
    }
    else if ((string) argv[iarg] == (string) "-solrad") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -solrad");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -solrad option");
      PhysCond::set_solrad(f);
    }
    else if ((string) argv[iarg] == (string) "-sterln") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -sterln");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -sterln option");
      PhysCond::set_sterln(f);
    }
    else if ((string) argv[iarg] == (string) "-ionicstr") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -ionicstr");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -ionicstr option");
      PhysCond::set_ionicstr(f);
    }
    else if ((string) argv[iarg] == (string) "-T") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -T");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -T option");
      PhysCond::set_T(f);
    }
    else if ((string) argv[iarg] == (string) "-kBolt") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -kBolt");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -kBolt option");
      PhysCond::set_kBolt(f);
    }
    else if ((string) argv[iarg] == (string) "-conconv") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -conconv");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -conconv option");
      PhysCond::set_conconv(f);
    }
    else if ((string) argv[iarg] == (string) "-econv") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -econv");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -econv option");
      PhysCond::set_econv(f);
    }
    else if ((string) argv[iarg] == (string) "-bohr_radius") {
      if (iarg+1>=argc)
	error ("main: ERROR, not enough args after -bohr_radius");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -bohr_radius option");
      PhysCond::set_bohr_radius(f);
    }
    else if ((string) argv[iarg] == (string) "-proton_charge") {
      if (iarg+1>=argc)
	error ("main: ERROR, not enough args after -proton_charge");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -proton_charge opt.");
      PhysCond::set_proton_charge(f);
    }
    // process the form "-membz float float"
    else if ((string) argv[iarg] == (string) "-membz") {
      ismemb = 1;
      if (iarg+2>=argc) error ("main: ERROR, not enough args after -membz");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read first value after -membz opt.");
      zbot = f;
      istringstream ist2(argv[++iarg]);
      ist2 >> f;
      if (!(ist2.good() || ist2.eof()))
	error ("main: FAILED trying to read second value after -membz opt.");
      ztop = f;
    }
    // process the form "-membhole float [float float]"
    else if ((string) argv[iarg] == (string) "-membhole") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -membhole");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read first value after -memhole opt.");
      holerad = f;
      if (iarg+2 < argc) { // There _could_ be two more floats...
	istringstream istx(argv[iarg+1]);
	istx >> f;
	if (istx.good() || istx.eof()) { // Yes, expect this float and another.
	  holex = f;
	  iarg+=2;
	  istringstream isty(argv[iarg]);
	  isty >> f;
	  if (!(isty.good() || isty.eof()))
	    error("main: ERROR, only two floats seen after -membhole.\n",
		  "One or three expected.\n");
	  else
	    holey = f;
	}
      }
    }
    else if ((string) argv[iarg] == (string) "-ReactionField")
      doReactionField = 1;
    else if ((string) argv[iarg] == (string) "-ProteinField")
      doProteinField = 1;
    else if ((string) argv[iarg] == (string) "-converge_oldway")
      FDGridLevel::set_use_fixed_maxrmsdiff(true);
    // Process -blab flags
    else if ((string) argv[iarg] == (string) "-blab1") {
      blab1pt = &cout;
    }
    else if ((string) argv[iarg] == (string) "-blab2") {
      blab1pt = &cout;
      blab2pt = &cout;
    }
    else if ((string) argv[iarg] == (string) "-blab3") {
      blab1pt = &cout;
      blab2pt = &cout;
      blab3pt = &cout;
    }
    else {
      cerr << "main: ERROR command option, " << argv[iarg]
	<< ", not recognized" << endl;
      error();
    }
  }

  if (molname == (string) "")
    error("main: ERROR, molname not given on command line");
  if (epsin1 == 0.0 && epsin2 == 0.0)
    ::error("main: ERROR, -epsin1 and -epsin2 values not specified.",
	    "   These flags are now mandatory.");
  else if (epsin1 == 0.0)
    ::error("main: ERROR, -epsin1 value not specified.",
	    "   This flag is now mandatory.");
  else if (epsin2 == 0.0)
    ::error("main: ERROR, -epsin2 value not specified.",
	    "   This flag is now mandatory.");

  // Done with comand line processing
  // Print summary of input info from the command line:
  cout << "Starting " << argv[0] << " for molecule named " << molname << "\n";
  cout << "using the following physical conditions:\n";
  PhysCond::print();
  cout << "  Solute interior dielectric, epsin1 = " << epsin1
    << "\n  Protein interior dielectric constant, epsin2 = " << epsin2
      << "\nVacuum dielectric constant, epsvac = " << epsvac
	<< "\n" << endl;
  if (ismemb) {
    cout << "A membrane is in the region from z = " << zbot << " to " << ztop;
    cout << ",\nwith a hole of radius " << holerad << " at (x,y) = ("
      << holex << "," << holey << ")\n";
    if (PhysCond::get_ionicstr() != 0)
      error ("Ability to have membrane and nonzero ionic strength not implemented\n");
  }
  if (blab3pt == &cout) cout << "Blab level set to " << 3 << endl;
  else if (blab2pt == &cout) cout << "Blab level set to " << 2 << endl;
  else if (blab1pt == &cout) cout << "Blab level set to " << 1 << endl;
  else cout << "No blab level set (so no blabbing)" << endl;

//  cout << "Malloc stats before reading atoms" << endl;
//  malloc_stats();
//  big_rigid_stats();
//  {

  AtomSet a1(molname);
  a1.read();
  ChargeDist rhogen(new AtomChargeSet(a1));

  AtomSet a2(protname);
  a2.read();
  ChargeDist rhofeel(new AtomChargeSet(a2));

  DielectricEnvironment_lett* delp = 0;
  if (ismemb) {
    Coord c(holex, holey, 0.0);
    delp = new ThreeValueDielMembAtomsAtoms(a1, PhysCond::get_solrad(),
					    a2, PhysCond::get_solrad(),
					    ztop, zbot, c, holerad,
					    epsin1, epsin2,
					    PhysCond::get_epsext());
  }
  else {
    delp = new ThreeValueDielectricByAtoms(a1, epsin1, a2, epsin2);
  }
  DielectricEnvironment eps(delp);

  // Make the union of a1 and a2 for defining electrolyte boundary
  // (FIXME I think STL provides a function to do this.)
  AtomSet u(a1);
  for (AtomSet::const_iterator i = a2.begin(); i!=a2.end(); ++i) {
    const Atom& a = i->second;
    AtomID k = i->first;
    if (a1.contains(k)) {
      cerr << argv[0] << ": WARNING:  The atom, \"" << k
	<< "\",\n   occurs in both the solute and protein atom sets.\n"
	  << "   For defining the Elelctrolyte boundary, the protein version\n"
	    << "   will be used."  << endl;
    }
    else {
      u.insert(a);
    }
  }
  ElectrolyteEnvironment ely(new ElectrolyteByAtoms(AtomSet(u)));

  FinDiffMethod fdm;
  string ogm_filename = molname + ".ogm";
  fdm.read(ogm_filename);
  Coord interesting (0,0,0);
  fdm.resolve(a1.geom_cent(), interesting);
  cout << "Using finite difference method with lattice levels:" << endl;
  cout << fdm;

  float prod_sol, protein_interaction;

  ElstatPot phi(fdm, eps, rhogen, ely);
  phi.solve();
  prod_sol = phi * rhogen;
  cout << "prod_sol = " << prod_sol << endl;
  protein_interaction = phi * rhofeel;
  protein_interaction *=  PhysCond::get_econv();
  cout << "Interaction of solute with protein charges = "
    << protein_interaction << "\n(probably in kcal/mole)" << endl;


  float safe_eps_sol = PhysCond::get_epsext(); // What for?
  PhysCond::set_epsext(epsvac);
  ElectrolyteEnvironment vac_ely;  // No electrolyte is the default
  DielectricEnvironment vac_eps(new TwoValueDielectricByAtoms(a1, epsin1));
  ElstatPot vac_phi(fdm, vac_eps, rhogen, vac_ely);
  vac_phi.solve();
  float prod_vac = vac_phi * rhogen;
  cout << "prod_vac = " << prod_vac << endl;
  float reac_energy = (prod_sol - prod_vac) / 2 * PhysCond::get_econv();
  cout << "\n\nReaction field component of solvation = " << reac_energy
    << "\n(probably in kcal/mole)" << endl;
  float solvation_energy = reac_energy + protein_interaction;
  cout << "\n\nSOLVATION ENERGY IN PROTEIN = " << solvation_energy
    << "\n(probably in kcal/mole)" << endl;

  if (doReactionField) {
    string fieldpoint_filename = molname + ".fpt" ;
    ifstream fpt (fieldpoint_filename.c_str());
    if (!fpt.good())
      ::error(argv[0], "main: failed to open for reading, ",
	      fieldpoint_filename.c_str());
    string reacfield_filename = molname + ".rf";
    ofstream rf(reacfield_filename.c_str());
    if (!rf.good())
      ::error(argv[0], "main: failed to open for writing, ",
	      reacfield_filename.c_str());
    cout << "Reaction field at points in " << fieldpoint_filename
      << " will be written to " << reacfield_filename << endl;
    Coord c;
    while (fpt >> c) {
      rf << phi.value(c) - vac_phi.value(c) << endl;
      if (!rf.good())
	::error(argv[0], "main: output failure while writing ",
		reacfield_filename.c_str());
    }
    // Expect fpt to be at EOF, else some kind of error or warning
    if (! fpt.eof()) {
      if (fpt.bad())
	::error(argv[0], "main: serious input failure while reading ",
		fieldpoint_filename.c_str());
      else if (fpt.fail()) {
	cerr << "WARNING: main: Reading of file, " << fieldpoint_filename
	     << ", ended with some unexpected condition, \n"
	     << "such as encountering an unexpected character while trying to read a float\n"
	     << "or a missing close parenthesis."
	     << endl;
      }
    }
  }

  if (doProteinField) {
    blab1 << "Protein Field requested, so calculating field due to"
      << "protein charges" << endl;
    PhysCond::set_epsext(safe_eps_sol);
    ElstatPot prot_phi(fdm, eps, rhofeel, ely);
    prot_phi.solve();
    string fieldpoint_filename = molname + ".fpt" ;
    ifstream fpt (fieldpoint_filename.c_str());
    if (!fpt.good())
      ::error(argv[0], "main: failed to open for reading, ",
	      fieldpoint_filename.c_str());
    string protfield_filename = molname + ".pf";
    ofstream pf(protfield_filename.c_str());
    if (!pf.good())
      ::error(argv[0], "main: failed to open for writing, ",
	      protfield_filename.c_str());
    cout << "Protein field at points in " << fieldpoint_filename
      << " will be written to " << protfield_filename << endl;
    Coord c;
    while (fpt >> c) {
      pf << prot_phi.value(c) << "\n";
      if (!pf.good())
	::error(argv[0], "main: output failure while writing ",
		protfield_filename.c_str());
    }
    // Expect fpt to be at EOF, else some kind of error or warning
    if (! fpt.eof()) {
      if (fpt.bad())
	::error(argv[0], "main: serious input failure while reading ",
		fieldpoint_filename.c_str());
      else if (fpt.fail()) {
	cerr << "WARNING: main: Reading of file, " << fieldpoint_filename
	     << ", ended with some unexpected condition, \n"
	     << "such as encountering an unexpected character while trying to read a float\n"
	     << "or a missing close parenthesis."
	     << endl;
      }
    }
  }

//  }
//  big_rigid_delete_all();
//  cout << "Malloc stats at end of Main" << endl;
//  malloc_stats();
//  big_rigid_stats();
  return 0;

}

// solinprot.cc ends here
