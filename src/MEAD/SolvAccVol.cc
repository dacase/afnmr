/* class for calculations relating to solvent-accessible volumes of molecules
    Copyright (c) 1993--1995 by Donald Bashford and Tony You.

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

$Id: SolvAccVol.cc,v 2.14 2004/12/06 17:58:39 bashford Exp $
*/

#include "MEAD/SolvAccVol.h"
#include "MEAD/PhysCond.h"
#include "MEAD/AtomSet.h"
#include "MEAD/Shell.h"
#include "MEAD/Sausage.h"
#include <sstream>
using std::ios;
#include <string>
using std::string;


SolvAccVol::SolvAccVol(const AtomSet& a)
{
  rep = new SolvAccVolRep(a, PhysCond::get_solrad());
}

SolvAccVol::SolvAccVol(const AtomSet& a, float probe_radius)
{
  rep = new SolvAccVolRep(a, probe_radius);
}

SolvAccVol::SolvAccVol(const SolvAccVol& s)
{
  s.rep->referenceCount++;
  rep = s.rep;
}

SolvAccVol&
SolvAccVol::operator= (const SolvAccVol& s)
{
  if (this != &s) {
    s.rep->referenceCount++;
    if (rep->referenceCount == 1)
      delete rep;
    else
      --rep->referenceCount;
    rep = s.rep;
  }
  return *this;
}

SolvAccVol::~SolvAccVol()
{
  if (rep->referenceCount == 1)
    delete rep;
  else
    --rep->referenceCount;
}

void
SolvAccVol::anal_calc()
{
  rep->anal_calc();
}

void
SolvAccVol::calc_cuberep(const CubeLatSpec& cls, AccTag* atg)
{
  rep->calc_cuberep(cls, atg);
}

void
SolvAccVol::tag_points(int npts, const Coord* pt, AccTag* acc_array)
{
  rep->tag_points(npts, pt, acc_array);
}

int
SolvAccVol::accessible(const Coord& r)
{
  return rep->accessible(r);
}

int
SolvAccVol::check_is_calculated() const
{
  return rep->check_is_calculated();
}

void SolvAccVol::write_top_in_binary() const
{
 rep->write_top_in_binary();
}

void SolvAccVol::write_top_in_ascii() const
{
 rep->write_top_in_ascii();
}

int SolvAccVol::read_top_in_binary()
{
 return rep->read_top_in_binary();
}

int SolvAccVol::read_top_in_ascii()
{
 return rep->read_top_in_ascii();
}

void SolvAccVol::write_top_in_binary(ostream& ost) const
{
 rep->write_top_in_binary(ost);
}

void SolvAccVol::write_top_in_ascii(ostream& ost) const
{
 rep->write_top_in_ascii(ost);
}

int SolvAccVol::read_top_in_binary(istream& ist)
{
 return rep->read_top_in_binary(ist);
}

int SolvAccVol::read_top_in_ascii(istream& ist)
{
 return rep->read_top_in_ascii(ist);
}

void SolvAccVol::write_top_in_binary(const string& filename) const
{
 rep->write_top_in_binary(filename);
}

void SolvAccVol::write_top_in_ascii(const string& filename) const
{
 rep->write_top_in_ascii(filename);
}

int SolvAccVol::read_top_in_binary(const string& filename)
{
 return rep->read_top_in_binary(filename);
}

int SolvAccVol::read_top_in_ascii(const string& filename)
{
 return rep->read_top_in_ascii(filename);
}

//-------- SolvAccVolRep, where the real work happens -----------

int SolvAccVolRep::instanceCount = 0;

SolvAccVolRep::SolvAccVolRep(const AtomSet& a, float probe_radius)
: Rprob(probe_radius)
{
  blab2 << "SolvAccVolRep constructor called with AtomSet" << endl;
  ++instanceCount;
  referenceCount = 1;
  sph_count = 0;
  sausage_count = 0;
  is_calculated = 0;
  sglist = 0;
  atsph = 0;
  // Create a name from the AtomSet and radius.
  // This will be used for save-file name so must be unique for a given
  // molecule coords, radii and probe radii.
  // Assume name of the AtomSet is already unique for coords and atom radii.
  {
    std::ostringstream ost;
    ost << a.get_name() << "_Rprob_" << Rprob << "_anal_SolvAccVol"
	<< std::flush;
    name = ost.str();
  }

  // Transfer atom data to the shell list.
  // Later in anal_calc it will be necessary to eliminate completely
  // buried atoms in a rewrite of shell list.
  if( a.size() == 0 ) { // No atoms is not a problem!
    blab2 << "SolvAccVol constructor No atoms so empty data in this SolvAccVol"
      << "\nAll space will be free"  << endl;
    sph_count = 0;
    is_calculated = 1;
    return;
  }
  atsph = new Shell [a.size()];
  int atcount = 0;
  for(AtomSet::const_iterator ip = a.begin(); ip!=a.end(); ++ip ) {
    Shell sp(ip->second, Rprob);
    if (sp.get_inner_rad() <= 0.0) continue;  // Zero rad atoms not seen.
    atsph[atcount] = sp;
    atcount++;
  }
  sph_count = atcount;
  blab3 << "Finishing copy atomic information from AtomSet\t# of atom ="
    << atcount << endl;;
}

SolvAccVolRep::~SolvAccVolRep()
{
  if (--referenceCount)
    cerr << "WARNING: SolvAccVolRep destructor called with referenceCount = "
      << referenceCount << endl;
  if (sph_count) delete [] atsph;
  if (sausage_count) delete [] sglist;
  --instanceCount;
}

int
SolvAccVolRep::check_is_calculated() const
{return is_calculated;}

int
SolvAccVolRep::instances()
{return instanceCount;}

void SolvAccVolRep::write_top_in_binary() const
{
  string filename = name + ".dat";
  write_top_in_binary(filename);
}

void SolvAccVolRep::write_top_in_binary(const string& filename_string) const
{
  const char *filename = filename_string.c_str();
  ofstream mol_top_in_binary(filename, ios::out | ios::binary);
  if( !mol_top_in_binary ) {
    cerr << "ERROR: file_name: " << filename << "\n";
    cerr <<"SolvAccVolRep::write_top_in_binary: couldn't open for witing\n";
  }
  else
    write_top_in_binary(mol_top_in_binary);
}

void SolvAccVolRep::write_top_in_binary(ostream& mol_top_in_binary) const
{
  if (!is_calculated) {
    cerr << "WARNING: SolvAccVolRep::write_top_in_binary called\n"
      << "when analytical representation not yet calculated" << endl;
    return;
  }
  mol_top_in_binary.write((char *) &Rprob, sizeof(float));
  mol_top_in_binary.write((char *) &sph_count, sizeof(float));
  int n;
  for(n=0; n<sph_count; n++) {
    atsph[n].write_top_in_binary(mol_top_in_binary);
  }
  mol_top_in_binary.write((char *) &sausage_count, sizeof(float));
  for( n=0; n<sausage_count; n++) {
    sglist[n].write_top_in_binary(mol_top_in_binary);
  }
  return;
}

void SolvAccVolRep::write_top_in_ascii() const
{
  string filename = name + ".txt";
  write_top_in_ascii(filename);
}

void SolvAccVolRep::write_top_in_ascii(const string& filename_string) const
{
  const char *filename = filename_string.c_str();
  ofstream mol_top_in_ascii(filename);
  if( !mol_top_in_ascii.good() ) {
    cout << "file_name: " << filename << "\n";
    cerr <<"SolvAccVolRep::write_trop_in_binary: couldn't open for witing\n";
  }
  else
    write_top_in_ascii(mol_top_in_ascii);
}
void SolvAccVolRep::write_top_in_ascii(ostream& mol_top_in_ascii) const
{
  if (!is_calculated) {
    cerr << "WARNING: SolvAccVolRep::write_top_in_binary called\n"
      << "when analytical representation not yet calculated" << endl;
    return;
  }
  mol_top_in_ascii << "Solvent Accessible Volume Elements of " << name << endl;
  mol_top_in_ascii << endl;
  mol_top_in_ascii << "Probe radius is " << Rprob << endl;
  mol_top_in_ascii << endl;
  mol_top_in_ascii << "Number of shells is " << sph_count << endl;
  mol_top_in_ascii << endl;
  int n;
  for(n=0; n<sph_count; n++)
     {
      mol_top_in_ascii << "Shell " << n+1 << "   ";
      atsph[n].write_top_in_ascii(mol_top_in_ascii);
      mol_top_in_ascii << endl;
     }
  mol_top_in_ascii << "Number of sausages is " << sausage_count << endl;
  for( n=0; n<sausage_count; n++)
     {
      mol_top_in_ascii << "Sausage " << n+1 << "   ";
      sglist[n].write_top_in_ascii(mol_top_in_ascii);
      mol_top_in_ascii << endl;
     }
}

int SolvAccVolRep::read_top_in_binary()
{
  string filename = name + ".dat";
  return read_top_in_binary(filename);
}
int SolvAccVolRep::read_top_in_binary(const string& filename_string)
{
  const char *filename = filename_string.c_str();
  blab2 << "SolvAccVol::read_top_in_binary(" << filename << ") called" << endl;
  ifstream ftop_file(filename, ios::in | ios::binary);
  if ( !ftop_file.good() )
    return 0;
  return read_top_in_binary(ftop_file);
}

int SolvAccVolRep::read_top_in_binary(istream& ftop_file)
{
 float sol_prob;
 ftop_file.read( (char*) &sol_prob, sizeof(float) );
 if( !ftop_file )
   {
    cerr << "ERROR : SolvAccVolRep::read_top_in_binary: " << endl;
    cerr << "Unable to read solvent probe radius" << endl;
    return 0;
   }
 if (Rprob != sol_prob)
   cerr << "ERROR: SolvAccVolRep::read_top_in_binary:\n"
     << "Solvent probe radius from file does not match current data.\n"
       << "Will ignore value from file and continue." << endl;
 ftop_file.read( (char*) &sph_count, sizeof(float) );
 blab2 << "sph_count = " << sph_count << endl;
 if( !ftop_file )
   {
    cerr << "ERROR : SolvAccVolRep::read_top_in_binary: " << endl;
    cerr << "Unable to read shell_count" << endl;
    return 0;
   }
 if( sph_count <= 0 )
   {
    cerr << "WARNING : SolvAccVolRep::read_top_in_binary: " << endl;
    cerr << "sph_count <= zero is detected." << endl;
   }
 Shell *shell_list = new Shell[sph_count];
 int n;
 for(n=0; n<sph_count; n++)
    {
     if( shell_list[n].read_top_in_binary(ftop_file) == 0 )
         return 0;
    }
 blab2 << "Number of shells read in = " << sph_count << endl;
 ftop_file.read( (char*) &sausage_count, sizeof(int) );
 if( !ftop_file )
   {
    cerr << "ERROR : SolvAccVolRep::read_top_in_binary: " << endl;
    cerr << "Unable to read sausage_count " << endl;
    return 0;
   }
 if( sausage_count <= 0 )
   {
    cerr << "WARNING : SolvAccVolRep::read_top_in_binary: " << endl;
    cerr << "sausage_count <= zero is detected." << endl;
   }
 Sausage *sausage_list = new Sausage[sausage_count];
 for( n=0; n<sausage_count; n++)
    {
     if( sausage_list[n].read_top_in_binary(ftop_file) == 0 )
         return 0;
    }
 blab2 << "Number of sausages read in = " << sausage_count << endl;
 delete [] atsph;
 delete [] sglist;
 atsph = shell_list;
 sglist = sausage_list;
 return is_calculated = 1;
}

int SolvAccVolRep::read_top_in_ascii()
{
  string filename = name + ".txt";
  return read_top_in_ascii(filename);
}
int SolvAccVolRep::read_top_in_ascii(const string& filename_string)
{
  const char *filename = filename_string.c_str();
  blab2 << "SolvAccVol::read_top_in_ascii(" << filename << ") called" << endl;
  ifstream ftop_file(filename, ios::in);
  if ( !ftop_file.good() )
    return 0;
  return read_top_in_ascii(ftop_file);
}

int SolvAccVolRep::read_top_in_ascii(istream& ftop_file)
{
 string word;
 int n;
 for(n=0; n<9; n++ )
      ftop_file >> word;
 float sol_prob;
 ftop_file >> sol_prob;
 if (Rprob != sol_prob)
   cerr << "ERROR: SolvAccVolRep::read_top_in_binary:\n"
     << "Solvent probe radius from file does not match current data.\n"
       << "Will ignore value from file and continue." << endl;
 for(n=0; n<4; n++ )
      ftop_file >> word;
 ftop_file >> sph_count;
 Shell *shell_list = new Shell[sph_count];
 for( n=0; n<sph_count; n++ )
    {
     if( shell_list[n].read_top_in_ascii(ftop_file) == 0 )
       {
        cerr << "ERROR : SolvAccVolRep::read_top_in_ascii : " << endl;
        return 0;
       }
    }
 blab2 << "Number of shells read in = " << sph_count << endl;
 for( n=0; n<4; n++ )
      ftop_file >> word;
 ftop_file >> sausage_count;
 Sausage *sausage_list = new Sausage[sausage_count];
 for( n=0; n<sausage_count; n++ )
    {
     if( sausage_list[n].read_top_in_ascii(ftop_file) == 0 )
       {
        cerr << "ERROR : SolvAccVolRep::read_top_in_ascii : " << endl;
        return 0;
       }
    }
 blab2 << "Number of sausages read in = " << sausage_count << endl;
 delete [] atsph;
 delete [] sglist;
 atsph = shell_list;
 sglist = sausage_list;
 return is_calculated = 1;
}

// SolvAccVol.cc ends here
