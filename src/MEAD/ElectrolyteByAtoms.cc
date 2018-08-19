#include "MEAD/ElectrolyteByAtoms.h"
#include "MEAD/CubeLatSpec.h"
#include "MEAD/SolvAccVol.h"
#include <iostream>

// First some utilities for implementation of ElectrolyteByAtoms --------

struct IShell { // A simplified version of class Shell
  Coord coord;
  float radius;
};

class CLIShell { // A simplified version of CLShell
public:
  inline CLIShell () {}
  CLIShell(const IShell& sph, const CubeLatSpec& cls);
  void zero_inside_points (int *array) const ;
private:
  Coord coord;
  float radius, radius_sq;
  int grid_dim;
  int i1, i2, j1, j2, k1, k2;  // grid box around center
  int is_in;
};

CLIShell::CLIShell(const IShell& sph, const CubeLatSpec& cls)
{
  float spacing = cls.get_spacing();
  grid_dim = cls.get_grid_dim();
  float grlen = (float) (grid_dim - 1);
  float halfgrlen = grlen/2;
  Coord grid_center_in_grid (halfgrlen, halfgrlen, halfgrlen);
  radius = sph.radius/spacing;
  radius_sq = radius*radius;
  coord = (sph.coord - cls.get_center())/spacing + grid_center_in_grid;
  // Now define corners of atom covering box to be scanned in grid-base coords.
  Coord diag(radius, radius, radius);
  Coord bottom_corner = coord - diag;
  Coord top_corner = coord + diag;
  is_in = 1;
  // define The max and min indicies for a box around sphere
  if( bottom_corner.x>=0 && bottom_corner.x<grid_dim)
    i1 = (int) bottom_corner.x;
  else
    if( bottom_corner.x < 0 )
      i1 = 0;
    else
      {is_in = 0; return;}
  if( bottom_corner.y>=0 && bottom_corner.y<grid_dim)
    j1 = (int) bottom_corner.y;
  else
    if( bottom_corner.y < 0 )
      j1 = 0;
    else
      {is_in = 0; return;}
  if( bottom_corner.z>=0 && bottom_corner.z<grid_dim)
    k1 = (int) bottom_corner.z;
  else
    if( bottom_corner.z < 0 )
      k1 = 0;
    else
      {is_in = 0; return;}
  if( top_corner.x>=0 && top_corner.x<grid_dim-1 )
    i2 = (int) top_corner.x+1;
  else
    if( top_corner.x >= grid_dim-1 )
      i2 = grid_dim-1;
    else
      {is_in = 0; return;}
  if( top_corner.y>=0 && top_corner.y<grid_dim-1 )
    j2 = (int) top_corner.y+1;
  else
    if( top_corner.y >= grid_dim-1 )
      j2 = grid_dim-1;
    else
      {is_in = 0; return;}
  if( top_corner.z>=0 && top_corner.z<grid_dim-1 )
    k2 = (int) top_corner.z+1;
  else
    if( top_corner.z >= grid_dim-1 )
      k2 = grid_dim-1;
    else
      {is_in = 0; return;}
}

inline void CLIShell::zero_inside_points (int *array) const
{
  if (!is_in) return;
  int nsq = grid_dim*grid_dim;
  for(int i=i1; i<=i2; i++ ) {
    int insq = i*nsq;
    float dxsq = (i-coord.x)*(i-coord.x);
    for(int j=j1; j<=j2; j++ ) {
      int jline = j*grid_dim;
      float dysq = (j-coord.y)*(j-coord.y);
      for(int k=k1; k<=k2; k++ ) {
	int h = insq + jline + k;
	if(array[h]) {
	  float dsq = (k-coord.z)*(k-coord.z);
	  dsq = dsq + dxsq + dysq;
	  if (dsq < radius_sq)
	    array[h] = 0;
	}
      }
    }
  }
}

struct IonAccVol { // A simplified version of SolvAccVol
  IonAccVol(const AtomSet&, float exclud_radius);
  void calc_cuberep(const CubeLatSpec&, int *array) const ;
  size_t sph_count;
  IShell * ishell;
};

IonAccVol::IonAccVol(const AtomSet& ats, float exclud_radius)
{
  sph_count = ats.size();
  ishell = new IShell [sph_count];
  int atcount = 0;
  for(AtomSet::const_iterator ip=ats.begin(); ip!=ats.end(); ++ip ) {
    const Atom& a = ip->second;
    ishell[atcount].coord = a.coord;
    ishell[atcount].radius = a.rad + exclud_radius;
    atcount++;
  }
  if (atcount != sph_count)
    cerr << "WARNING from IonAccVol constructor: number of atoms conflict.\n";
}

void IonAccVol::calc_cuberep(const CubeLatSpec& cls, int* array) const
{
  int grid_dim = cls.get_grid_dim();
  int ncube = grid_dim*grid_dim*grid_dim;
  for(int i=0; i<ncube; ++i) array[i] = 1;
  CLIShell *cli = new CLIShell [sph_count];
  unsigned s;
  for (s = 0; s<sph_count; ++s)
    cli[s] = CLIShell(ishell[s], cls);
  for (s = 0; s<sph_count; ++s)
    cli[s].zero_inside_points(array);
  delete [] cli;
}

// ---------------- Now the member funcs of EBA itself ---------


ElectrolyteByAtoms::ElectrolyteByAtoms
(const AtomSet& ats, float ionic_strnth, float exclus_radius)
: ElectrolyteEnvironment_lett(), ionic_str(ionic_strnth)
{
  blab2 << "Contructing a ElectrolyteByAtoms" << endl;
  ion_acc_vol = new IonAccVol(ats, exclus_radius);
  sav = 0;
}

ElectrolyteByAtoms::ElectrolyteByAtoms (const AtomSet& ats,
					float ionic_strnth,
					const SolvAccVol& savarg)
: ElectrolyteEnvironment_lett(), ionic_str(ionic_strnth)
{
  blab2 << "Contructing a ElectrolyteByAtoms using an SAV" << endl;
  ion_acc_vol = 0;
  sav = new SolvAccVol(savarg);
}


ElyCubeRep* ElectrolyteByAtoms::get_cuberep(const CubeLatSpec& cls) const
{
  blab2 << "Calculating an ElyCubeRep for this ElectrolyteByAtoms"
    << endl;
  ElyCubeRep *p = new ElyCubeRep(cls);
  if (sav == 0 && ion_acc_vol != 0) {
    ion_acc_vol->calc_cuberep(cls, p->isgarr);
  }
  else if (sav != 0 && ion_acc_vol == 0) {
    int ngrid = cls.get_grid_dim();
    int ncube = ngrid*ngrid*ngrid;
    AccTag* accarr = new AccTag[ncube];
    sav->calc_cuberep(cls, accarr);
    int errcount = 0;
    for (int h=0; h<ncube; ++h) {
      if (accarr[h] == interior)
	p->isgarr[h] = 0;
      else if (accarr[h] == exterior)
	p->isgarr[h] = 1;
      else
	++errcount;
    }
    if (errcount)
      cerr << "ERROR: ElectrolyteByAtoms::get_cuberep:\n"
	<< "The AccTag array returned by the SolvAccVol had " << errcount
	  << "\nelements that were neither interior nor exterior.\n"
	    << "Continuing anyway." << endl;
  }
  else {
    ::error("INTERNAL ERROR IN ElectrolyteByAtoms::get_cuberep:\n",
	    "sav and ion_acc_vol pointers in a senseless state.\n");
  }
  return p;
}

// ElectrolyteByAtoms.cc ends here
