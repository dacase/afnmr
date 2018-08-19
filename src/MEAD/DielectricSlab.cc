#include "MEAD/DielectricSlab.h"
#include "MEAD/globals.h"

DielectricSlab::DielectricSlab(float eslab, float eext, float zl, float zu)
 : DielectricEnvironment_lett(),
   epsslab(eslab), epsext(eext), zlower(zl), zupper(zu)
{
  // Sanity checks
  if (zlower >= zupper) {
    cerr << "ERROR: DielectricSlab ctor called with zlower >= zupper\n"
      << "(" << zlower << " >= " << zupper << ")" << endl;
    ::error();
  }
  if (epsslab < 0.99) {
    cerr << "ERROR: DielectricSlab ctor called with epsslab < 0.99\n"
      << "(" << epsslab << " < 0.99)" << endl;
    ::error();
  }
  if (epsext < 0.99) {
    cerr << "ERROR: DielectricSlab ctor called with epsext < 0.99\n"
      << "(" << epsext << " < 0.99)" << endl;
    ::error();
  }
}


DielCubeRep DielectricSlab::get_cuberep (const CubeLatSpec& cls)
{
  blab3 << "DielectricSlab::get_cuberep called.\n"
    << "\nSetting all cubic dielectric grid points between z = " << zupper
      << " and z = " << zlower << " to a value of " << epsslab
	<< "\nand other points to " << epsext << endl;

  int grid_dim = cls.get_grid_dim();
  int nsq = grid_dim*grid_dim;
  float spacing = cls.get_spacing();
  Coord grid_center_in_space = cls.get_center();

  float grlen = (float) (grid_dim - 1);
  float halfgrlen = grlen/2;

  Coord grid_center_in_grid (halfgrlen, halfgrlen, halfgrlen);

  float zupper_in_grid = (zupper - grid_center_in_space.z)
    / spacing + grid_center_in_grid.z;
  float zlower_in_grid = (zlower - grid_center_in_space.z)
    / spacing + grid_center_in_grid.z;


  DielCubeRep dcr(cls);
  for (int i = 0; i<grid_dim; ++i) {
    int iterm = i*nsq;
    for (int j = 0; j<grid_dim; ++j) {
      int jterm = j*grid_dim;
      int ijterms = iterm + jterm;
      for (int k = 0; k<grid_dim; ++k) {
	int h = ijterms + k;
	float z = (float) k;
	dcr[h] = z > zlower_in_grid && z < zupper_in_grid ? epsslab : epsext;
      }
    }
  }

  dcr.declare_defined();
  return dcr;
}

// DielectricSlab.cc ends here
