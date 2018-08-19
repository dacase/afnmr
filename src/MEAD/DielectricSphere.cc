#include "MEAD/DielectricSphere.h"
#include "MEAD/DielCubeRep.h"


DielectricSphere::DielectricSphere (float ein, float eext,
				    float rad, Coord cntr)
 : DielectricEnvironment_lett(),
   epsin (ein), epsext (eext), radius(rad), center(cntr)
{
}

DielCubeRep
DielectricSphere::get_cuberep (const CubeLatSpec& cls)
{
  blab3 << "DielectricSphere::get_cuberep called.\n"
    << "Setting all cubic dielectric grid points within " << radius
      << " of center to a value of " << epsin
	<< "\nand other points to " << epsext << endl;

// FIXME?  If I was sure inlining of the cls functions produced just
// as fast code I could use cls.get.. () directly instead of these
// local variables.

  int grid_dim = cls.get_grid_dim();
  int nsq = grid_dim*grid_dim;
  float spacing = cls.get_spacing();
  Coord grid_center_in_space = cls.get_center();

  float grlen = (float) (grid_dim - 1);
  float halfgrlen = grlen/2;

  Coord grid_center_in_grid (halfgrlen, halfgrlen, halfgrlen);
  Coord sphere_center_in_grid = (center - grid_center_in_space)
    / spacing + grid_center_in_grid;
  float radius_squared = radius / spacing;
  radius_squared *= radius_squared;

  DielCubeRep dcr(cls);
  for (int i = 0; i<grid_dim; ++i) {
    int iterm = i*nsq;
    float dxsq = ((float) i) - sphere_center_in_grid.x;
    dxsq *= dxsq;
    for (int j = 0; j<grid_dim; ++j) {
      int jterm = j*grid_dim;
      int ijterms = iterm + jterm;
      float dysq = ((float) j) - sphere_center_in_grid.y;
      dysq *= dysq;
      float dxsq_plus_dysq = dxsq + dysq;
      for (int k = 0; k<grid_dim; ++k) {
	int h = ijterms + k;
	float dzsq = ((float) k) - sphere_center_in_grid.z;
	dzsq *= dzsq;
	float rsq = dxsq_plus_dysq + dzsq;
	dcr[h] = rsq > radius_squared ? epsext : epsin;
      }
    }
  }

  dcr.declare_defined();
  return dcr;
}

// DielectricSphere.cc ends here
