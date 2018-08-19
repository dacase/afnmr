#include "MEAD/UniformDielectric.h"
#include "MEAD/DielCubeRep.h"


UniformDielectric::UniformDielectric(float e)
: DielectricEnvironment_lett(), eps(e)
{
}

DielCubeRep
UniformDielectric::get_cuberep (const CubeLatSpec& cls)
{
  blab3 << "UniformDielectric::get_cuberep: called"
	<< "\nSet the whole cubic dielectric grid to " << eps << endl;
  DielCubeRep dcr(cls);
  int grid_dim = cls.get_grid_dim();
  int ncube = grid_dim*grid_dim*grid_dim;
  for (int h = 0; h<ncube; ++h)
    dcr[h] = eps;
  dcr.declare_defined();
  return dcr;
}

// UniformDielectric.cc ends here
