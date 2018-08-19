#include "MEAD/UniformElectrolyte.h"
#include "MEAD/CubeLatSpec.h"
#include <iostream>


UniformElectrolyte::UniformElectrolyte(float ionic)
 : ElectrolyteEnvironment_lett(), ionic_str(ionic)
{
}

ElyCubeRep* UniformElectrolyte::get_cuberep(const CubeLatSpec& cls) const
{
  ElyCubeRep* p = new ElyCubeRep(cls);
  int n = cls.get_grid_dim();
  n=n*n*n;
  for(int i=0; i<n; ++i) p->isgarr[i] = 1;
  return p;
}

// UniformElectrolyte.cc ends here
