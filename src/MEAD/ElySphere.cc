#include "MEAD/ElySphere.h"
#include "MEAD/CubeLatSpec.h"


ElySphere::ElySphere(float ionic, Coord c, float r)
 : ElectrolyteEnvironment_lett(), ionic_str(ionic)
{
  center = c;
  radius = r;
}
ElyCubeRep* ElySphere::get_cuberep(const CubeLatSpec& cls) const
{
  ElyCubeRep* p = new ElyCubeRep(cls);
  int n = cls.get_grid_dim();
  n=n*n*n;
  for(int i=0; i<n; ++i) p->isgarr[i] = 1;
  return p;
}

// ElySphere.cc ends here
