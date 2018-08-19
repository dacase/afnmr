/*
$Id: ChargeCubeRep.cc,v 2.4 2007/05/28 01:26:42 bashford Exp $
*/

#include "MEAD/globals.h"
#include "MEAD/ChargeCubeRep.h"
#include "MEAD/CubeLatSpec.h"

SparseChargeCubeRep::SparseChargeCubeRep(CubeLatSpec c, size_t num_charged_points)
{
  cls = c;
  int d = cls.get_grid_dim();
  lat_cube = d*d*d;
  if (num_charged_points > lat_cube)
    ::error("ERROR: SparseChargeCubeRep constructor: too many charged points");
  chrp = new ChargedPoint[lat_cube];
  last_added = last_indexed = chrp-1;
  highest_allowed = chrp + num_charged_points - 1;
  iscomplete = 0;
  is_sprs = 1;
}

SparseChargeCubeRep::~SparseChargeCubeRep()
{
  delete [] chrp;
}

// ChargeCubeRep.cc ends here
