/* Dielectric on a cubic lattice.

    Copyright (c) 1993--1995 by Donald Bashford

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

$Id: DielCubeRep.cc,v 2.7 2008/11/06 20:38:51 bashford Exp $
*/
#include "MEAD/globals.h"
#include "MEAD/Bigmem.h"
#include "MEAD/DielCubeRep.h"
#include <iostream>
#include <cstring>
using std::memset;

DielCubeRep::DielCubeRep(const CubeLatSpec& cls)
{
  grid_dim = cls.get_grid_dim();
  nsq = grid_dim*grid_dim;
  ncube = grid_dim*grid_dim*grid_dim;
  eps_def =  0;
  rep = new DielCubeRepRep;
  rep->referenceCount = 1;
  rep->eps_array = (float *) big_rigid_malloc((sizeof (float))*ncube);
  memset ((char *) rep->eps_array, 0, ncube * sizeof (float));
}

DielCubeRep::DielCubeRep(const DielCubeRep& dcr)
{
  grid_dim = dcr.grid_dim;
  nsq = dcr.nsq;
  ncube = dcr.ncube;
  if (dcr.eps_def) { // The old one is not to change so safe to point to it.
    rep = dcr.rep;
    rep->referenceCount++;
  }
  else { // The old one or this one may change so make a copy
    rep = new DielCubeRepRep;
    rep->referenceCount = 1;
    rep->eps_array = new float [ncube];
    for (int h = 0; h<ncube; ++h)
      rep->eps_array[h] = dcr.rep->eps_array[h];
  }
  eps_def = dcr.eps_def;
}

DielCubeRep& DielCubeRep::operator=(const DielCubeRep& dcr)
{
  if (eps_def)
    ::error("Tryed to assign to a DieCubeRep with a fixed eps_array");
  else {
    if (rep->referenceCount != 1)
      ::error("ERROR: DielCubeRep::operator=: this object has eps_def",
	      " but referenceCount != 1");
    grid_dim = dcr.grid_dim;
    nsq = dcr.nsq;
    // The usual decrement and test of referenceCount is not needed
    // since we know it is 1.
    if (dcr.eps_def) { // The old one is not to change so safe to point to it.
      big_rigid_free(rep->eps_array);
      delete rep;
      rep = dcr.rep;
      rep->referenceCount++;
    }
    else { // The old one or this one may change so make a copy
      // Use the old space if it is the right size
      if(dcr.ncube != ncube) {
	big_rigid_free (rep->eps_array);
	rep->eps_array = new float [dcr.ncube];
      }
      for (int h = 0; h<dcr.ncube; ++h)
	rep->eps_array[h] = dcr.rep->eps_array[h];
    }
    ncube = dcr.ncube;
    eps_def = dcr.eps_def;
  }
  return *this;
}

void DielCubeRep::declare_defined()
{
  eps_def = 1;
}


DielCubeRep::~DielCubeRep()
{
  if(--rep->referenceCount == 0) {
    big_rigid_free(rep->eps_array);
    delete rep;
  }
}

// DielCubeRep.cc ends here
