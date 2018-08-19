/* The value of the dielectric constant on a cubic lattice. -*- C++ -*-

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

$Id: DielCubeRep.h,v 2.5 2007/05/28 01:26:42 bashford Exp $
*/
#ifndef _DielCubeRep_h
#define _DielCubeRep_h 1

#include "MEAD/CubeLatSpec.h"
#include "MEAD/globals.h"

class DielCubeRep {
public:
  DielCubeRep(const CubeLatSpec&);
  DielCubeRep(const DielCubeRep&);
  DielCubeRep& operator= (const DielCubeRep&);
  ~DielCubeRep();

  int eps_array_defined() const {return eps_def;}
  void declare_defined();
  const float& eps(int i, int j, int k) const;
  float operator[](int h) const;
  float& operator[](int h);
  float eps_neg_z(int i, int j, int k);
  float eps_pos_z(int i, int j, int k);
  float eps_neg_y(int i, int j, int k);
  float eps_pos_y(int i, int j, int k);
  float eps_neg_x(int i, int j, int k);
  float eps_pos_x(int i, int j, int k);
private:
  int grid_dim;
  int nsq;
  int ncube;
  int eps_def;
  struct DielCubeRepRep {
    unsigned int referenceCount;
    float *eps_array;
  };
  DielCubeRepRep *rep;
};

inline float DielCubeRep::operator[] (int h) const
{
  return rep->eps_array[h];
}

inline float& DielCubeRep::operator[] (int h)
{
  if (eps_def) {
    ::error("ERROR: illegal application of operator[] to a DielCubeRep",
	    "with eps_array defined");
  }
  return rep->eps_array[h];
}

inline const float& DielCubeRep::eps(int i, int j, int k) const
{
  return rep->eps_array[i*nsq + j*grid_dim + k];
}

inline float DielCubeRep::eps_neg_z(int i, int j, int k)
{
  return 4.0F / (1.0F/eps(i,j,k) + 1.0F/eps(i,j+1,k)
		+ 1.0F/eps(i+1,j,k) + 1.0F/eps(i+1,j+1,k));
}

inline float DielCubeRep::eps_pos_z(int i, int j, int k)
{
  return 4.0F / (1.0F/eps(i,j,k+1) + 1.0F/eps(i,j+1,k+1)
		+ 1.0F/eps(i+1,j,k+1) + 1.0F/eps(i+1,j+1,k+1));
}

inline float DielCubeRep::eps_neg_y(int i, int j, int k)
{
  return 4.0F / (1.0F/eps(i,j,k) + 1.0F/eps(i,j,k+1)
		+ 1.0F/eps(i+1,j,k) + 1.0F/eps(i+1,j,k+1));
}

inline float DielCubeRep::eps_pos_y(int i, int j, int k)
{
  return 4.0F / (1.0F/eps(i,j+1,k) + 1.0F/eps(i,j+1,k+1)
		+ 1.0F/eps(i+1,j+1,k) + 1.0F/eps(i+1,j+1,k+1));
}

inline float DielCubeRep::eps_neg_x(int i, int j, int k)
{
  return 4.0F / (1.0F/eps(i,j,k) + 1.0F/eps(i,j,k+1)
		+ 1.0F/eps(i,j+1,k) + 1.0F/eps(i,j+1,k+1));
}

inline float DielCubeRep::eps_pos_x(int i, int j, int k)
{
  return 4.0F / (1.0F/eps(i+1,j,k) + 1.0F/eps(i+1,j,k+1)
		+ 1.0F/eps(i+1,j+1,k) + 1.0F/eps(i+1,j+1,k+1));
}


#endif

// DielCubeRep.h ends here
