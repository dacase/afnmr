// This is -*- C++ -*-
/*
Iteration through charge lattice with hooks into potential lattice
   for finite-difference procedure.

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

$Id: FDChargeIterator.h,v 2.4 2004/12/06 17:58:20 bashford Exp $
*/
#ifndef _FDChargeIterator_h
#define _FDChargeIterator_h 1

#include <iostream>
#include "MEAD/CubeLatSpec.h"
#include "MEAD/ChargeCubeRep.h"
#include "MEAD/globals.h"

class FDChargeIterator {
public:
  FDChargeIterator(ChargeCubeRep* ccr, float* phieven, float* phiodd);
  virtual ~FDChargeIterator();
  virtual void complete_point(int h, float epsave, float inv_six_plus_kappa,
			      float *phieven, float*phiodd)
    {rep->complete_point(h, epsave, inv_six_plus_kappa, phieven, phiodd);}
  virtual void fd_iterate_odd (float omega) {rep->fd_iterate_odd(omega);}
  virtual void fd_iterate_even (float omega) {rep->fd_iterate_even(omega);}
  virtual ostream& print(ostream& ost) const {return rep->print(ost);}
protected:
  FDChargeIterator() {rep = 0;}
  FDChargeIterator *rep;
};

inline ostream & operator<< (ostream& ost, const FDChargeIterator& fdc)
{return fdc.print(ost);}

class SparseFDChargeIterator : public FDChargeIterator {
public:
// FIXME phieven and phiodd not needed in ctor?
  SparseFDChargeIterator(ChargeCubeRep* ccrp, float* phieven, float* phiodd);
  void complete_point(int h, float epsave, float inv_six_plus_kappa,
		      float *phieven, float*phiodd);
  virtual ~SparseFDChargeIterator();
  virtual void fd_iterate_odd (float omega);
  virtual void fd_iterate_even (float omega);
  virtual ostream& print(ostream& ost) const ;
private:
  CubeLatSpec cls;
  int numcheven, numchodd;
  struct ChargedPoint {
    union {
      int h;
      float *phiptr;
    };
    float charge_term;
  };
  ChargedPoint *cheven, *chodd, *nextevenp, *nextoddp, *higheven, *highodd;
  float *phieven, *phiodd;
  int nextevenh, nextoddh;
  float fourpi_over_spacing;
};

inline void SparseFDChargeIterator::complete_point(int h, float epsave,
					    float inv_six_plus_kappa,
					    float *phievenp, float*phioddp)
{
  if (phieven) {
    if (phievenp != phieven)
      ::error ("ERROR: SparseFDChargeIterator::complete_point: ",
	       "phieven does not match value from previous call\n");
  }
  else
    phieven = phievenp;
  if (phiodd) {
    if (phioddp != phiodd)
      ::error ("ERROR: SparseFDChargeIterator::complete_point: ",
	       "phiodd does not match value from previous call\n");
  }
  else
    phiodd = phioddp;

  if (h%2) {
    if (h==nextoddh) {
      nextoddp->phiptr = phiodd + h/2;
      nextoddp->charge_term *= (inv_six_plus_kappa * fourpi_over_spacing
				/ epsave);
      ++nextoddp;
      if (nextoddp <= highodd)
	nextoddh = nextoddp->h;
    }
  }
  else {
    if (h==nextevenh) {
      nextevenp->phiptr = phieven + h/2;
      nextevenp->charge_term *= (inv_six_plus_kappa * fourpi_over_spacing
				/ epsave);
      ++nextevenp;
      if (nextevenp <= higheven)
	nextevenh = nextevenp->h;
    }
  }
}


#endif

// FDChargeIterator.h ends here
