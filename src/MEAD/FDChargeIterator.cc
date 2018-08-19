/* Iteration through charge lattice with hooks into potential lattice
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


$Id: FDChargeIterator.cc,v 2.3 2004/12/06 17:58:19 bashford Exp $
*/

#include "MEAD/globals.h"
#include "MEAD/FDChargeIterator.h"

// FIXME!!  The semantics of complete_point are more complex than
// they need to be since phieven and phiodd are already given here!

FDChargeIterator::FDChargeIterator(ChargeCubeRep* ccrp,
				    float* phieven, float* phiodd)
{
  if (ccrp->is_sparse())
    rep = new SparseFDChargeIterator(ccrp, phieven, phiodd);
  else ::error("Non-sparse case not implemented");
}

FDChargeIterator::~FDChargeIterator()
{
  if (rep) // then this is an envelope
    delete rep;
}

SparseFDChargeIterator::SparseFDChargeIterator(ChargeCubeRep* ccrp,
					       float* phievenp, float* phioddp)
{
  cls = ccrp->get_cubelatspec();
  // Do counts for even and odd arrays
  numcheven = numchodd = 0;
  int h;
  for(h = ccrp->firstindex(); ccrp->valid(h); h = ccrp->nextindex())
    h%2 ? ++numchodd : ++numcheven;
  // Create and fill even and odd arrays
  if (numcheven) {
    cheven = new ChargedPoint[numcheven];
    higheven = cheven + numcheven - 1;
  }
  else {
    cheven = 0;
    higheven = 0;
    nextevenh = -1;
  }
  if (numchodd) {
    chodd = new ChargedPoint[numchodd];
    highodd = chodd + numchodd - 1;
  }
  else {
    chodd = 0;
    highodd = 0;
    nextoddh = -1;
  }
  nextoddp = chodd;
  nextevenp = cheven;
  for(h = ccrp->firstindex(); ccrp->valid(h); h = ccrp->nextindex()) {
    if (h%2) {
      nextoddp->charge_term = ccrp->get(h);
      nextoddp->h = h;
      ++nextoddp;
    }
    else {
      nextevenp->charge_term = ccrp->get(h);
      nextevenp->h = h;
      ++nextevenp;
    }
  }
  if (chodd) nextoddh = chodd[0].h;
  nextoddp = chodd;
  if (cheven) nextevenh = cheven[0].h;
  nextevenp = cheven;
  fourpi_over_spacing = 4 * pi / (ccrp->get_cubelatspec()).get_spacing();
  phieven = phiodd = 0;
}

SparseFDChargeIterator::~SparseFDChargeIterator()
{
  if (cheven)
    delete [] cheven;
  if (chodd)
    delete [] chodd;
}

void SparseFDChargeIterator::fd_iterate_odd (float omega)
{
  if (chodd==0) return;
  if (nextoddp <= highodd)
    ::error("ERROR: SparseFDChargeIterator::fd_interate_odd invoked on",
	    "incomplete iterator");
  for (ChargedPoint *cpt = chodd; cpt<=highodd; ++cpt) {
    *cpt->phiptr += omega * cpt->charge_term;
  }
}

void SparseFDChargeIterator::fd_iterate_even (float omega)
{
  if (cheven==0) return;
  if (nextevenp <= higheven)
    ::error("ERROR: SparseFDChargeIterator::fd_interate_even invoked on",
	    "incomplete iterator");
  for (ChargedPoint *cpt = cheven; cpt<=higheven; ++cpt) {
    *cpt->phiptr += omega * cpt->charge_term;
  }
}

ostream & SparseFDChargeIterator::print (ostream& ost) const
{
  if (cheven==0)
    ost << "No even points" << endl;
  else {
    if (nextevenp <= higheven)
      ost << "The even points are incomplete" << endl;
    for (ChargedPoint *cpt = cheven; cpt<=higheven; ++cpt) {
      if (cpt < nextevenp)
	ost << (cpt->phiptr - phieven) << "  " << cpt->charge_term << endl;
      else
	ost << cpt->h << "  " << cpt->charge_term << " (incomplete)" << endl;
    }
  }
  if (chodd==0)
    ost << "No odd points" << endl;
  else {
    if (nextoddp <= highodd)
      ost << "The odd points are incomplete" << endl;
    for (ChargedPoint *cpt = chodd; cpt<=highodd; ++cpt) {
      if (cpt < nextoddp)
	ost << (cpt->phiptr - phiodd) << "  " << cpt->charge_term << endl;
      else
	ost << cpt->h << "  " << cpt->charge_term << " (incomplete)" << endl;
    }
  }
  return ost;
}

// FDChargeIterator.cc ends here
