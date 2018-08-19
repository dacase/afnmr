/* Successive over-relaxation to solve P--B eqn. by Fin. Diff.

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

$Id: sor.cc,v 2.6 2007/05/28 01:26:43 bashford Exp $
*/

#include "MEAD/FDGridLevel.h"
#include "MEAD/ChargeDist.h"
#include "MEAD/ChargeCubeRep.h"
#include "MEAD/FDChargeIterator.h"
#include "MEAD/NonUnif.h"
#include "MEAD/Bigmem.h"
#include <math.h>

bool FDGridLevel::use_fixed_maxrmsdiff = false;

void FDGridLevel::sor() // Do Successive overrelaxation.
{
  blab2 << "Entering FDGridLevel::sor." << endl;
  const float rjac = 0.989F;  // The Jacobi radius.

  // Some indices, offsets and pointers for speed and convenience

  // Starting points for runs through phi arrays for even calc.
  int e_init_index = (nsq+1)/2;
  float *e_init = phieven + e_init_index;    // The first even non i=0 point
  float *e_fence = phieven + (ncube-nsq)/2;  // The first even i=nside point
  float *ekappa_init = inv_six_plus_kappa_even + e_init_index;
  float *ekm_init = phiodd + e_init_index + ekm_off;
  float *ekp_init = phiodd + e_init_index + ekp_off;
  float *ejm_init = phiodd + e_init_index + ejm_off;
  float *ejp_init = phiodd + e_init_index + ejp_off;
  float *eim_init = phiodd + e_init_index + eim_off;
  float *eip_init = phiodd + e_init_index + eip_off;

  // Starting points for runs through phi arrays for odd calc.
  int o_init_index = (nsq+1)/2;
  float *o_init = phiodd + o_init_index;    // The first even non i=0 point
  float *o_fence = phiodd + (ncube-nsq+1)/2;  // The first odd i=nside point
  float *okappa_init = inv_six_plus_kappa_odd + o_init_index;
  float *okm_init = phieven + o_init_index + okm_off;
  float *okp_init = phieven + o_init_index + okp_off;
  float *ojm_init = phieven + o_init_index + ojm_off;
  float *ojp_init = phieven + o_init_index + ojp_off;
  float *oim_init = phieven + o_init_index + oim_off;
  float *oip_init = phieven + o_init_index + oip_off;

  NonUnif * nueven_fence = nueven+numnueven;
  NonUnif * nuodd_fence = nuodd+numnuodd;
  OuterBoundPoint *outbou_even_fence = phiboueven + num_outbou_even;
  OuterBoundPoint *outbou_odd_fence = phibouodd + num_outbou_odd;

  // The overrelaxation parameter.  USE CHEBYSHEV ACC?
  float omega = 1.0;

  // For convergence checks ...
  float *oldphi = (float*) big_rigid_malloc ((sizeof(float))*(ncube/2+1));

  int minits = (grid_dim-1)/2 * 3;
  int maxits = minits*10;
  const float maxrmsdiff = get_use_fixed_maxrmsdiff() ?
                              1.0e-5F
                            : 2.0e-5F / float(grid_dim-1);
  float rmsdiff = 0.0;
  float maxdiff = 0.0;
  for (int iteration = 0; (iteration < maxits &&
			   (iteration < minits || rmsdiff > maxrmsdiff));
       ++iteration) {

    int h, hfence;
    float sumsqdiff = 0.0;
    float maxdiffsq = 0.0;
    if (!(iteration%10)) {

      hfence = ncube/2 + 1;
      for (h = 0 ; h<hfence; ++h)
	oldphi[h] = phieven[h];
    }

    // Even
    // Uniform space calc for all points (this is really a for loop)
    float *phipt = e_init;       // The phi point being updated this iteration.
    float const *kappt = ekappa_init;  // Its 1/(6+spacing^2 * kappa^2) value.
    float const *km = ekm_init;        // Its i,j,k-1 neighbor.
    float const *kp = ekp_init;        // Its i,j,k+1 neighbor... etc.
    float const *jm = ejm_init;
    float const *jp = ejp_init;
    float const *im = eim_init;
    float const *ip = eip_init;

    while (phipt < e_fence) {
      float target = *kappt * (*km + *kp + *jm + *jp + *im + *ip);
      *phipt += omega * (target - *phipt);
      ++phipt;
      ++kappt;
      ++km; ++kp; ++jm; ++jp; ++im; ++ip;

    }

    // Correction for points with non-uniform epsilon
    NonUnif *pt;
    for (pt=nueven; pt<nueven_fence; ++pt) {
      float deltarget = pt->fac1 * *pt->ptr1
	+ pt->fac2 * *pt->ptr2
	  + pt->fac3 * *pt->ptr3
	    + pt->fac4 * *pt->ptr4
              + pt->fac5 * *pt->ptr5
		+ pt->fac6 * *pt->ptr6;
      *pt->phiptr += omega * deltarget;
    }

    // Correction for points with charges;
    fdchargeitr->fd_iterate_even(omega);

    // Correction for j and k boundary points wrongly included in
    // the uniform space calculation.
    OuterBoundPoint *outbou_ptr;
    for (outbou_ptr = phiboueven;
	 outbou_ptr<outbou_even_fence; ++outbou_ptr) {
      *outbou_ptr->phiptr = outbou_ptr->value;
    }

    if (!(iteration%10)) {
      for (h=0; h<hfence; ++h) {
	float diff = phieven[h] - oldphi[h];
	float diffsq = diff*diff;
	maxdiffsq = diffsq > maxdiffsq ? diffsq : maxdiffsq;
	sumsqdiff += diffsq;
      }
      hfence = ncube/2;
      for (h = 0 ; h<hfence; ++h)
	oldphi[h] = phiodd[h];
    }

    // Now do the same for odd...

    // Uniform space calc for all points (this is really a for loop)
    phipt = o_init;       // The phi point being updated this iteration.
    kappt = okappa_init;  // Its 1/(6+spacing^2 * kappa^2) value.
    km = okm_init;        // Its i,j,k-1 neighbor.
    kp = okp_init;        // Its i,j,k+1 neighbor... etc.
    jm = ojm_init;
    jp = ojp_init;
    im = oim_init;
    ip = oip_init;

    while (phipt < o_fence) {
      float target = *kappt * (*km + *kp + *jm + *jp + *im + *ip);
      *phipt += omega * (target - *phipt);
      ++phipt;
      ++kappt;
      ++km; ++kp; ++jm; ++jp; ++im; ++ip;
    }

    // Correction for points with non-uniform epsilon
    for (pt=nuodd; pt<nuodd_fence; ++pt) {
      float deltarget = pt->fac1 * *pt->ptr1
	+ pt->fac2 * *pt->ptr2
	  + pt->fac3 * *pt->ptr3
	    + pt->fac4 * *pt->ptr4
              + pt->fac5 * *pt->ptr5
		+ pt->fac6 * *pt->ptr6;
      *pt->phiptr += omega * deltarget;
    }

    // Correction for points with charges;
    fdchargeitr->fd_iterate_odd(omega);

    // Correction for j and k boundary points wrongly included in
    // the uniform space calculation.
    for (outbou_ptr = phibouodd; outbou_ptr<outbou_odd_fence;
	 ++outbou_ptr) {
      *outbou_ptr->phiptr = outbou_ptr->value;
    }

    if (!(iteration%10)) {
      for (h=0; h<hfence; ++h) {
	float diff = phiodd[h] - oldphi[h];
	float diffsq = diff*diff;
	maxdiffsq = diffsq > maxdiffsq ? diffsq : maxdiffsq;
	sumsqdiff += diffsq;
      }
      rmsdiff = sqrt (sumsqdiff / float(ncube));
      maxdiff = sqrt(maxdiffsq);
      blab2 << "iteration " << iteration
	<< " rms = " << rmsdiff << " maxdiff = " << maxdiff << endl;
    }

    // Chebyshev acceleration
    if (iteration ==1)
      omega = 1.0F / (1.0F - rjac*rjac/2.0F);
    else
      omega = 1.0F / (1.0F - rjac*rjac * omega / 4.0F);

  }
  if (rmsdiff > maxrmsdiff)
    blab2 << "PBEquation::sor WARNING: Didn't converge, rmsdiff > "
      << maxrmsdiff << endl;
  big_rigid_free (oldphi);
  blab2 << "Exiting FDGridLevel::sor." << endl;
}

// sor.cc ends here
