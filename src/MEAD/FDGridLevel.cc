/* A level of grid spacing/size in finite difference procedure.

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

$Id: FDGridLevel.cc,v 2.16 2011/02/10 23:22:33 agiammon Exp $
*/

#include <vector>
using std::vector;
#include "MEAD/NonUnif.h"

#include "MEAD/ElstatPot.h"
#include "MEAD/FDGridLevel.h"
#include "MEAD/CubeLatSpec.h"
#include "MEAD/ChargeDist.h"
#include "MEAD/ChargeCubeRep.h"
#include "MEAD/FDChargeIterator.h"
#include "MEAD/DielectricEnvironment.h"
#include "MEAD/ElectrolyteEnvironment.h"
#include "MEAD/DielCubeRep.h"

#include "MEAD/PhysCond.h"
#include "MEAD/Bigmem.h"

bool FDGridLevel::epsave_oldway = false;
bool FDGridLevel::ZeroForOutOfRange = false;

FDGridLevel::FDGridLevel(const CubeLatSpec& s, ChargeDist_lett* r,
			 DielectricEnvironment_lett* ep,
			 ElectrolyteEnvironment_lett* el,
			 AnalyticEP *aapr)
{
  coarser = 0;
  rho = r;
  eps = ep;
  electrolyte = el;
  analytic_approx = aapr;
  center = s.get_center();
  grid_dim = s.get_grid_dim();
  spacing = s.get_spacing();
  lat = s;
  solved = 0;
  misc_settings();
}

// FIXME?  How can I avoid requiring an AnalyticEP for this ctor
// when AnalyticEP::AnalyticEP() is protected.  The FDGridLevel here
// doesnt need analytic_approx since in has the previous grid.
// Use pointers instead?  Make a public AnalyticEP::AnalyticEP()??

FDGridLevel::FDGridLevel(const CubeLatSpec& s, ChargeDist_lett* r,
			 DielectricEnvironment_lett* ep,
			 ElectrolyteEnvironment_lett* el,
			 FDGridLevel * coarser_fdl)
{
  coarser = coarser_fdl;
  rho = r;
  eps = ep;
  electrolyte = el;
  center = s.get_center();
  grid_dim = s.get_grid_dim();
  spacing = s.get_spacing();
  lat = s;
  solved = 0;
  misc_settings();
}

FDGridLevel::FDGridLevel (const CubeLatSpec& s, float *phi)
{
  coarser = 0;
  rho = 0;
  eps = 0;
  electrolyte = 0;
  analytic_approx = 0;
  center = s.get_center();
  grid_dim = s.get_grid_dim();
  spacing = s.get_spacing();
  lat = s;
  misc_settings();  // some of which are useless for this case.
  phieven = (float *) big_rigid_malloc ((sizeof (float))*(ncube/2 + 1));
  phiodd = (float *) big_rigid_malloc ((sizeof (float))*(ncube/2));
  float *evenp = phieven;
  float *oddp = phiodd;
  for (int h=0; h<ncube; ++h) {
    if (h%2)
      *oddp++ = phi[h];
    else
      *evenp++ = phi[h];
  }
  solved = 1;
}

FDGridLevel::FDGridLevel (const CubeLatSpec& s, float *phi,
			  DielectricEnvironment_lett* e)
{
  coarser = 0;
  rho = 0;
  eps = e;
  electrolyte = 0;
  analytic_approx = 0;
  center = s.get_center();
  grid_dim = s.get_grid_dim();
  spacing = s.get_spacing();
  lat = s;
  misc_settings();  // some of which are useless for this case.
  phieven = (float *) big_rigid_malloc ((sizeof (float))*(ncube/2 + 1));
  phiodd = (float *) big_rigid_malloc ((sizeof (float))*(ncube/2));
  float *evenp = phieven;
  float *oddp = phiodd;
  for (int h=0; h<ncube; ++h) {
    if (h%2)
      *oddp++ = phi[h];
    else
      *evenp++ = phi[h];
  }
  solved = 1;
}

FDGridLevel::FDGridLevel (const CubeLatSpec& s, float *phi, FDGridLevel* c)
{
  coarser = c;
  rho = 0;
  eps = 0;
  electrolyte = 0;
  analytic_approx = 0;
  center = s.get_center();
  grid_dim = s.get_grid_dim();
  spacing = s.get_spacing();
  lat = s;
  misc_settings();  // some of which are useless for this case.
  phieven = (float *) big_rigid_malloc ((sizeof (float))*(ncube/2 + 1));
  phiodd = (float *) big_rigid_malloc ((sizeof (float))*(ncube/2));
  float *evenp = phieven;
  float *oddp = phiodd;
  for (int h=0; h<ncube; ++h) {
    if (h%2)
      *oddp++ = phi[h];
    else
      *evenp++ = phi[h];
  }
  solved = 1;
}

FDGridLevel::FDGridLevel (const CubeLatSpec& s, float *phi,
			  DielectricEnvironment_lett* e, FDGridLevel* c)
{
  coarser = c;
  rho = 0;
  eps = e;
  electrolyte = 0;
  analytic_approx = 0;
  center = s.get_center();
  grid_dim = s.get_grid_dim();
  spacing = s.get_spacing();
  lat = s;
  misc_settings();  // some of which are useless for this case.
  phieven = (float *) big_rigid_malloc ((sizeof (float))*(ncube/2 + 1));
  phiodd = (float *) big_rigid_malloc ((sizeof (float))*(ncube/2));
  float *evenp = phieven;
  float *oddp = phiodd;
  for (int h=0; h<ncube; ++h) {
    if (h%2)
      *oddp++ = phi[h];
    else
      *evenp++ = phi[h];
  }
  solved = 1;
}



void FDGridLevel::misc_settings()
{
  // Zero data structues that are not yet in use.
  numnueven = numnuodd = num_outbou_even = num_outbou_odd = 0;
  phieven = phiodd = inv_six_plus_kappa_even = inv_six_plus_kappa_odd = 0;
  cubecharge = 0;
  fdchargeitr = 0;
  dcrp = 0;
  phiboueven = phibouodd = 0;

  // Set up various convenience varaibles
  nsq = grid_dim*grid_dim;
  ncube = nsq*grid_dim;
  // phieven[i]'s j-1 neighbor is phiodd[i+ejm_off] etc.
  ekm_off = -1;
  ekp_off = 0;
  ejm_off = -grid_dim/2 - 1;
  ejp_off = grid_dim/2;
  eim_off = -nsq/2 - 1;
  eip_off = nsq/2;
  // phiodd[i]'s j-1 neighbor is phieven[i+ojm_off] etc.
  okm_off = 0;
  okp_off = 1;
  ojm_off = -grid_dim/2;
  ojp_off = grid_dim/2 + 1;
  oim_off = -nsq/2;
  oip_off = nsq/2 + 1;

  float halfside = (grid_dim-1)/2 * spacing;
  xlo = center.x - halfside;
  xhi = center.x + halfside;
  ylo = center.y - halfside;
  yhi = center.y + halfside;
  zlo = center.z - halfside;
  zhi = center.z + halfside;
}

FDGridLevel::~FDGridLevel()
{
  // coarser?
  // rho, eps, electrolyte?
  big_rigid_free(phieven);
  big_rigid_free(phiodd);
  big_rigid_free(inv_six_plus_kappa_even);
  big_rigid_free(inv_six_plus_kappa_odd);
  delete [] nueven;
  delete [] nuodd;
  delete cubecharge;
  delete fdchargeitr;
  delete [] phiboueven;
  delete [] phibouodd;
  delete dcrp;
}

void FDGridLevel::solve()
{
  if (solved) {
    cerr << "WARNING: FDGridLevel::solve called for already solved object.\n"
      << "Returning without doing anything.\n" << endl;
    return;
  }
  blab2 << "FDGridLevel::solve called.  Getting cubic reps:" << endl;

  if (rho==0 || eps==0 || electrolyte==0)
    error("ERROR: FDGridLevel::solve: This object lacks one or more of the\n",
	  "DielectricEnvironment, ElectrolyteEnvironment or ChargeDist \n",
	  "needed to define a solution.\n");
  phieven = (float *) big_rigid_malloc ((sizeof (float))*(ncube/2 + 1));
  phiodd = (float *) big_rigid_malloc ((sizeof (float))*(ncube/2));
  num_outbou_even = (nsq/2 + 1)*4;
  num_outbou_odd = (nsq/2)*4;
  phiboueven = new OuterBoundPoint [num_outbou_even];
  phibouodd = new OuterBoundPoint [num_outbou_odd];

  if (coarser) {
    blab2 << "Initializing from coarser grid" << endl;
    phi_from_prev();
  }
  else {
    if (analytic_approx==0)
      error ("ERROR: FDGridLevel::solve: This object lacks analytical \n",
	     "approximation (AnalyticEP) for the boundary potential.\n");
    blab2 << "Initilizing phi by analytical formula" << endl;
    phi_init();
  }

// FIXME should have a special case if no salt.
  set_up_sor();
  sor();
  solved = 1;

  big_rigid_free (inv_six_plus_kappa_even);
  big_rigid_free (inv_six_plus_kappa_odd);
  delete [] nueven;
  delete [] nuodd;
  delete [] phiboueven;
  delete [] phibouodd;
  nueven = nuodd = 0;
  inv_six_plus_kappa_even = inv_six_plus_kappa_odd = 0;
  phiboueven = phibouodd = 0;
}

void FDGridLevel::solve_using_initial(const string& initial_fieldname)
{
  if (solved) {
    cerr << "WARNING: FDGridLevel::solve called for already solved object.\n"
      << "Returning without doing anything.\n" << endl;
    return;
  }
  blab2 << "FDGridLevel::solve called.  Getting cubic reps:" << endl;

  if (rho==0 || eps==0 || electrolyte==0)
    error("ERROR: FDGridLevel::solve: This object lacks one or more of the\n",
	  "DielectricEnvironment, ElectrolyteEnvironment or ChargeDist \n",
	  "needed to define a solution.\n");
  phieven = (float *) big_rigid_malloc ((sizeof (float))*(ncube/2 + 1));
  phiodd = (float *) big_rigid_malloc ((sizeof (float))*(ncube/2));
  num_outbou_even = (nsq/2 + 1)*4;
  num_outbou_odd = (nsq/2)*4;
  phiboueven = new OuterBoundPoint [num_outbou_even];
  phibouodd = new OuterBoundPoint [num_outbou_odd];

  blab2 << "Initializing from AVS field, " << initial_fieldname << "."
    << endl;
  if (read(lat, initial_fieldname)) {
    blab2 << "Reading in AVS field succeded." << endl;
    init_sides_from_phi();
  }
  else {
    blab2 << "Initialization from AVS field failed.  Will try analytic..."
      << endl;
    if (analytic_approx==0)
      error ("ERROR: FDGridLevel::solve: This object lacks analytical \n",
	     "approximation (AnalyticEP) for the boundary potential.\n");
    blab2 << "Initilizing phi by analytical formula" << endl;
    phi_init();
  }

// FIXME should have a special case if no salt.
  set_up_sor();
  sor();
  solved = 1;

  big_rigid_free (inv_six_plus_kappa_even);
  big_rigid_free (inv_six_plus_kappa_odd);
  delete [] nueven;
  delete [] nuodd;
  delete [] phiboueven;
  delete [] phibouodd;
  nueven = nuodd = 0;
  inv_six_plus_kappa_even = inv_six_plus_kappa_odd = 0;
  phiboueven = phibouodd = 0;
}

float FDGridLevel::potint (float xp, float yp, float zp) const
{
  if (!solved)
    ::error("FDGridLevel::potint called for unsolved object");
  if (xp>xlo && xp<xhi && yp>ylo && yp<yhi && zp>zlo && zp<zhi) {
    float halflen = ((float) (grid_dim - 1)) / 2.0F;
    float x = (xp - center.x)/spacing + halflen;
    float y = (yp - center.y)/spacing + halflen;
    float z = (zp - center.z)/spacing + halflen;
    int i0 = (int) x;
    // Avoid going out of array bounds for boundary point:
    if (i0 == grid_dim) --i0;
    int i1 =  i0 + 1;
    int j0 = (int) y;
    if (j0 == grid_dim) --j0;
    int j1 = j0 + 1;
    int k0 = (int) z;
    if (k0 == grid_dim) --k0;
    int k1 = k0 + 1;
    float xr = x - (float) i0;
    float yr = y - (float) j0;
    float zr = z - (float) k0;
    float phi000 = phiof_ijk(i0,j0,k0);
    float phi001 = phiof_ijk(i0,j0,k1);
    float phi010 = phiof_ijk(i0,j1,k0);
    float phi011 = phiof_ijk(i0,j1,k1);
    float phi100 = phiof_ijk(i1,j0,k0);
    float phi101 = phiof_ijk(i1,j0,k1);
    float phi110 = phiof_ijk(i1,j1,k0);
    float phi111 = phiof_ijk(i1,j1,k1);
    return (1.0F-xr) * ( (1.0F-yr) * ( (1.0F-zr)*phi000 + zr*phi001 )
		       + yr * ( (1.0F-zr)*phi010 + zr*phi011 ) )
      + xr * ( (1.0F-yr) * ( (1.0F-zr)*phi100 + zr*phi101 )
	      + yr * ( (1.0F-zr)*phi110 + zr*phi111 ) );
  }
  else if (!coarser){
    if (get_ZeroForOutOfRange())
      return 0.0;
    else 
      error ("potint: point out of range\n");
  }
  return coarser->potint(xp,yp,zp);
}

// Return the electric field vector at point p.
Coord FDGridLevel::fieldint (const Coord& p) const
{
  if (!solved)
    ::error("FDGridLevel::potint called for unsolved object");
  // Translate p into lattice coordinates.
  float len = ((float) (grid_dim - 1));
  float halflen = len / 2.0F;
  float x = (p.x - center.x)/spacing + halflen;
  float y = (p.y - center.y)/spacing + halflen;
  float z = (p.z - center.z)/spacing + halflen;
  float toplim = len - 0.5F;
  if (x>0.5 && x<toplim && y>0.5 && y<toplim && z>0.5 && z<toplim) {
    // Indices for surrounding phi points
    int i0 = (int) x;
    int i1 = i0 + 1;
    int j0 = (int) y;
    int j1 = j0 + 1;
    int k0 = (int) z;
    int k1 = k0 + 1;
    // X component of the field:
    // The lattice defines the X component of the field at points along
    // a line parallel to the X-axis midway between phi points.  So find
    // The eight nearest such points and interpolate between them.
    // These will be points half integer in X and integer and Y and Z.
    int ic = i0; // possible central phi lattice plane just behind us.
    float xr = x - (float) ic + 0.5F;  // Distance from rear midway plane
    if (xr > 1.0) { // Too far!  Next one should be used.
      xr -= 1.0;
      ++ic;
    }
    int ip =  ic + 1;
    int im =  ic - 1;
    // Integer values for y and z interpolation points
    float yr = y - (float) j0;
    float zr = z - (float) k0;

    // The values to interpolate:
    float e000 = phiof_ijk(im,j0,k0) - phiof_ijk(ic,j0,k0);
    float e001 = phiof_ijk(im,j0,k1) - phiof_ijk(ic,j0,k1);
    float e010 = phiof_ijk(im,j1,k0) - phiof_ijk(ic,j1,k0);
    float e011 = phiof_ijk(im,j1,k1) - phiof_ijk(ic,j1,k1);
    float e100 = phiof_ijk(ic,j0,k0) - phiof_ijk(ip,j0,k0);
    float e101 = phiof_ijk(ic,j0,k1) - phiof_ijk(ip,j0,k1);
    float e110 = phiof_ijk(ic,j1,k0) - phiof_ijk(ip,j1,k0);
    float e111 = phiof_ijk(ic,j1,k1) - phiof_ijk(ip,j1,k1);
    // Interpolation for the X component of the field:
    float ex = (1.0F-xr) * ( (1.0F-yr) * ( (1.0F-zr)*e000 + zr*e001 )
		       + yr * ( (1.0F-zr)*e010 + zr*e011 ) )
      + xr * ( (1.0F-yr) * ( (1.0F-zr)*e100 + zr*e101 )
	      + yr * ( (1.0F-zr)*e110 + zr*e111 ) );

    // Y component:
    int jc = j0;
    yr = y - (float) jc + 0.5F;
    if (yr > 1.0) {
      yr -= 1.0;
      ++jc;
    }
    int jp =  jc + 1;
    int jm =  jc - 1;
    xr = x - (float) i0;
    zr = z - (float) k0;
    e000 = phiof_ijk(i0,jm,k0) - phiof_ijk(i0,jc,k0);
    e001 = phiof_ijk(i0,jm,k1) - phiof_ijk(i0,jc,k1);
    e010 = phiof_ijk(i0,jc,k0) - phiof_ijk(i0,jp,k0);
    e011 = phiof_ijk(i0,jc,k1) - phiof_ijk(i0,jp,k1);
    e100 = phiof_ijk(i1,jm,k0) - phiof_ijk(i1,jc,k0);
    e101 = phiof_ijk(i1,jm,k1) - phiof_ijk(i1,jc,k1);
    e110 = phiof_ijk(i1,jc,k0) - phiof_ijk(i1,jp,k0);
    e111 = phiof_ijk(i1,jc,k1) - phiof_ijk(i1,jp,k1);
    // Interpolation for the X component of the field:
    float ey = (1.0F-xr) * ( (1.0F-yr) * ( (1.0F-zr)*e000 + zr*e001 )
		       + yr * ( (1.0F-zr)*e010 + zr*e011 ) )
      + xr * ( (1.0F-yr) * ( (1.0F-zr)*e100 + zr*e101 )
	      + yr * ( (1.0F-zr)*e110 + zr*e111 ) );

    // Z component:
    int kc = k0;
    zr = z - (float) kc + 0.5F;
    if (zr > 1.0) {
      zr -= 1.0;
      ++kc;
    }
    int kp =  kc + 1;
    int km =  kc - 1;
    xr = x - (float) i0;
    yr = y - (float) j0;
    e000 = phiof_ijk(i0,j0,km) - phiof_ijk(i0,j0,kc);
    e001 = phiof_ijk(i0,j0,kc) - phiof_ijk(i0,j0,kp);
    e010 = phiof_ijk(i0,j1,km) - phiof_ijk(i0,j1,kc);
    e011 = phiof_ijk(i0,j1,kc) - phiof_ijk(i0,j1,kp);
    e100 = phiof_ijk(i1,j0,km) - phiof_ijk(i1,j0,kc);
    e101 = phiof_ijk(i1,j0,kc) - phiof_ijk(i1,j0,kp);
    e110 = phiof_ijk(i1,j1,km) - phiof_ijk(i1,j1,kc);
    e111 = phiof_ijk(i1,j1,kc) - phiof_ijk(i1,j1,kp);
    float ez = (1.0F-xr) * ( (1.0F-yr) * ( (1.0F-zr)*e000 + zr*e001 )
		       + yr * ( (1.0F-zr)*e010 + zr*e011 ) )
      + xr * ( (1.0F-yr) * ( (1.0F-zr)*e100 + zr*e101 )
	      + yr * ( (1.0F-zr)*e110 + zr*e111 ) );

    return Coord(ex, ey, ez)/spacing;
  }

  else if (!coarser)
    error ("fieldint: point out of range\n");
  return coarser->fieldint(p);
}

// Return the electric displacement vector at point p.
Coord FDGridLevel::displacement_int (const Coord& p)
{
  if (!solved)
    ::error("FDGridLevel::potint called for unsolved object");
  // Translate p into lattice coordinates.
  float len = ((float) (grid_dim - 1));
  float halflen = len / 2.0F;
  float x = (p.x - center.x)/spacing + halflen;
  float y = (p.y - center.y)/spacing + halflen;
  float z = (p.z - center.z)/spacing + halflen;
  float toplim = len - 0.5F;
  if (x>0.5 && x<toplim && y>0.5 && y<toplim && z>0.5 && z<toplim) {
    // The DielCubeRep must be present and ready.  Make sure it is.
    if (!dcrp) {
      Coord half_offset (spacing/2.0F, spacing/2.0F, spacing/2.0F);
      CubeLatSpec offset_lattice(grid_dim, spacing, center - half_offset);
      dcrp = new DielCubeRep(offset_lattice);
    }
    if (!dcrp->eps_array_defined()) {
      if (eps == 0)
	::error("ERROR: FDGridLevel::displacement_int called for an",
		"object with no defined DielectricEnviroment\n",
		"and thus, no defined electric displacement.\n");

      Coord half_offset (spacing/2.0F, spacing/2.0F, spacing/2.0F);
      CubeLatSpec offset_lattice(grid_dim, spacing, center - half_offset);
      *dcrp = eps->get_cuberep(offset_lattice);
    }
    // Indices for surrounding phi points
    int i0 = (int) x;
    int i1 = i0 + 1;
    int j0 = (int) y;
    int j1 = j0 + 1;
    int k0 = (int) z;
    int k1 = k0 + 1;
    // X component of the displacement:
    // The lattice defines the X component of the displacement at points along
    // a line parallel to the X-axis midway between phi points.  So find
    // The eight nearest such points and interpolate between them.
    // These will be points half integer in X and integer and Y and Z.
    int ic = i0; // possible central phi lattice plane just behind us.
    float xr = x - (float) ic + 0.5F;  // Distance from rear midway plane
    if (xr > 1.0) { // Too far!  Next one should be used.
      xr -= 1.0;
      ++ic;
    }
    int ip =  ic + 1;
    int im =  ic - 1;
    // Integer values for y and z interpolation points
    float yr = y - (float) j0;
    float zr = z - (float) k0;



    // The values to interpolate:
    float e000 = (phiof_ijk(im,j0,k0) - phiof_ijk(ic,j0,k0)) * dcrp->eps_neg_x(ic,j0,k0);
    float e001 = (phiof_ijk(im,j0,k1) - phiof_ijk(ic,j0,k1)) * dcrp->eps_neg_x(ic,j0,k1);
    float e010 = (phiof_ijk(im,j1,k0) - phiof_ijk(ic,j1,k0)) * dcrp->eps_neg_x(ic,j1,k0);
    float e011 = (phiof_ijk(im,j1,k1) - phiof_ijk(ic,j1,k1)) * dcrp->eps_neg_x(ic,j1,k1);
    float e100 = (phiof_ijk(ic,j0,k0) - phiof_ijk(ip,j0,k0)) * dcrp->eps_neg_x(ip,j0,k0);
    float e101 = (phiof_ijk(ic,j0,k1) - phiof_ijk(ip,j0,k1)) * dcrp->eps_neg_x(ip,j0,k1);
    float e110 = (phiof_ijk(ic,j1,k0) - phiof_ijk(ip,j1,k0)) * dcrp->eps_neg_x(ip,j1,k0);
    float e111 = (phiof_ijk(ic,j1,k1) - phiof_ijk(ip,j1,k1)) * dcrp->eps_neg_x(ip,j1,k1);
    // Interpolation for the X component of the displacement:
    float ex = (1.0F-xr) * ( (1.0F-yr) * ( (1.0F-zr)*e000 + zr*e001 )
		       + yr * ( (1.0F-zr)*e010 + zr*e011 ) )
      + xr * ( (1.0F-yr) * ( (1.0F-zr)*e100 + zr*e101 )
	      + yr * ( (1.0F-zr)*e110 + zr*e111 ) );

    // Y component:
    int jc = j0;
    yr = y - (float) jc + 0.5F;
    if (yr > 1.0) {
      yr -= 1.0;
      ++jc;
    }
    int jp =  jc + 1;
    int jm =  jc - 1;
    xr = x - (float) i0;
    zr = z - (float) k0;
    e000 = (phiof_ijk(i0,jm,k0) - phiof_ijk(i0,jc,k0)) * dcrp->eps_neg_y(i0,jc,k0);
    e001 = (phiof_ijk(i0,jm,k1) - phiof_ijk(i0,jc,k1)) * dcrp->eps_neg_y(i0,jc,k1);
    e010 = (phiof_ijk(i0,jc,k0) - phiof_ijk(i0,jp,k0)) * dcrp->eps_neg_y(i0,jp,k0);
    e011 = (phiof_ijk(i0,jc,k1) - phiof_ijk(i0,jp,k1)) * dcrp->eps_neg_y(i0,jp,k1);
    e100 = (phiof_ijk(i1,jm,k0) - phiof_ijk(i1,jc,k0)) * dcrp->eps_neg_y(i1,jc,k0);
    e101 = (phiof_ijk(i1,jm,k1) - phiof_ijk(i1,jc,k1)) * dcrp->eps_neg_y(i1,jc,k1);
    e110 = (phiof_ijk(i1,jc,k0) - phiof_ijk(i1,jp,k0)) * dcrp->eps_neg_y(i1,jp,k0);
    e111 = (phiof_ijk(i1,jc,k1) - phiof_ijk(i1,jp,k1)) * dcrp->eps_neg_y(i1,jp,k1);
    // Interpolation for the X component of the displacement:
    float ey = (1.0F-xr) * ( (1.0F-yr) * ( (1.0F-zr)*e000 + zr*e001 )
		       + yr * ( (1.0F-zr)*e010 + zr*e011 ) )
      + xr * ( (1.0F-yr) * ( (1.0F-zr)*e100 + zr*e101 )
	      + yr * ( (1.0F-zr)*e110 + zr*e111 ) );

    // Z component:
    int kc = k0;
    zr = z - (float) kc + 0.5F;
    if (zr > 1.0) {
      zr -= 1.0;
      ++kc;
    }
    int kp =  kc + 1;
    int km =  kc - 1;
    xr = x - (float) i0;
    yr = y - (float) j0;
    e000 = (phiof_ijk(i0,j0,km) - phiof_ijk(i0,j0,kc)) * dcrp->eps_neg_z(i0,j0,kc);
    e001 = (phiof_ijk(i0,j0,kc) - phiof_ijk(i0,j0,kp)) * dcrp->eps_neg_z(i0,j0,kp);
    e010 = (phiof_ijk(i0,j1,km) - phiof_ijk(i0,j1,kc)) * dcrp->eps_neg_z(i0,j1,kc);
    e011 = (phiof_ijk(i0,j1,kc) - phiof_ijk(i0,j1,kp)) * dcrp->eps_neg_z(i0,j1,kp);
    e100 = (phiof_ijk(i1,j0,km) - phiof_ijk(i1,j0,kc)) * dcrp->eps_neg_z(i1,j0,kc);
    e101 = (phiof_ijk(i1,j0,kc) - phiof_ijk(i1,j0,kp)) * dcrp->eps_neg_z(i1,j0,kp);
    e110 = (phiof_ijk(i1,j1,km) - phiof_ijk(i1,j1,kc)) * dcrp->eps_neg_z(i1,j1,kc);
    e111 = (phiof_ijk(i1,j1,kc) - phiof_ijk(i1,j1,kp)) * dcrp->eps_neg_z(i1,j1,kp);
    float ez = (1.0F-xr) * ( (1.0F-yr) * ( (1.0F-zr)*e000 + zr*e001 )
		       + yr * ( (1.0F-zr)*e010 + zr*e011 ) )
      + xr * ( (1.0F-yr) * ( (1.0F-zr)*e100 + zr*e101 )
	      + yr * ( (1.0F-zr)*e110 + zr*e111 ) );

    return Coord(ex, ey, ez)/spacing;
  }

  else if (!coarser)
    error ("displacement_int: point out of range\n");
  return coarser->displacement_int(p);
}


// Derived from the old PBEquation::phi_init, but departs from it in
// that potentials are to be calculated directly from the "native"
// representation of the charges contained in the AnalyticEP rather
// than from the ChargeCubeRep.
// Although this is more general it may be slow since AnalyticEP:value()
// must do sqrt on all grid points, etc.
// It might be more efficient if this function where preformed by AnalyticEP,
// but that means putting grid capability into AnalyticEPs.  Not clean!

void FDGridLevel::phi_init()
{
  blab2 << "Entering FDGridLevel::phi_init" << endl;

  int halfh;
  for (halfh=0; halfh<ncube/2; ++halfh)
    phiodd[halfh] = phieven[halfh] = 0.0;
  phieven[halfh] = 0.0;  // one more even point than odd.

  float grlen = (float) (grid_dim - 1);
  float spacelen = grlen * spacing;
  float halfgrlen = grlen/2;
  Coord gridcenter_in_grid(halfgrlen, halfgrlen, halfgrlen);
  Coord& gridcenter_in_space = center;
  Coord sidepoint_in_space;
  analytic_approx->solve();

  // j=0, j=grid_dim-1 slabs
  int oddbouind = 0;
  int evenbouind = 0;
  int nsq_minus_grid_dim_over_2 = (nsq - grid_dim)/2;
  int i,j,k;
  for (i=0; i<grid_dim; ++i) {
    sidepoint_in_space.x =
      spacing * ((float) i - gridcenter_in_grid.x)
      + gridcenter_in_space.x;
    int iterm = i*nsq;
    for (k=0; k<grid_dim; ++k) {
      int h = iterm + k;
      sidepoint_in_space.z = spacing * ((float) k - gridcenter_in_grid.z)
	+ gridcenter_in_space.z;
      sidepoint_in_space.y =  spacing * ( - gridcenter_in_grid.y)
	+ gridcenter_in_space.y;
      if (h%2) { // odd case
	float *phiptr = phiodd + h/2;
	phibouodd[oddbouind].phiptr = phiptr;
	phibouodd[oddbouind].value = analytic_approx->value(sidepoint_in_space);
	++oddbouind;
        sidepoint_in_space.y += spacelen;
	phiptr += nsq_minus_grid_dim_over_2;
	phibouodd[oddbouind].phiptr = phiptr;
	phibouodd[oddbouind].value = analytic_approx->value(sidepoint_in_space);
	++oddbouind;
      }
      else { //even case
	float *phiptr = phieven + h/2;
	phiboueven[evenbouind].phiptr = phiptr;
	phiboueven[evenbouind].value = analytic_approx->value(sidepoint_in_space);
        sidepoint_in_space.y += spacelen;
	++evenbouind;
	phiptr += nsq_minus_grid_dim_over_2;
	phiboueven[evenbouind].phiptr = phiptr;
	phiboueven[evenbouind].value = analytic_approx->value(sidepoint_in_space);
	++evenbouind;
      }
    }
  }

  // k=0, k=grid_dim-1 slabs
  int grid_dim_minus_1_over_2 = (grid_dim - 1)/2;
  for (i=0; i<grid_dim; ++i) {
    sidepoint_in_space.x =
      spacing * ((float) i - gridcenter_in_grid.x)
      + gridcenter_in_space.x;
    int iterm = i*nsq;
    for (j=0; j<grid_dim; ++j) {
      sidepoint_in_space.y = spacing * ((float) j - gridcenter_in_grid.y)
	+ gridcenter_in_space.y;
      sidepoint_in_space.z =  spacing * ( - gridcenter_in_grid.z)
	+ gridcenter_in_space.z;
      int h = iterm + j*grid_dim;
      if (h%2) { // odd case
	float *phiptr = phiodd + h/2;
	phibouodd[oddbouind].phiptr = phiptr;
	phibouodd[oddbouind].value = analytic_approx->value(sidepoint_in_space);
	++oddbouind;
        sidepoint_in_space.z += spacelen;
	phiptr += grid_dim_minus_1_over_2;
	phibouodd[oddbouind].phiptr = phiptr;
	phibouodd[oddbouind].value = analytic_approx->value(sidepoint_in_space);
	++oddbouind;
      }
      else { //even case
	float *phiptr = phieven + h/2;
	phiboueven[evenbouind].phiptr = phiptr;
	phiboueven[evenbouind].value = analytic_approx->value(sidepoint_in_space);
	++evenbouind;
        sidepoint_in_space.z += spacelen;
	phiptr += grid_dim_minus_1_over_2;
	phiboueven[evenbouind].phiptr = phiptr;
	phiboueven[evenbouind].value = analytic_approx->value(sidepoint_in_space);
	++evenbouind;
      }
    }
  }

// The i=0, i=ngrid slabs
  int ncube_minus_nsq_over2 = (ncube - nsq) / 2;
  for (j=0; j<grid_dim; ++j) {
    sidepoint_in_space.y = spacing * ((float) j - gridcenter_in_grid.y)
      + gridcenter_in_space.y;
    int jterm = j*grid_dim;
    for (k=0; k<grid_dim; ++k) {
      int h = jterm + k;
      sidepoint_in_space.z = spacing * ((float) k - gridcenter_in_grid.z)
	+ gridcenter_in_space.z;
// FIXME here and similarly above, the r.h.s could be calculated outside loop
      sidepoint_in_space.x =  spacing * ( - gridcenter_in_grid.x)
	+ gridcenter_in_space.x;
      if (h%2) { // odd case
	phiodd[h/2] = analytic_approx->value(sidepoint_in_space);
        sidepoint_in_space.x += spacelen;
	phiodd[h/2+ncube_minus_nsq_over2]
	  = analytic_approx->value(sidepoint_in_space);
      }
      else { //even case
	phieven[h/2] = analytic_approx->value(sidepoint_in_space);
        sidepoint_in_space.x += spacelen;
	phieven[h/2+ncube_minus_nsq_over2]
	  = analytic_approx->value(sidepoint_in_space);
      }
    }
  }

  for (int iob=0; iob<num_outbou_odd; ++iob) {
    *phibouodd[iob].phiptr = phibouodd[iob].value;
  }
  for (int ieb=0; ieb<num_outbou_even; ++ieb) {
    *phiboueven[ieb].phiptr = phiboueven[ieb].value;
  }


  blab2 << "Exiting FDGridLevel::phi_init" << endl;
}



void FDGridLevel::phi_from_prev()
{
  blab2 << "Entering FDGridLevel::phi_from_prev" << endl;


  float grlen = (float) (grid_dim - 1);
  float spacelen = grlen * spacing;
  float halfgrlen = grlen/2;
  Coord gridcenter_in_grid(halfgrlen, halfgrlen, halfgrlen);
  Coord& gridcenter_in_space = center;


  if (1) { // initialize ENTIRE grid (not just sides) from the previous grid.

    // FIXME!  This is time consuming and it is not clear that it is
    // necessary (although it affects SOR convergence and things like
    // background interactions sligtly).  Perhaps it should be an option
    // controlled by a static data member setable from the command line.

    // FIXME This is slower than it needs to be since we use real space
    // coords as an intermediate between coarse and fine grid coords,
    // and because the bound checking of potint is a waste of time if
    // we already know that this grid is fully inside the previous one.

    Coord gridpoint_in_space;
    for (int i=0; i<grid_dim; ++i) {
      gridpoint_in_space.x =
	spacing * ((float) i - gridcenter_in_grid.x)
	  + gridcenter_in_space.x;
      int iterm = i*nsq;
      for (int j=0; j<grid_dim; ++j) {
	gridpoint_in_space.y = spacing * ((float) j - gridcenter_in_grid.y)
	  + gridcenter_in_space.y;
	int jterm = j*grid_dim;
	for (int k=0; k<grid_dim; ++k) {
	  gridpoint_in_space.z = spacing * ((float) k - gridcenter_in_grid.z)
	    + gridcenter_in_space.z;
	  int h = iterm + jterm + k;
	  if (h%2)
	    phiodd[h/2] = coarser->potint(gridpoint_in_space);
	  else
	    phieven[h/2] = coarser->potint(gridpoint_in_space);
	}
      }
    }
  }
  else {
    int halfh;
    for (halfh=0; halfh<ncube/2; ++halfh)
      phiodd[halfh] = phieven[halfh] = 0.0;
    phieven[halfh] = 0.0;  // one more even point than odd.
  }

// Set up phiboueven and phibouodd.
// FIXME.  Now should be able to use init_sides_from_phi() instead of
// this.  It will be faster.
  Coord sidepoint_in_space;

  // j=0, j=grid_dim-1 slabs
  int oddbouind = 0;
  int evenbouind = 0;
  int nsq_minus_grid_dim_over_2 = (nsq - grid_dim)/2;
  int i,j,k;
  for (i=0; i<grid_dim; ++i) {
    sidepoint_in_space.x =
      spacing * ((float) i - gridcenter_in_grid.x)
      + gridcenter_in_space.x;
    int iterm = i*nsq;
    for (k=0; k<grid_dim; ++k) {
      int h = iterm + k;
      sidepoint_in_space.z = spacing * ((float) k - gridcenter_in_grid.z)
	+ gridcenter_in_space.z;
      sidepoint_in_space.y =  spacing * ( - gridcenter_in_grid.y)
	+ gridcenter_in_space.y;
      if (h%2) { // odd case
	float *phiptr = phiodd + h/2;
	phibouodd[oddbouind].phiptr = phiptr;
	phibouodd[oddbouind].value = coarser->potint(sidepoint_in_space);
	++oddbouind;
        sidepoint_in_space.y += spacelen;
	phiptr += nsq_minus_grid_dim_over_2;
	phibouodd[oddbouind].phiptr = phiptr;
	phibouodd[oddbouind].value = coarser->potint(sidepoint_in_space);
	++oddbouind;
      }
      else { //even case
	float *phiptr = phieven + h/2;
	phiboueven[evenbouind].phiptr = phiptr;
	phiboueven[evenbouind].value = coarser->potint(sidepoint_in_space);
        sidepoint_in_space.y += spacelen;
	++evenbouind;
	phiptr += nsq_minus_grid_dim_over_2;
	phiboueven[evenbouind].phiptr = phiptr;
	phiboueven[evenbouind].value = coarser->potint(sidepoint_in_space);
	++evenbouind;
      }
    }
  }

  // k=0, k=grid_dim-1 slabs
  int grid_dim_minus_1_over_2 = (grid_dim - 1)/2;
  for (i=0; i<grid_dim; ++i) {
    sidepoint_in_space.x =
      spacing * ((float) i - gridcenter_in_grid.x)
      + gridcenter_in_space.x;
    int iterm = i*nsq;
    for (j=0; j<grid_dim; ++j) {
      sidepoint_in_space.y = spacing * ((float) j - gridcenter_in_grid.y)
	+ gridcenter_in_space.y;
      sidepoint_in_space.z =  spacing * ( - gridcenter_in_grid.z)
	+ gridcenter_in_space.z;
      int h = iterm + j*grid_dim;
      if (h%2) { // odd case
	float *phiptr = phiodd + h/2;
	phibouodd[oddbouind].phiptr = phiptr;
	phibouodd[oddbouind].value = coarser->potint(sidepoint_in_space);
	++oddbouind;
        sidepoint_in_space.z += spacelen;
	phiptr += grid_dim_minus_1_over_2;
	phibouodd[oddbouind].phiptr = phiptr;
	phibouodd[oddbouind].value = coarser->potint(sidepoint_in_space);
	++oddbouind;
      }
      else { //even case
	float *phiptr = phieven + h/2;
	phiboueven[evenbouind].phiptr = phiptr;
	phiboueven[evenbouind].value = coarser->potint(sidepoint_in_space);
	++evenbouind;
        sidepoint_in_space.z += spacelen;
	phiptr += grid_dim_minus_1_over_2;
	phiboueven[evenbouind].phiptr = phiptr;
	phiboueven[evenbouind].value = coarser->potint(sidepoint_in_space);
	++evenbouind;
      }
    }
  }

// FIXME?? Is this redundant since all of phi already initalized??
// The i=0, i=ngrid slabs
  int ncube_minus_nsq_over2 = (ncube - nsq) / 2;
  for (j=0; j<grid_dim; ++j) {
    sidepoint_in_space.y = spacing * ((float) j - gridcenter_in_grid.y)
      + gridcenter_in_space.y;
    int jterm = j*grid_dim;
    for (k=0; k<grid_dim; ++k) {
      int h = jterm + k;
      sidepoint_in_space.z = spacing * ((float) k - gridcenter_in_grid.z)
	+ gridcenter_in_space.z;
// FIXME here and similarly above, the r.h.s could be calculated outside loop
      sidepoint_in_space.x =  spacing * ( - gridcenter_in_grid.x)
	+ gridcenter_in_space.x;
      if (h%2) { // odd case
	phiodd[h/2] = coarser->potint(sidepoint_in_space);
        sidepoint_in_space.x += spacelen;
	phiodd[h/2+ncube_minus_nsq_over2]
	  = coarser->potint(sidepoint_in_space);
      }
      else { //even case
	phieven[h/2] = coarser->potint(sidepoint_in_space);
        sidepoint_in_space.x += spacelen;
	phieven[h/2+ncube_minus_nsq_over2]
	  = coarser->potint(sidepoint_in_space);
      }
    }
  }


  blab2 << "Exiting FDGridLevel::phi_from_prev" << endl;
}


void FDGridLevel::init_sides_from_phi()
{
// Set up phiboueven and phibouodd.

  Coord sidepoint_in_space;

  // j=0, j=grid_dim-1 slabs
  int oddbouind = 0;
  int evenbouind = 0;
  int nsq_minus_grid_dim_over_2 = (nsq - grid_dim)/2;
  int i,j,k;
  for (i=0; i<grid_dim; ++i) {
    int iterm = i*nsq;
    for (k=0; k<grid_dim; ++k) {
      int h = iterm + k;
      if (h%2) { // odd case
	float *phiptr = phiodd + h/2;
	phibouodd[oddbouind].phiptr = phiptr;
	phibouodd[oddbouind].value = *phiptr;
	++oddbouind;
	phiptr += nsq_minus_grid_dim_over_2;
	phibouodd[oddbouind].phiptr = phiptr;
	phibouodd[oddbouind].value = *phiptr;
	++oddbouind;
      }
      else { //even case
	float *phiptr = phieven + h/2;
	phiboueven[evenbouind].phiptr = phiptr;
	phiboueven[evenbouind].value = *phiptr;
	++evenbouind;
	phiptr += nsq_minus_grid_dim_over_2;
	phiboueven[evenbouind].phiptr = phiptr;
	phiboueven[evenbouind].value = *phiptr;
	++evenbouind;
      }
    }
  }

  // k=0, k=grid_dim-1 slabs
  int grid_dim_minus_1_over_2 = (grid_dim - 1)/2;
  for (i=0; i<grid_dim; ++i) {
    int iterm = i*nsq;
    for (j=0; j<grid_dim; ++j) {
      int h = iterm + j*grid_dim;
      if (h%2) { // odd case
	float *phiptr = phiodd + h/2;
	phibouodd[oddbouind].phiptr = phiptr;
	phibouodd[oddbouind].value = *phiptr;
	++oddbouind;
	phiptr += grid_dim_minus_1_over_2;
	phibouodd[oddbouind].phiptr = phiptr;
	phibouodd[oddbouind].value = *phiptr;
	++oddbouind;
      }
      else { //even case
	float *phiptr = phieven + h/2;
	phiboueven[evenbouind].phiptr = phiptr;
	phiboueven[evenbouind].value = *phiptr;
	++evenbouind;
	phiptr += grid_dim_minus_1_over_2;
	phiboueven[evenbouind].phiptr = phiptr;
	phiboueven[evenbouind].value = *phiptr;
	++evenbouind;
      }
    }
  }
}

// Translates a big array of epsilon values (eps_array) into a
// (hopefully) much smaller array of NonUnif, directly suitable for
// use in a finite difference algorithm.  The epsilon array is
// presumed offset from the phi array by (-1/2, -1/2, -1/2) so that
// the eight epsilon points (i,j,k), (i,j,k+1), ... (i+1,j+1,k+1) form
// a cube with phi point (i,j,k) at the center.  NonUnif "fac" fields
// are multiplicative factors that include the dielectric constant
// "between" two phi points, which is obtained by averaging over cube
// faces.  Here the "inverse averaging" scheme is used.

void FDGridLevel::set_up_sor ()
{
  blab2 << "Entering FDGridLevel::set_up_sor." << endl;

  if (! coarser) // This is the coarsest grid, so warn about outside charges.
    cubecharge = rho->get_cuberep(lat, true);
  else
    cubecharge = rho->get_cuberep(lat, false);

  fdchargeitr = new FDChargeIterator(cubecharge, phieven, phiodd);

  Coord half_offset (spacing/2.0F, spacing/2.0F, spacing/2.0F);
  CubeLatSpec offset_lattice(grid_dim, spacing, center - half_offset);
  dcrp = new DielCubeRep(offset_lattice);
  *dcrp = eps->get_cuberep(offset_lattice);
  const DielCubeRep &cubeeps = *dcrp;
  ElyCubeRep *cubeely = electrolyte->get_cuberep(lat);

  if (!cubeeps.eps_array_defined()) {
    ::error ("ERROR FDGridLevel::set_up_sor: eps_array undefined.  Can't continue.");
  }

  float spacingsq = spacing*spacing;

  vector<NonUnif> nuevenplex;
  vector<NonUnif> nuoddplex;

  inv_six_plus_kappa_even =
    (float *) big_rigid_malloc ((sizeof (float))*(ncube/2 + 1));
  inv_six_plus_kappa_odd =
    (float *) big_rigid_malloc ((sizeof (float))*(ncube/2));

  // It is necessary to put some sane value into these arrays in advance
  // because the big set-up loop below excludes all box boundary points
  // but the sor() main loops exclude only the +-i boundaries.  sor()
  // may get arithmatic exceptions if +-jk boundary points have "crazy"
  // values.  However, the final result for phi should be independent
  // of these boundary 1/(6+kappa) values, so 1/6 is OK even in salt.
  int hh;
  for (hh=0; hh<ncube/2; ++hh) {
    inv_six_plus_kappa_even[hh] = 1.0F/6.0F;
    inv_six_plus_kappa_odd[hh] = 1.0F/6.0F;
  }
  inv_six_plus_kappa_even[hh] = 1.0F/6.0F;


  float hueck = float(8 * pi * PhysCond::get_conconv()
    * electrolyte->ionic_strength()
    / PhysCond::get_kBolt() / PhysCond::get_T());

/*
Convention for designating corners and sides of epsilon box:
 Corners: (used in e, ek, etc.)
  i,   j,   k   :
  i,   j,   k+1 : k
  i,   j+1, k   : j
  i,   j+1, k+1 : jk
  i+1, j,   k   : i
  i+1, j,   k+1 : ik
  i+1, j+1, k   : ij
  i+1, j+1, k+1 : ijk
  Sides: (used in variable names, eps1, mask1, ...)
  i,j,k  to  i,   j,   k-1  :1
  i,j,k  to  i,   j,   k+1  :2
  i,j,k  to  i,   j-1, k    :3
  i,j,k  to  i,   j+1, k    :4
  i,j,k  to  i-1, j,   k,   :5
  i,j,k  to  i+1, j,   k,   :6

*/

  // This loop must exclude outer boundary points.
  int nm1 = grid_dim - 1;
  int i;
  for (i = 1; i<nm1; ++i) {
    int iterm = i*nsq;
    for (int j=1; j<nm1; ++j) {
      int jterm = j*grid_dim;
      for (int k=1; k<nm1; ++k) {
	int h = iterm + jterm + k;
	float *inv_six_plus_kappa = h%2 ? inv_six_plus_kappa_odd :
	   inv_six_plus_kappa_even;
	float epsave;
	// Is this a non uniform point?
	float e = cubeeps[h];
	if (e != cubeeps[h+1]
	    || e != cubeeps[h+grid_dim] || e != cubeeps[h+grid_dim+1]
	    || e != cubeeps[h+nsq] || e != cubeeps[h+nsq+1]
	    || e != cubeeps[h+nsq+grid_dim]
	    || e != cubeeps[h+nsq+grid_dim+1]) {
	  float ek = cubeeps[h+1];
	  float ej = cubeeps[h+grid_dim];
	  float ejk = cubeeps[h+grid_dim+1];
	  float ei = cubeeps[h+nsq];
	  float eik = cubeeps[h+nsq+1];
	  float eij = cubeeps[h+nsq+grid_dim];
	  float eijk = cubeeps[h+nsq+grid_dim+1];

	  float ep1;
	  float ep2;
	  float ep3;
	  float ep4;
	  float ep5;
	  float ep6;
	  if (epsave_oldway) {
	    ep1 = (e + ej + ei + eij) / 4.0F;
	    ep2 = (ek + ejk + eik + eijk) / 4.0F;
	    ep3 = (e + ek + ei + eik) / 4.0F;
	    ep4 = (ej + ejk + eij + eijk) / 4.0F;
	    ep5 = (e + ek + ej + ejk) / 4.0F;
	    ep6 = (ei + eik + eij + eijk) / 4.0F;
	  }
	  else {
	    ep1 = 4.0F / (1.0F/e + 1.0F/ej + 1.0F/ei + 1.0F/eij);
	    ep2 = 4.0F / (1.0F/ek + 1.0F/ejk + 1.0F/eik + 1.0F/eijk);
	    ep3 = 4.0F / (1.0F/e + 1.0F/ek + 1.0F/ei + 1.0F/eik);
	    ep4 = 4.0F / (1.0F/ej + 1.0F/ejk + 1.0F/eij + 1.0F/eijk);
	    ep5 = 4.0F / (1.0F/e + 1.0F/ek + 1.0F/ej + 1.0F/ejk);
	    ep6 = 4.0F / (1.0F/ei + 1.0F/eik + 1.0F/eij + 1.0F/eijk);
	  }
	  epsave = (ep1 + ep2 + ep3 + ep4 + ep5 + ep6) / 6.0F;

	  int halfh = h/2;
	  if (cubeely->is_guoy(h)) {
	    inv_six_plus_kappa[halfh]
	      = 1.0F / (6.0F + spacingsq * hueck / epsave);
	  }
	  else {
	    inv_six_plus_kappa[halfh] = 1.0F / 6.0F;
	  }

	  // Need some fast way of adding new elements to nu lists...
	  NonUnif newnu;
	  if (h%2) { // odd case
	    newnu.phiptr = phiodd + halfh;
	    float kappafac = inv_six_plus_kappa_odd[halfh];
	    newnu.ptr1 = phieven + halfh;  // same as + okm_off
	    newnu.fac1 = kappafac * (ep1/epsave - 1.0F);
	    newnu.ptr2 = phieven + halfh + 1;  // same as + okp_off
	    newnu.fac2 = kappafac * (ep2/epsave - 1.0F);
	    newnu.ptr3 = phieven + halfh + ojm_off;
	    newnu.fac3 = kappafac * (ep3/epsave - 1.0F);
	    newnu.ptr4 = phieven + halfh + ojp_off;
	    newnu.fac4 = kappafac * (ep4/epsave - 1.0F);
	    newnu.ptr5 = phieven + halfh + oim_off;
	    newnu.fac5 = kappafac * (ep5/epsave - 1.0F);
	    newnu.ptr6 = phieven + halfh + oip_off;
	    newnu.fac6 = kappafac * (ep6/epsave - 1.0F);
	    nuoddplex.push_back(newnu);
	  }
	  else {  // even case
	    newnu.phiptr = phieven + halfh;
	    float kappafac = inv_six_plus_kappa_even[halfh];
	    newnu.ptr1 = phiodd + halfh - 1;  // same as + ekm_off
	    newnu.fac1 = kappafac * (ep1/epsave - 1.0F);
	    newnu.ptr2 = phiodd + halfh;  // same as + ekp_off
	    newnu.fac2 = kappafac * (ep2/epsave - 1.0F);
	    newnu.ptr3 = phiodd + halfh + ejm_off;
	    newnu.fac3 = kappafac * (ep3/epsave - 1.0F);
	    newnu.ptr4 = phiodd + halfh + ejp_off;
	    newnu.fac4 = kappafac * (ep4/epsave - 1.0F);
	    newnu.ptr5 = phiodd + halfh + eim_off;
	    newnu.fac5 = kappafac * (ep5/epsave - 1.0F);
	    newnu.ptr6 = phiodd + halfh + eip_off;
	    newnu.fac6 = kappafac * (ep6/epsave - 1.0F);
	    nuevenplex.push_back(newnu);
	  }
	}
	else {
	  epsave = e;
	  if (cubeely->is_guoy(h)) {
	    inv_six_plus_kappa[h/2] = 1.0F / (6.0F + spacingsq * hueck / epsave);
	  }
	  else {
	    inv_six_plus_kappa[h/2] = 1.0F / 6.0F;
	  }
	}

	// Complete the charge data;
	fdchargeitr->complete_point(h, epsave, inv_six_plus_kappa[h/2],
				    phieven, phiodd);
      }
    }
  }

  // Set up the NonUnif arrays for use by sor().  (Plexs get deleted.)
  numnueven = nuevenplex.size();
  nueven = new NonUnif [numnueven];
  std::vector<NonUnif>::iterator p;
  for (i = 0, p = nuevenplex.begin();
       p != nuevenplex.end();
       i++, p++)
    nueven[i] = *p;
  numnuodd = nuoddplex.size();
  nuodd = new NonUnif [numnuodd];
  for (i = 0, p = nuoddplex.begin();
       p != nuoddplex.end();
       i++, p++)
    nuodd[i] = *p;
  delete cubeely;
  delete dcrp;
  dcrp = 0;
  blab2 << "Exiting FDGridLevel::set_up_sor." << endl;
}

void FDGridLevel::grid_to_c_array (float* data_array)
{
  if (!(phieven && phiodd))
    ::error ("INTERNAL ERROR: FDGridLevel::gred_to_c_array: \n",
	     "arrays phieven and/or phiodd not initialized\n");
  int h=0;
  for (int i=0; i<grid_dim; ++i) {
    for (int j=0; j<grid_dim; ++j) {
      for (int k=0; k<grid_dim; ++k) {
	int cind = i*nsq + j*grid_dim + k;
	if (cind!=h)
	  cerr << "INDICES OUT OF WHACK!" << endl;
	if (h%2)
	  data_array[cind] = phiodd[h/2];
	else
	  data_array[cind] = phieven[h/2];
	++h;
      }
    }
  }
}


// More specialization of AvsScalarField:

float FDGridLevel::avsScaleFactor = 1.0;

bool
FDGridLevel::put_val_array (const float* data_array)
     // Arg is a FORTRAN-ordered array!
{
  if (!(phieven && phiodd))
    ::error ("INTERNAL ERROR: FDGridLevel::put_val_array: \n",
	     "arrays phieven and/or phiodd not initialized\n");
  int h=0;
  for (int i=0; i<grid_dim; ++i) {
    for (int j=0; j<grid_dim; ++j) {
      for (int k=0; k<grid_dim; ++k) {
	int fortind = k*nsq + j*grid_dim + i;
	if (h%2)
	  phiodd[h/2] = data_array[fortind]/avsScaleFactor;
	else
	  phieven[h/2] = data_array[fortind]/avsScaleFactor;
	++h;
      }
    }
  }
  return true;
}


bool
FDGridLevel::get_val_array (float* data_array) const
     // Arg is a FORTRAN-ordered array!
{
  if (!(phieven && phiodd))
    ::error ("INTERNAL ERROR: FDGridLevel::get_val_array: \n",
	     "arrays phieven and/or phiodd not initialized\n");
  int h=0;
  for (int i=0; i<grid_dim; ++i) {
    for (int j=0; j<grid_dim; ++j) {
      for (int k=0; k<grid_dim; ++k) {
	int fortind = k*nsq + j*grid_dim + i;
	if (h%2)
	  data_array[fortind] = phiodd[h/2];
	else
	  data_array[fortind] = phieven[h/2];
	++h;
      }
    }
  }
  if (avsScaleFactor != 1.0) {
    const float * enddat = data_array + ncube;
    for (float *p = data_array; p < enddat; ++p)
      *p *= avsScaleFactor;
  }
  return true;
}

// FDGridLevel.cc ends here
