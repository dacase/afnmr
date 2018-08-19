/* A level of grid spacing/size in finite difference procedure. -*- C++ -*-

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

$Id: FDGridLevel.h,v 2.19 2011/02/10 23:24:30 agiammon Exp $
*/
#ifndef _FDGridLevel_h
#define _FDGridLevel_h 1

#include <string>
using std::string;
#include "MEAD/Coord.h"
#include "MEAD/CubeLatSpec.h"
#include "MEAD/AvsScalarField.h"

class AnalyticEP;

#include "MEAD/NonUnif.h"

class ChargeDist_lett;
class DielectricEnvironment_lett;
class ElectrolyteEnvironment_lett;
class ChargeCubeRep;
class FDChargeIterator;
class DielCubeRep;

class FDGridLevel : public AvsScalarField {
friend class FinDiffElstatPot;
public:
  // There is no coarser grid
  FDGridLevel(const CubeLatSpec&, ChargeDist_lett*,
	      DielectricEnvironment_lett*,
	      ElectrolyteEnvironment_lett*, AnalyticEP*);
  // Will initialize from coarser grid
  FDGridLevel(const CubeLatSpec&, ChargeDist_lett*,
	      DielectricEnvironment_lett*,
	      ElectrolyteEnvironment_lett*, FDGridLevel * coarser_fdl);
  // Ctors for installing a pre-determined grid that is presumed to be
  // the solution.  Giving coarser grid allows focusing and giving
  // Diel allows electric displacement vectors.
  // solve() should not be called; potint, etc will use values given here.
  FDGridLevel (const CubeLatSpec&, float *phi);
  FDGridLevel (const CubeLatSpec&, float *phi, DielectricEnvironment_lett*);
  FDGridLevel (const CubeLatSpec&, float *phi, FDGridLevel* coarser);
  FDGridLevel (const CubeLatSpec&, float *phi, DielectricEnvironment_lett*,
	       FDGridLevel* coarser);
  ~FDGridLevel();
  virtual void phi_init();
  virtual void phi_from_prev();
  virtual void solve();

  // Provide the interface to AvsScalarField including required value func.
  virtual void solve_using_initial(const string& initial_fieldname);
  virtual float value(Coord c) const {return potint(c);}
  virtual bool putval (Coord, float) {return false;} // Basically, a no-op.
  virtual void write(const CubeLatSpec& cls, const string& s) const
    {write_uniform_native(cls, s);}
  virtual void write(const string& s) const
    {write_uniform_native(lat, s);}
  // Instead of value and putval, *_val_array functions
  // to overload those of AvsScalarField
  virtual void grid_to_c_array (float* data_array);
private:
  virtual bool get_val_array (float*) const ; // ditto
  virtual bool put_val_array (const float*); // Arg is a FORTRAN-ordered array!
public:

  virtual float phiof_ijk(int i, int j, int k) const;
  virtual float potint (float xp, float yp, float zp) const;
  virtual float potint (const Coord& c) const {return potint(c.x, c.y, c.z);}
  virtual Coord fieldint (const Coord& c) const;
  virtual Coord displacement_int (const Coord& p); // const in principle but might rebuild dcrp in practice FIXME?
  static bool get_epsave_oldway() {return epsave_oldway;}
  static void set_epsave_oldway(bool b) {epsave_oldway = b;}
  static bool get_use_fixed_maxrmsdiff() {return use_fixed_maxrmsdiff;}
  static void set_use_fixed_maxrmsdiff(bool b) {use_fixed_maxrmsdiff = b;}
  static float get_avsScaleFactor() {return avsScaleFactor;}
  static void set_avsScaleFactor(float f) {avsScaleFactor = f;}
  static bool get_ZeroForOutOfRange() {return ZeroForOutOfRange;};
  static void set_ZeroForOutOfRange(bool b) {ZeroForOutOfRange = b;};
protected:
  static bool ZeroForOutOfRange;
  static bool epsave_oldway;
  static bool use_fixed_maxrmsdiff;
  static float avsScaleFactor; // grid is multiplied by this before output
                               // and after input.
  virtual void misc_settings();
  virtual void init_sides_from_phi();
  virtual void set_up_sor();
  virtual void sor();

  FDGridLevel *coarser;
  CubeLatSpec lat;
  Coord center;
  float spacing;
  int grid_dim;
  DielCubeRep *dcrp;

  ChargeDist_lett* rho;
  DielectricEnvironment_lett* eps;
  ElectrolyteEnvironment_lett* electrolyte;
  AnalyticEP *analytic_approx;

  int solved;
  float *phieven, *phiodd;
  float *inv_six_plus_kappa_even;
  float *inv_six_plus_kappa_odd;
  size_t numnueven, numnuodd;
  NonUnif *nueven, *nuodd;

  ChargeCubeRep *cubecharge;
  FDChargeIterator *fdchargeitr;

 // for restoring correct phi's at outer boundary.
  int num_outbou_even, num_outbou_odd;
  struct OuterBoundPoint {
    float *phiptr;
    float value;
  } *phiboueven, *phibouodd;

  // Some variables for speed and convenience
  int nsq;
  int ncube;
  float xhi, xlo, yhi, ylo, zhi, zlo;
  // phieven[i]'s j-1 neighbor is phiodd[i+ejm_off] etc.
  int ekm_off, ekp_off, ejm_off, ejp_off, eim_off, eip_off;
  // phiodd[i]'s j-1 neighbor is phieven[i+ojm_off] etc.
  int okm_off, okp_off, ojm_off, ojp_off, oim_off, oip_off;

};

#include "MEAD/globals.h"

inline float FDGridLevel::phiof_ijk(int i, int j, int k) const
{
  if (!solved)
    ::error("ERROR: FDGridLevel::phiof_ijk called for unsolved object\n");
  int h = i*nsq + j*grid_dim + k;
  return h%2 ? phiodd[h/2]  : phieven[h/2];
}


#endif

// FDGridLevel.h ends here
