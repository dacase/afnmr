/* Base class for charge distributions and envelope for subclasses -*- C++ -*-

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

$Id: ChargeDist.h,v 2.21 2007/05/28 01:26:42 bashford Exp $
*/
#ifndef _ChargeDist_h
#define _ChargeDist_h 1
#include "MEAD/globals.h"
#include "MEAD/Coord.h"
#include "MEAD/PointCharge.h"


class CubeLatSpec;
class ChargeCubeRep;
class ElstatPot;
class Coord;
class ChargeDist_lett;

// FIXME? Maybe define ChargeDist::operator*(ChargeDist) to calculate the
// coulombic interation of two ChargeDists.

class CDL_const_iterator;

//DON'T !wrap!
class ChargeDist {
  //!nowrap!+
friend class ElstatPot;
  //!nowrap!-
public:
  //!nowrap!+
  ChargeDist();
  ChargeDist(const ChargeDist & c);
  ChargeDist& operator = (const ChargeDist&);
  //!nowrap!-
  ChargeDist(ChargeDist_lett *);
  //!nowrap!+
  virtual ChargeCubeRep* get_cuberep(const CubeLatSpec& cls,
				     bool warn_outside=false) ;
  virtual float total_charge() const ;
  virtual int has_charges() const ;
  virtual size_t number_points() const ;
  typedef PointCharge_const_iterator const_iterator;
  const_iterator begin() const;
  const_iterator end() const;

  virtual float vacuum_coulomb(const Coord& c) const ;

  virtual ChargeDist addChargeDist (const ChargeDist& c) const ;
  virtual bool dispatch_addChargeDist (const ChargeDist& c,
				       ChargeDist &result) const ;

  static int instances() {return instanceCount;}
  virtual ~ChargeDist();
  //!nowrap!-

private:
  ChargeDist_lett *rep;
  static int instanceCount;
};

/* These members (or something) should be implemented. FIXME!
class ChargeDist {
public:

  virtual ChargeDist operator + (const ChargeDist& c) const
    {return rep->operator+ (c);}
  virtual ChargeDist operator - (const ChargeDist& c) const
    {return rep->operator- (c);}

  virtual float debyeMultiply(const Debye& phi) const
    {return rep->debyeMultiply(phi);}
  virtual float analySphereMultiply(const AnalySphere& phi) const
    {return rep->analySphereMultiply(phi);}
  virtual float analySlabMultiply(const AnalySlab& phi) const
    {return rep->analySlabMultiply(phi);}
  virtual float finDiffElstatPotMultiply(const FinDiffElstatPot& phi) const
    {return rep->finDiffElstatPotMultiply(phi);}

};
*/

// FIXME! Commented out for now to prevent too much file interdependency
//float operator* (const ChargeDist& c, const ElstatPot& e);

// THE PURE VIRTUAL LETTER

//!wrap!
class ChargeDist_lett {
  //!nowrap!+
friend class ChargeDist;
  //!nowrap!-
public:
  //!nowrap!+
  ChargeDist_lett();
  ChargeDist_lett(const ChargeDist_lett& a);
  ChargeDist_lett& operator=(const ChargeDist_lett& acs);
  virtual ~ChargeDist_lett();

  virtual ChargeCubeRep* get_cuberep(const CubeLatSpec&,
				     bool warn_outside=false) = 0;
  virtual float total_charge() const = 0;
  virtual int has_charges() const = 0;
  virtual float vacuum_coulomb(const Coord& c) const = 0;
  virtual size_t number_points() const = 0;
  virtual ChargeDist::const_iterator pc_begin() const =0;
  virtual ChargeDist::const_iterator pc_end() const =0;

  virtual bool dispatch_addChargeDist (const ChargeDist_lett& c,
				       ChargeDist &result) const = 0;
  static int instances() {return instanceCount;}
  void operator delete (void *p, size_t s);
  //!nowrap!-

private:
  static int instanceCount;
  short referenceCount;
};

inline ChargeCubeRep* ChargeDist::get_cuberep(const CubeLatSpec& cls,
					      bool warn_outside)
{return rep->get_cuberep(cls, warn_outside);}

inline float ChargeDist::total_charge() const
{return rep->total_charge();}

inline int ChargeDist::has_charges() const
{return rep->has_charges();}

inline size_t ChargeDist::number_points() const
{return rep->number_points();}

inline ChargeDist::const_iterator ChargeDist::begin() const
{return rep->pc_begin();}

inline ChargeDist::const_iterator ChargeDist::end() const
{return rep->pc_end();}

inline float ChargeDist::vacuum_coulomb(const Coord& c) const
{return rep->vacuum_coulomb(c);}


inline bool ChargeDist::dispatch_addChargeDist (const ChargeDist& c,
						     ChargeDist &result) const
{ return rep->dispatch_addChargeDist(*c.rep, result); }


/*
template<class Container, class ValConverter>
class CDL_Const_Iter : public CDL_const_iterator,
		     public Container::const_iterator {
  typedef typename Container::const_iterator base_iter;
public:
  CDL_Const_Iter(const typename Container::const_iterator& i)
    : Container::const_iterator(i) {}
  virtual const PointCharge& operator*()
    {val = ValConverter(base_iter::operator*()); return val;}
  virtual CDL_const_iterator& operator++()
    {base_iter::operator++(); return *this;}
  virtual bool operator!=(const CDL_const_iterator& x) const
    {return base_iter::operator!=(dynamic_cast<const CDL_Const_Iter&>(x));}
private:
  PointCharge val;
};
*/

// The HP C++ compiler needs function template forward declaration
// here, and template definition elsewhere (try_addChargeDist.cc).
#if defined(__hpux) && ! defined(__GNUG__)
template<class C1, class C2>
bool try_addChargeDist (const C1& dummy, const C2& c2,
			const ChargeDist_lett& c, ChargeDist& result);
#else
template<class C1, class C2>
bool try_addChargeDist (const C1& dummy, const C2& c2,
			const ChargeDist_lett& c, ChargeDist& result)
{
  const C1 *c1ptr = dynamic_cast<const C1*>(&c);
  if (c1ptr) {
    cerr << "About to do result = c2.addChargeDist(*c1ptr)" << endl;
    result = c2.addChargeDist(*c1ptr);
    cerr << "Just finished result = c2.addChargeDist(*c1ptr)" << endl;
    return true;
  }
  else
    return false;
}
#endif

#endif

// ChargeDist.h ends here
