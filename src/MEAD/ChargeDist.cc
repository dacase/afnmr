/* Base class for charge distributions, envelope for subclasses.

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

$Id: ChargeDist.cc,v 2.13 2007/05/28 01:26:42 bashford Exp $
*/
#include "MEAD/ChargeDist.h"
#include "MEAD/globals.h"


// FIXME this include can be left out as long as operator* is left out
//#include "ElstatPot.h"


int ChargeDist::instanceCount = 0;

#if defined(__hpux) && ! defined(__GNUG__)

// The HP C++ compiler that came with HPUX 10.20 wants template
// definitions to be in a file with the same basename as the .h
// file in which the template is declared.  (Sheesh!)

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

class NoCharge : public ChargeDist_lett {
public:
  NoCharge() : ChargeDist_lett() {}
  NoCharge(const NoCharge& a)
    : ChargeDist_lett(a) {}
  NoCharge& operator=(const NoCharge& acs)
    { ChargeDist_lett::operator=(acs);
      return *this;}
  virtual ChargeCubeRep* get_cuberep(const CubeLatSpec&,
				     bool warn_outside=false) { return 0; }
  virtual float total_charge() const  { return 0.0; }
  virtual int has_charges() const  { return 0; }
  virtual float vacuum_coulomb(const Coord& c) const { return 0; }
  virtual size_t number_points() const { return 0; }
  class const_iterator {
  public:
    const_iterator() {}
    virtual const PointCharge& operator*()
    { error("NoCharge::const_iterotor::operator* should never be called\n",
	    "since a NoCharge has no points");
      return *new PointCharge();
    }
    virtual const_iterator& operator++() {return *this;}
    virtual bool operator!=(const const_iterator& x) const
      {return false;}
  private:
    PointCharge _pc;
  };
  ChargeDist::const_iterator pc_begin() const
    {return make_PointCharge_const_iterator(the_itr);}
  ChargeDist::const_iterator pc_end() const
    {return make_PointCharge_const_iterator(the_itr);}

  virtual bool dispatch_addChargeDist (const ChargeDist_lett& c,
				       ChargeDist &result) const ;


private:
  const_iterator the_itr;
};


bool NoCharge::dispatch_addChargeDist (const ChargeDist_lett& c,
				       ChargeDist &result) const
{
  result = ChargeDist((ChargeDist_lett*) &c);  // Casting away const!
  return true;
}

// Constructs an envelope with a ChargDist letter, that is, a blank letter.
ChargeDist::ChargeDist()
{

  rep = new NoCharge();
}

ChargeDist::ChargeDist(const ChargeDist & c)
{
  if(c.rep) { // Arg is an envelope.
    c.rep->referenceCount++;
    rep = c.rep;
  }
  else
    error("ERROR: ChargeDist copy constructor:\n",
	  "Arugments, `rep' pointer is zero\n");
  ++instanceCount;
}

ChargeDist& ChargeDist::operator = (const ChargeDist& c)
{
  if (this != &c) {
    if (rep) {
      if (c.rep) {
        c.rep->referenceCount++;
        if (rep->referenceCount == 1)
          delete rep;
        else
          --rep->referenceCount;
        rep = c.rep;
      }
      else {
        error("ERROR: ElstatPot::operator=:\n", "argument has zero rep");
      }
    }
    else {
      error("ERROR ElstatPot::operator=: this->rep is zero");
    }
  }
  return *this;
}

// This will be an envelope containing *cdp.
ChargeDist::ChargeDist(ChargeDist_lett * cdp)
{
  cdp->referenceCount++;
  rep = cdp;
  instanceCount++;
}

ChargeDist::~ChargeDist()
{
  blab2 << "ChargeDist destructor called with rep->referenceCount = " << rep->referenceCount << endl;
  if (rep) {
    if (rep->referenceCount == 1)
      delete rep;
    else
      --rep->referenceCount;
  }
  instanceCount--;
}


// FIXME! Commented out for now to prevent too much file interdependency
// at compile time.
//float operator* (const ChargeDist& c, const ElstatPot& e)
//{return e.operator* (c);}

ChargeDist ChargeDist::addChargeDist (const ChargeDist& c) const
{
  blab2 << "ChargeDist::addChargeDist (const ChargeDist&) called" << endl;
  ChargeDist result;
  if (! dispatch_addChargeDist(c, result) ) {
    blab2 << "ChargeDist::dispatch_addChargeDist failed the first way, try other" << endl;
    if (! c.dispatch_addChargeDist(*this, result) ) {
      cerr << "ERROR: ChargeDist::dispatch_addChargeDist failed both ways"
	<< endl;
    }
  }
  return result;
}

int ChargeDist_lett::instanceCount = 0;


ChargeDist_lett::ChargeDist_lett()
{
  referenceCount=1;
  ++instanceCount;
}


ChargeDist_lett::ChargeDist_lett(const ChargeDist_lett&)
{
  referenceCount = 1;
  ++instanceCount;
}

// Don't delete this letter until all references have been released
void ChargeDist_lett::operator delete(void *p, size_t s)
{
  ChargeDist_lett *del = reinterpret_cast<ChargeDist_lett*> (p);
  blab2 << "ChargeDist_lett operator delete called with referenceCount = " << del->referenceCount << endl;
  if (del->referenceCount == 0) {
    blab2 << "ChargeDist_lett operator delete frees object" << endl;
    ::operator delete(p);
    }
}

ChargeDist_lett::~ChargeDist_lett()
{
  blab2 << "ChargeDist_lett destructor called with referenceCount = " << referenceCount << endl;
  if (--referenceCount)
    // It's now OK for this to happen
    blab2 << "WARNING: ChargeDist_lett destructor called\n"
      << "on object with referenceCount = " << referenceCount+1
	<< "\nas if the object is still in use somewhere" << endl;
  --instanceCount;
}

ChargeDist_lett& ChargeDist_lett::operator=(const ChargeDist_lett&)
{return *this;}

// ChargeDist.cc ends here
