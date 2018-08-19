// -*- C++ -*-
#ifndef _PointCharge_h
#define _PointCharge_h

/*
  A PointCharge has a position (coord) and a charge (charge).
  PointCharge_const_iterator is a polymorphic iterator that can
  be created for any kind of container containing objects capable
  of being converted into PointCharges by a function
  convert_to_PointCharge(const whatever&), where "whatever" might be,
  for example, Atom.

  For details see the notes and code examples in the directory,
  t-stuff.
*/

#include "MEAD/Coord.h"

//!wrap!
struct PointCharge {
  PointCharge();
  PointCharge(float charge, const Coord& coord);
  PointCharge(const Coord& coord, float charge);
  ~PointCharge() {}
  bool operator==(const PointCharge& pc) const;
  float charge;
  Coord coord;
};

inline PointCharge::PointCharge()
    : charge(0.0), coord(Coord()) {}
inline PointCharge::PointCharge(float charge, const Coord& coord)
    : charge(charge), coord(coord) {}
inline PointCharge::PointCharge(const Coord& coord, float charge)
    : charge(charge), coord(coord) {}

inline bool PointCharge::operator==(const PointCharge& pc) const
{ return (coord == pc.coord && charge == pc.charge); }

/*
  This is an "envelope" base class that manages the memory of,
  and forwards operations to a "letter" which is of a derived class.
  Polymorphism comes through the virtual function mechanism
  with respect to the pointer to the letter.
  The templated ctor, PointCharge_const_iterator(const FromIter&),
  provides a way of encapsulating other iterators.
  Example usage:

  f(const vector<Atom>& atch)
  {
    for (PointCharge_const_iterator j = atch.begin(); j != atch.end(); ++j) {
      const PointCharge& pc = *j;
      // ... use pc ...
  }
*/

class PointCharge_const_iterator {
public:
  // Standard canonical ctors and dtor and assign
  PointCharge_const_iterator() :itrptr(0), refcount(0) {}
  PointCharge_const_iterator(const  PointCharge_const_iterator& i)
    : itrptr(i.itrptr), refcount(0)
    {if (itrptr) itrptr->refcount++; }
  PointCharge_const_iterator(PointCharge_const_iterator* p)
    : itrptr(p), refcount(0)
      { itrptr->refcount = 1; }
  virtual ~PointCharge_const_iterator()
    {
      if (itrptr) {
	itrptr->refcount--;
	if (itrptr->refcount == 0) delete itrptr;
      }
    }
  PointCharge_const_iterator& operator=(const PointCharge_const_iterator& i)
    {
      if (itrptr) {
	itrptr->refcount--;
	if (itrptr->refcount == 0) delete itrptr;
      }
      itrptr = i.itrptr;
      if (itrptr) itrptr->refcount++;
      return *this;
    }
  // ctor encapsulating a concrete iterator polymorphically in this object.
  // iterator operations, forwarded to itrptr
  virtual const PointCharge& operator*()
    {return **itrptr;}
  virtual PointCharge_const_iterator& operator++()
    {++(*itrptr); return *this;}
  virtual bool operator!=(const PointCharge_const_iterator& x) const
    {return (*itrptr) != (*x.itrptr);}
private:
  PointCharge_const_iterator *itrptr;
  short unsigned refcount;
};

#include <assert.h>

template<class FromIter>
class L_Const_Iter : public PointCharge_const_iterator {
public:
  L_Const_Iter(const FromIter& i)
    :  _itr(i) {}
  virtual const PointCharge& operator*()
    {val = convert_to_PointCharge(*_itr); return val;}
  virtual PointCharge_const_iterator& operator++()
    {++_itr; return *this;}
  virtual bool operator!=(const PointCharge_const_iterator& x) const
    { // The real type of the arg had better be the same is this!
      const L_Const_Iter* xc = dynamic_cast<const L_Const_Iter*>(&x);
      assert(xc!=0);
      return _itr != xc->_itr;
    }
private:
  FromIter _itr;
  PointCharge val;
};

template<class FromIter>
PointCharge_const_iterator
make_PointCharge_const_iterator(const FromIter i)
{
  return PointCharge_const_iterator(new L_Const_Iter<FromIter> (i));
}

// Duh!  In case the source iterator already points to PointCharges:
inline PointCharge convert_to_PointCharge (const PointCharge& pc)
{return pc;}


#endif

// PointCharge.h ends here
