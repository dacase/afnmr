#include "MEAD/ChargeDist.h"
#include "MEAD/ElstatPot.h"

// These work for all letter types, but more specialized versions
// would be faster.  FIXME.

typedef ChargeDist::const_iterator cdcitr;

float operator* (const ChargeDist& c, const ElstatPot& e)
{
  float product = 0;
  for (cdcitr i = c.begin(); i!=c.end(); ++i) {
    PointCharge pc = *i;
    product += pc.charge * e.value(pc.coord);
  }
  return product;
}

float operator* (const ChargeDist& c, const ElstatPot_lett& e)
{
  float product = 0;
  for (cdcitr i = c.begin(); i!=c.end(); ++i) {
    PointCharge pc = *i;
    product += pc.charge * e.value(pc.coord);
  }
  return product;
}

float operator* (const ChargeDist_lett& c, const ElstatPot& e)
{
  float product = 0;
  for (cdcitr i = c.pc_begin(); i!=c.pc_end(); ++i) {
    PointCharge pc = *i;
    product += pc.charge * e.value(pc.coord);
  }
  return product;
}

float operator* (const ChargeDist_lett& c, const ElstatPot_lett& e)
{
  float product = 0;
  for (cdcitr i = c.pc_begin(); i!=c.pc_end(); ++i) {
    PointCharge pc = *i;
    product += pc.charge * e.value(pc.coord);
  }
  return product;
}

// multPotxCharge.cc ends here
