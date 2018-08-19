/* Class for electrostatic potential combinations -*- C++ -*-

    Copyright (c) 1993--2001 by Donald Bashford

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

$Id: ElstatPotCombination.cc,v 2.3 2007/05/28 01:26:42 bashford Exp $
*/

#include "MEAD/ElstatPotCombination.h"

ElstatPotCombination::ElstatPotCombination() : global_scale(1.0) { }

ElstatPotCombination::ElstatPotCombination(const ElstatPotCombination& epc)
{
  global_scale = epc.get_scale();
  for (ElstatPotCombination::const_iterator ind = epc.begin(); ind != epc.end(); ++ind) {
    push_back(*ind);
  }
}

ElstatPotCombination::ElstatPotCombination(ElstatPot_lett& epl, float scale) : global_scale(scale)
{
  push_back(pair<ElstatPot, float>(ElstatPot(&epl), 1.0));
}

ElstatPotCombination::ElstatPotCombination(ElstatPot_lett& epl1, ElstatPot_lett& epl2)
{
  global_scale = 1.0;
  push_back(pair<ElstatPot, float>(ElstatPot(&epl1), 1.0));
  push_back(pair<ElstatPot, float>(ElstatPot(&epl2), 1.0));
}

ElstatPotCombination::ElstatPotCombination(ElstatPot& ep, float scale) : global_scale(scale)
{
  push_back(pair<ElstatPot, float>(ep, 1.0));
}

ElstatPotCombination::ElstatPotCombination(ElstatPot& ep1, ElstatPot& ep2)
{
  global_scale = 1.0;
  push_back(pair<ElstatPot, float>(ep1, 1.0));
  push_back(pair<ElstatPot, float>(ep2, 1.0));
}

ElstatPotCombination::ElstatPotCombination(ElstatPot& ep, ElstatPot_lett& epl)
{
  global_scale = 1.0;
  push_back(pair<ElstatPot, float>(ep, 1.0));
  push_back(pair<ElstatPot, float>(ElstatPot(&epl), 1.0));
}

ElstatPotCombination::ElstatPotCombination(ElstatPot_lett& epl, ElstatPot& ep)
{
  global_scale = 1.0;
  push_back(pair<ElstatPot, float>(ep, 1.0));
  push_back(pair<ElstatPot, float>(ElstatPot(&epl), 1.0));
}

ElstatPotCombination& ElstatPotCombination::operator = (const ElstatPotCombination& epc) {
  list<pair<ElstatPot, float> >::operator=(epc);
  global_scale = epc.get_scale();
  return *this;
}

void ElstatPotCombination::solve() {
  for (ElstatPotCombination::const_iterator ind = begin(); ind != end(); ++ind) {
    pair<ElstatPot, float> p = *ind;
    p.first.solve();
  }
}

float ElstatPotCombination::value(Coord x) const {
  float v = 0.0;
  for (ElstatPotCombination::const_iterator ind = begin(); ind != end(); ++ind) {
    pair<ElstatPot, float> p = *ind;
    v += p.first.value(x) * p.second;
  }
  return v * global_scale;
}

Coord ElstatPotCombination::field(Coord x) const {
  Coord f(0.0, 0.0, 0.0);
  for (ElstatPotCombination::const_iterator ind = begin(); ind != end(); ++ind) {
    pair<ElstatPot, float> p = *ind;
    f += p.first.field(x) * p.second;
  }
  return f * global_scale;
}

Coord ElstatPotCombination::displacement(Coord x) const {
  Coord d(0.0, 0.0, 0.0);
  for (ElstatPotCombination::const_iterator ind = begin(); ind != end(); ++ind) {
    pair<ElstatPot, float> p = *ind;
    d += p.first.displacement(x) * p.second;
  }
  return d * global_scale;
}

ElstatPotCombination& ElstatPotCombination::operator += (const ElstatPotCombination& epc) {
  if (global_scale != epc.global_scale) {
    error("Can't add ElstatPotCombinations with different scale factors");
  }
  for (ElstatPotCombination::const_iterator ind = epc.begin(); ind != epc.end(); ++ind) {
    push_back(*ind);
  }
  return *this;
}

ElstatPotCombination& ElstatPotCombination::operator += (ElstatPot_lett& epl) {
  push_back(pair<ElstatPot, float>(ElstatPot(&epl), 1.0));
  return *this;
}

ElstatPotCombination& ElstatPotCombination::operator += (ElstatPot& ep) {
  push_back(pair<ElstatPot, float>(ep, 1.0));
  return *this;
}

ElstatPotCombination& ElstatPotCombination::operator *= (float scale) {
  global_scale *= scale;
  return *this;
}

ElstatPotCombination& ElstatPotCombination::operator /= (float scale) {
  if (scale == 0.0) {
    error("Attempt to divide ElstatPotCombination by zero");
  }
  global_scale /= scale;
  return *this;
}

ElstatPotCombination operator + (const ElstatPotCombination& epc1, const ElstatPotCombination& epc2) {
  if (epc1.get_scale() != epc2.get_scale()) {
    error("Can't add ElstatPotCombinations with different scale factors");
  }
  ElstatPotCombination sum(epc1);
  return sum += epc2;
}

ElstatPotCombination operator + (const ElstatPotCombination& epc, ElstatPot_lett& epl) {
  ElstatPotCombination sum(epc);
  return sum += epl;
}

ElstatPotCombination operator + (ElstatPot_lett& epl, const ElstatPotCombination& epc) {
  return epc + epl;
}

ElstatPotCombination operator + (const ElstatPotCombination& epc, ElstatPot& ep) {
  ElstatPotCombination sum(epc);
  return sum += ep;
}

ElstatPotCombination operator + (ElstatPot& ep, const ElstatPotCombination& epc) {
  return epc + ep;
}

ElstatPotCombination operator * (const ElstatPotCombination& epc, float scale) {
  ElstatPotCombination prod(epc);
  return prod *= scale;
}

ElstatPotCombination operator * (float scale, const ElstatPotCombination& epc) {
  return epc * scale;
}

ElstatPotCombination operator / (const ElstatPotCombination& epc, float scale) {
  if (scale == 0.0) {
    error("Attempt to divide ElstatPotCombination by zero");
  }
  ElstatPotCombination div(epc);
  return div /= scale;
}

ElstatPotCombination operator / (float scale, const ElstatPotCombination& epc) {
  return epc / scale;
}

float operator * (const ChargeDist& cd, const ElstatPotCombination& epc) {
  float energy = 0.0;
  for (ElstatPotCombination::const_iterator ind = epc.begin(); ind != epc.end(); ++ind) {
    pair<ElstatPot, float> p = *ind;
    energy += cd * p.first;
    energy *= p.second;
  }
  return energy * epc.get_scale();
}

float operator * (const ElstatPotCombination& epc, const ChargeDist& cd) {
  return cd * epc;
}

float operator * (const ChargeDist_lett& cd, const ElstatPotCombination& epc) {
  float energy = 0.0;
  for (ElstatPotCombination::const_iterator ind = epc.begin(); ind != epc.end(); ++ind) {
    pair<ElstatPot, float> p = *ind;
    energy += cd * p.first;
    energy *= p.second;
  }
  return energy * epc.get_scale();
}

float operator * (const ElstatPotCombination& epc, const ChargeDist_lett& cd) {
  return cd * epc;
}

ElstatPotCombination operator* (ElstatPot& ep, float scale)
{ ElstatPotCombination epc(ep, scale); return epc; }

ElstatPotCombination operator* (float scale, ElstatPot& ep)
{ ElstatPotCombination epc(ep, scale); return epc; }

ElstatPotCombination operator* (ElstatPot_lett& epl, float scale)
{ ElstatPotCombination epc(epl, scale); return epc; }

ElstatPotCombination operator* (float scale, ElstatPot_lett& epl)
{ ElstatPotCombination epc(epl, scale); return epc; }

ElstatPotCombination operator/ (ElstatPot& ep, float scale)
{ ElstatPotCombination epc(ep, 1.0F/scale); return epc; }

ElstatPotCombination operator/ (ElstatPot_lett& epl, float scale)
{ ElstatPotCombination epc(epl, 1.0F/scale); return epc; }

ElstatPotCombination operator+ (ElstatPot_lett& epl1, ElstatPot_lett& epl2)
{ ElstatPotCombination epc(epl1, epl2); return epc; }

ElstatPotCombination operator+ (ElstatPot_lett& epl, ElstatPot& ep)
{ ElstatPotCombination epc(epl, ep); return epc; }

ElstatPotCombination operator+ (ElstatPot& ep, ElstatPot_lett& epl)
{ ElstatPotCombination epc(ep, epl); return epc; }

ElstatPotCombination operator+ (ElstatPot& ep1, ElstatPot& ep2)
{ ElstatPotCombination epc(ep1, ep2); return epc; }

