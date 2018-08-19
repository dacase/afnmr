/* Shell is sphere with inner and outer radius and set of Cones  -*- C++ -*-
    Copyright (c) 1993--1995 by Donald Bashford and Tony You.

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

$Id: Shell.h,v 2.6 2004/12/06 17:58:38 bashford Exp $
*/
#ifndef _Shell_h
#define _Shell_h 1

#include "MEAD/globals.h"
#include "MEAD/Atom.h"
#include "MEAD/Coord.h"
#include "MEAD/Cone.h"
#include <iostream>

class Shell {
friend class ShellConeIterator;
public:

  Shell (const Coord& c, float inner_r, float probe_r );
  Shell (const Atom& a, float probe_r );

  Shell ();
  Shell (const Shell& );
  Shell& operator = (const Shell& );

  inline Coord get_coord() const {return coord;}
  inline float get_inner_rad() const {return inner_rad;}
  inline float get_outer_rad() const {return outer_rad;}
  inline float get_outer_rad_sq() const {return outer_rad_sq;}
  inline int is_buried() const {return (flag == buried);}
  inline int is_partially_buried() const {return (flag==partially_buried);}
  inline int is_free() const {return (flag==free);}
  inline int get_numcones() const {return numcones;}
  inline void bury() {flag=buried;}
  inline void partially_bury() {flag=partially_buried;}
  void add_cone(const Cone&);
  ostream& write_top_in_binary(ostream&);
  ostream& write_top_in_ascii(ostream&);
  int read_top_in_binary(istream&);
  int read_top_in_ascii(istream&);
  ~Shell( );

protected:
  void delete_conelist();

  Coord coord;
  float inner_rad, outer_rad, outer_rad_sq;
  enum Flag {free=0, partially_buried=1, buried=2 };
  Flag flag;
  int numcones;
  Cone *cone_head;
};

class ShellConeIterator {
public:
  inline ShellConeIterator(const Shell& s) {cur = s.cone_head;}
  int next();
  inline const Cone* cur_cone_pt() const {return (const Cone *) cur;}
  inline int valid() const {return cur ? 1 : 0;}
private:
  Cone *cur;
  // Bogus functions meant to prevent empty construction and copying
  ShellConeIterator(){}
  ShellConeIterator(const ShellConeIterator&){}
  ShellConeIterator& operator=(const ShellConeIterator&){return *this;}
};



inline
Shell::Shell ()
:coord(0.0, 0.0, 0.0), inner_rad(0), outer_rad(0), outer_rad_sq(0),
 numcones(0), flag(free), cone_head(0)
{}

inline
Shell::~Shell( )
{
  if (numcones)
    delete_conelist();
}

inline int ShellConeIterator::next()
{
  if (!cur)
    ::error("ERROR ShellConeIterator ran off end");
  cur = cur->next;
  int retval = cur ? 1 : 0;
  return retval;
}


#endif

// Shell.h ends here
