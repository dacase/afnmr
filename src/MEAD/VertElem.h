// This is -*- C++ -*-
#ifndef _VertElem_h

#define _VertElem_h 1

/* Provides proper sementics for adding and deleting vertices from rings.
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

   $Id: VertElem.h,v 1.7 2000/05/27 01:25:13 ttnttn Exp $
*/

#include "MEAD/Vertex.h"

class VertElem {

  public :

  VertElem ( const Coord& c );
  VertElem ( const VertElem& v );
  VertElem& operator = ( const VertElem& v );
  ~VertElem( );

  float angle;
  Vertex *vpt;
  VertElem *next;
  VertElem *prev;
#ifdef POPULATION_COUNT
  static population_report() {
    cout <<   "VertElem current population = " << population;
    cout << "\n             number created = " << num_created;
    cout << "\n             number deleted = " << num_deleted;
    cout << "\n  maximum population so far = " << max_population << endl;
  }
#endif

  private :
  // for debugging FIXME!
#ifdef POPULATION_COUNT
  static int population;
  static int max_population;
  static int num_created;
  static int num_deleted;
#endif

  VertElem () : next(0),angle(0),vpt(0) { };
};

inline VertElem :: ~VertElem( )
{
  // cout << "deleting VertElem with angle = " << angle << endl;
  if( --vpt->count <= 0 )
    delete vpt;
#ifdef POPULATION_COUNT
  ++num_deleted;
  --population;
#endif
}

inline VertElem :: VertElem( const Coord& c )
{
  vpt = new Vertex (c);
  next = 0;
  prev = 0;
#ifdef POPULATION_COUNT
  ++num_created;
  ++population;
  if (population > max_population)
    max_population = population;
#endif
}

inline VertElem :: VertElem( const VertElem& v )
{
  vpt = v.vpt;
  vpt -> count++;
  next = v.next;
  prev = v.prev;
#ifdef POPULATION_COUNT
  ++num_created;
  ++population;
  if (population > max_population)
    max_population = population;
#endif
}

// NOTE: Assignment does not alter next and prev pointers so
// position in linked list remains that of original l.h.s.
inline VertElem& VertElem :: operator = ( const VertElem& v )
{
  v.vpt->count++;
  if( --vpt->count <= 0 )
    delete vpt;
  vpt = v.vpt;
  return *this;
}


#endif

// VertElem.h ends here
