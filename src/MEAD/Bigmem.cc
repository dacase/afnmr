/* For multiple allocate-free cycles of big arrays of the same size.

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

$Id: Bigmem.cc,v 2.6 2005/01/25 03:23:55 bashford Exp $
*/

/*
big_rigid_malloc and big_rigid_free are for big arrays that are
repeatedly allocated and freed in a recurring cycle, always with the
same size.  To use them for anything else is a sure road to hell.  They
are really dumb.  They are intended to solve a memory fragmentation
problem that caused sbrks to go on and memory requirements to grow
even though the total memory malloced remained the same from cycle to
cycle.
*/

#include "MEAD/globals.h"
#include <iostream>

struct MemNode {
  unsigned size;
  int inuse;
  char* ptr;
  MemNode* bigger;
};

static MemNode* biggest = 0;
static unsigned totalsize = 0;
static unsigned numnodes = 0;

void * big_rigid_malloc (unsigned sz)
{
#ifdef BIGMEM_DEBUG
  cout << "big_rigid_malloc called with sz = " << sz << endl;
#endif
  if (sz < 50000)
    blab3 << "ADVICE: big_rigid_malloc: " << sz << " isn't so big,\n"
      << "maybe you should user regular malloc instead.\n"
	<< "(Proceeding with big_malloc_anyway" << endl;
  if (biggest != 0) {
    MemNode* p = biggest->bigger;
    for (;;) {
      if (p->size > sz) break;
      if (!p->inuse && p->size == sz) {
	// We'll take it!
	p->inuse = 1;
#ifdef BIGMEM_DEBUG
  cout << "big_rigid_malloc returning pointer " << p->ptr << endl;
#endif
	return p->ptr;
      }
      if (p == biggest) break;
      p = p->bigger;
    }
  }
  // Fooey!  Didn't find a space.  Create one, 5% larger than needed.
  MemNode *n = new MemNode;
  n->ptr = new char[sz];
  n->size = sz;
  n->inuse = 1;
  if (biggest != 0) {
    MemNode* p = biggest->bigger;  // that is, the smallest
    MemNode* prev = biggest;
    for (;;) {
      if (p->size > sz) {
	//Insert it here to keep list size ordered
	prev->bigger = n;
	n->bigger = p;
	break;
      }
      if (p == biggest) {
        // The new one will be the biggest
	n->bigger = p->bigger;
	p->bigger = n;
	biggest = n;
	break;
      }
      prev = p;
      p = p->bigger;
    }
  }
  else
    biggest = n->bigger = n;

  totalsize += sz;
  ++numnodes;
#ifdef BIGMEM_DEBUG
  cout << "big_rigid_malloc returning pointer " << n->ptr << endl;
#endif
  return n->ptr;
}

void big_rigid_free (void* pt)
{
#ifdef BIGMEM_DEBUG
  cout << "big_rigid_free called with pt = " << pt << endl;
#endif
  if (pt == 0) return;  // analogue of harmless "delete 0"
  if (!biggest)
    error ("big_rigid_free: There are no nodes allocated yet, can't free\n");
  MemNode* p = biggest->bigger;
  for (;;) {
    if (p->ptr == pt) {
      if (p->inuse) {
	p->inuse=0;
	return;
      }
      else
	error ("big_rigid_free: This node was already freed\n");
    }
    if (p == biggest) break;
    p = p->bigger;
  }
  error ("big_rigid_free:  Couldn't find this address in the list\n");
}

void big_rigid_delete_all()
{
#ifdef BIGMEM_DEBUG
  cout << "big_rigid_delete_all called" << endl;
#endif
  if (biggest == 0)
    return;
  MemNode *next_victim = biggest->bigger;
  while (next_victim) {
    MemNode *p = next_victim;
    if (next_victim==biggest)
      next_victim = 0;
    else
      next_victim = p->bigger;
    p = p->bigger;
    if (p->inuse)
      ::error("ERROR: big_rigid_delete_all called with a block still in use");
    delete [] p->ptr;
    delete p;
  }
  biggest = 0;
}

#ifdef MALLOC_STATS

void big_rigid_stats ()
{
  if (biggest != 0) {
    cout << numnodes  << "\nBlocks currently maintained by Bigmem code:"
      << endl;
    unsigned long int sused = 0;
    unsigned long int sfree = 0;
    MemNode* p = biggest->bigger;
    for (;;) {
      cout << "Block of size " << p->size;
      if (p->inuse) {
	sused += p->size;
	cout << " is in use" << endl;
      }
      else {
	sfree += p->size;
	cout << " is free" << endl;
      }
      if (p == biggest) break;
      p = p->bigger;
    }
    cout << "\nTotal free block size = " << sfree
      << "\nTotal used block size = " << sused
	<< "\nGrand total under Bigmem control = " << (sused + sfree) << endl;
  }
  else
    cout << "No memory currently maintained by Bigmem code" << endl;
}

#endif

// Bigmem.cc ends here
