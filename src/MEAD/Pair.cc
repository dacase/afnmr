/* A Pair of shells, including the associated ring of vertexes.
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

   $Id: Pair.cc,v 1.10 2004/12/10 20:36:11 bashford Exp $
*/


#include "MEAD/Pair.h"
#include <math.h>


/* For each pair, there is a ring of vertexes along the overlapping
   circle.  Each triplex for which a solution to spcor3 exists gives
   rise to a pair vertexes which cause some section of the ring to
   become "closed" (vertexes sterically disallowed in such regions)
   while leaving the possibility of the rest of the ring being "open."

   We maintain a linked list of vertexes which are not yet dissallowed
   by other vertexes for the Pair.  Each vertex is a boundary between
   an open and a closed region.  By convention, the list starts with a
   vertex that begins an open region.

   Check_vertex_pair is the main function responsible for maintaining
   this list.  Given a new pair it makes the proper adjustments to the
   list, possibly adding one or both new vertexes, possibly deleting
   existing vertexes.
   */


  // For debugging purposes...
#ifdef DEBUG
#define VERTLIST_SANITY_CHECK  if (vert_head) {				   \
    VertElem *cur = vert_head;						   \
    if (cur->prev) {							   \
      cerr << "ERROR: vert_head->prev non-zero on call # " << call_number  \
	<< endl;							   \
    }									   \
    float ang = cur->angle;						   \
    VertElem *prev = cur;						   \
    cur = cur->next;							   \
    int vert_count = 1;							   \
    for ( ; cur; cur=cur->next) {					   \
      if (prev != cur->prev) {						   \
	cerr << "ERROR: prev pointer screwed up on call # " << call_number \
	  << endl;							   \
      }									   \
      if (cur->angle < ang) {						   \
	cerr << "ERROR: angles out of order on call # " << call_number	   \
	  << endl;							   \
      }									   \
      ++vert_count;							   \
      prev = cur;							   \
      ang = cur->angle;							   \
    }									   \
    if (vert_count%2) {							   \
	cerr << "ERROR: odd number of vertexes on call # " << call_number  \
	  << endl;							   \
      }									   \
  }
#endif

void
Pair::Check_vertex_pair(const VertElem& vert1, const VertElem& vert2)
{
  static int call_number = 0;
  ++call_number;

  /* On entry, vert1 is known to be the HEAD, that is, the start, in a
     postive angle sense, of an open region.  vert2 is the TAIL, that
     is, the end of this open region.  */

  if( vert_head == 0 ) {             // List is empty so add the pair.
    ref_vector = vert1.vpt->coord - overlap_center;
    VertElem *v1 = new VertElem(vert1);
    VertElem *v2 = new VertElem(vert2);
    v1->angle = 0;
    v2->angle = vertex_angle(vert2);
    vert_head = v1;
    v1->prev = 0;
    v1->next = v2;
    v2->prev = v1;
    v2->next = 0;
    flag = partially_buried;

#ifdef DEBUG
  VERTLIST_SANITY_CHECK
#endif

    return;
  }

  /* List is non-emty so scan it to see where the new vertices can go:

     Set current1 to point to the first current vertex whose angle is
     greater than than the angle of new head vertex and set current2
     to point to the first current vertex whose angle is greater than
     than the angle of new tail vertex.  That is, current1 or 2 will
     end up pointing to the vertex following the region conaining the
     head or tail, respectively.  However, if either the head or tail
     angle is greater than the largest current angle (if head or tail
     is in the last region) then the corresponding current# will be
     null.

     Also, the indices n1 and n2 give the ordinality, starting with
     the first region after vert_head being number 1, of the regions
     designated by current1 and current2, respectively.  Since the
     regions go: open, closed, open, ...; an even or odd index
     indicates a closed or open region, respectively.

     There is a special case however, the space between the ref_vector
     and vert_head (such a space may occur when the original vert_head
     is replaced) is intially numbered zero, but it is really the same
     region as the last region which is given a high (and always even)
     index.  This is harmless!  Since this is always a closed region,
     the correct behavior when a new vertex falls in either section
     of this region is to just delete other vertices according to
     whether the new vertexes invalidate them.

     There is a tricky question of relating to wether angle checking
     should be done with inequalites like > or >=.  It is the question
     of whether a vertex which forms a region boundary "belongs" to
     one region or the other.  I (D. B.) have decided that the
     vertexes must belong to the OPEN region that they bound.  I think
     this convention is necessery otherwise a vertex would invalidate
     itself.  (It would fail to belong to its "own" open region.)
     this is the reason for the code:

         (i%2 ? (cur->angle>=head_ang) : (cur->angle>head_ang))

     in the conditions below.  i%2 is the condition that the region
     before cur->angle is an open one; if so a more inclusive
     condition is used.

  */

  float head_ang = vertex_angle(vert1);
  float tail_ang = vertex_angle(vert2);

  VertElem *current1 = 0;
  VertElem *current2 = 0;
  int n1 = 0;
  int n2 = 0;
  int hflag = 0;
  int tflag = 0;
  int i;
  VertElem *cur;
  for(cur = vert_head, i=0; cur; cur=cur->next, ++i) {
    if (( hflag==0)
	&& (i%2 ? (cur->angle>=head_ang) : (cur->angle>head_ang))) {
      n1 = i;
      hflag = 1;
      current1 = cur;
    }
    if( tflag==0
	&& (i%2 ? (cur->angle>=tail_ang) : (cur->angle>tail_ang))) {
      n2 = i;
      tflag = 1;
      current2 = cur;
    }
    if (tflag && hflag) break;
  }
  if (current1==0)
    n1 = i;
  if (current2==0)
    n2 = i;

  // Adjust the linked list (see below)...
  if( n1==n2 ) {
    if( n1%2 ) {
      /* Here, both vertices are in the same open region so both must
	 be added to the list.  If head_ang > tail_ang, a small new
	 closed region is to be formed inside this open region, so
	 insert two new elements into list and leave the rest of the
	 list alone.  But if head_ang < tail_ang, a small new open
	 region is to REPLACE this old open region and will become the
	 ONLY open region, so replace entire list with these two new
	 elements. */
      VertElem *v1 = new VertElem(vert1);
      VertElem *v2 = new VertElem(vert2);
      v1->angle = head_ang;
      v2->angle = tail_ang;

      if( head_ang>tail_ang ) {
	VertElem*  current = current1->prev;
	current->next = v2;
	v2->prev = current;
	v2->next = v1;
	v1->prev = v2;
	v1->next = current1;
	current1->prev = v1;
      }
      else {
	delete_vertex_ring();
	vert_head = v1;
	v1->prev = 0;
	v1->next = v2;
	v2->prev = v1;
	v2->next = 0;
      }
    }

    /* Both vertices in the same closed region.  if head_ang <
       tail_ang, it means that that everything outside this region is
       closed and since region itself is already closed, this whole
       Pair is closed (i.e. a buried Pair).  */

       else {
	 if( head_ang<tail_ang ) {
	   delete_vertex_ring();
	 }
       }
  }
  else                   //if two vertexes are in different regions
    {
      /* vert1 in an open region different from vert2 means that the
	 old head of this region needs to be replaced by vert1.
	 Likewise, vert2 in an open region different from vert1 means
	 that the old tail of this region needs to be replaced by
	 vert2.  A subsequent cleanup will take care of other old
	 vertices now invalidated.
	 */
      /* This would seem to lose of current1 is zero, but in that
         case, we're in the last region, which is odd numbered and closed
         (provided ring is sane), so the problem shouldn't happen. */
      if(n1%2) {
	*current1->prev = vert1;
	current1->prev->angle = head_ang;
      }
      if(n2%2) {
	*current2 = vert2;
	current2->angle = tail_ang;
      }
      /* Here is the "subsequent cleanup" mentionded above.  The
	 vertex list is scanned again and any vertexes outside the
	 range from the new head to the new tail are deleted.  The
	 exact criteria are a bit different according to whether
	 head_ang < tail_ang.  */

      if( head_ang < tail_ang ) {
	/* This for loop and the one below must be carfully done so
	   that it can continue sensibly even after deletion of cur
	   during the loop.  BUG PRONE!
	 */
	VertElem *cur = vert_head;
	while(cur) {
	  float ang = cur->angle;
	  if( ang<head_ang || ang>tail_ang ) {
	    if (cur->prev)
	      cur->prev->next = cur->next;
	    else
	      vert_head = cur->next;
	    if (cur->next)
	      cur->next->prev = cur->prev;
	    VertElem *doomed_cur = cur;
	    cur = cur->next; // iteration ahead of the delete
	    delete doomed_cur;
	    continue;
	  }
	  cur = cur->next; // normal iteration
	}
      }
      else {
	VertElem *cur = vert_head;
	while(cur) {
	  float ang = cur->angle;
	  if( ang<head_ang && ang>tail_ang ) {
	    if (cur->prev)
	      cur->prev->next = cur->next;
	    else
	      vert_head = cur->next;
	    if (cur->next)
	      cur->next->prev = cur->prev;
	    VertElem *doomed_cur = cur;
	    cur = cur->next; // iteration ahead of the delete
	    delete doomed_cur;
	    continue;
	  }
	  cur = cur->next; // normal iteration
	}
      }
    }

  if( vert_head==0 )
    flag = buried;
  else
    flag = partially_buried;

#ifdef DEBUG
  VERTLIST_SANITY_CHECK
#endif

  return;
}

// Pair.cc ends here
