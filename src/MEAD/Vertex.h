// This is -*- C++ -*-
#ifndef _Vertex_h
#define _Vertex_h 1

/* A vertex occurs where three shells intersect.
   $Id: Vertex.h,v 2.3 2000/05/27 01:25:11 ttnttn Exp $
*/

#include "MEAD/Coord.h"

class Vertex{

  public :

  inline Vertex( const Coord& c ) { coord=c; count=1;
#ifdef POPULATION_COUNT
				    ++num_created;
				    ++population;
				    if (population > max_population)
				      max_population = population;
#endif
				  };
#ifdef POPULATION_COUNT
  static population_report() {
    cout <<   "Vertex current population = " << population;
    cout << "\n           number created = " << num_created;
    cout << "\n           number deleted = " << num_deleted;
    cout << "\n   max. population so far = " << max_population << endl;
  }
  ~Vertex() { ++num_deleted;
	      --population;
	    };
#endif
  Coord coord;
  int   count;

  private :
  Vertex( ) { };

#ifdef POPULATION_COUNT
  static int population;
  static int max_population;
  static int num_created;
  static int num_deleted;
#endif
};

#endif

// Vertex.h ends here
