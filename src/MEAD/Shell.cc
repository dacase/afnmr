/* Shell is a sphere with inner and outer radius and a set of Cones
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

$Id: Shell.cc,v 1.7 2004/12/06 17:58:37 bashford Exp $
*/

#include <string>
#include "MEAD/Shell.h"


Shell::Shell (const Coord& a, float inner_r, float probe_r)
{
  coord = a;
  inner_rad = inner_r;
  outer_rad = inner_rad + probe_r;
  outer_rad_sq = outer_rad*outer_rad;
  numcones = 0;
  flag = free;
  cone_head = 0;
}


Shell::Shell (const Atom& a, float probe_r) {
  coord = a.coord;
  inner_rad = a.rad;
  outer_rad = inner_rad + probe_r;
  outer_rad_sq = outer_rad*outer_rad;
  flag = free;
  numcones = 0;
  cone_head = 0;
}


Shell::Shell (const Shell& a)
{
  coord = a.coord;
  inner_rad = a.inner_rad;
  outer_rad = a.outer_rad;
  outer_rad_sq = a.outer_rad_sq;
  flag = a.flag;
  numcones = a.numcones;
  if (a.cone_head) {
    cone_head = new Cone(*a.cone_head);
    Cone *tc, *ac;
    for (tc = cone_head, ac = a.cone_head;
	 ac->next; ac=ac->next, tc=tc->next)
      tc->next = new Cone (*ac->next);
    tc->next=0;
  }
}

Shell&
Shell::operator = (const Shell& a)
{
  coord = a.coord;
  flag = a.flag;
  inner_rad = a.inner_rad;
  outer_rad = a.outer_rad;
  outer_rad_sq = a.outer_rad_sq;
  numcones = a.numcones;
  if (a.cone_head) {
    if (cone_head && cone_head != a.cone_head)
      delete_conelist();
    cone_head = new Cone(*a.cone_head);
    Cone *tc, *ac;
    for (tc = cone_head, ac = a.cone_head;
	 ac->next; ac=ac->next, tc=tc->next)
      tc->next = new Cone(*ac->next);
    tc->next=0;
  }
  return *this;
}

void
Shell::add_cone(const Cone& c)
{
  if (flag == partially_buried) {
    Cone *nc = new Cone(c);
    nc->next = cone_head;
    cone_head = nc;
    ++numcones;
  }
  else
    ::error("Shell::add_cone: operation allowed only if partially buried");
}

void
Shell::delete_conelist()
{
  for (Cone *tc = cone_head; tc; ) {
    Cone *tmp = tc->next;
    delete tc;
    tc = tmp;
  }
  cone_head=0;
  numcones=0;
}

ostream& Shell::write_top_in_binary(ostream& bfile)
{
 float *topology_data = new float[4];
 if( !bfile )
   {
    cerr<<"Shell::write_top_in_binary : can not open the file "<<endl;
    return bfile;
   }
 else
   {
    topology_data[0] = coord.x;
    topology_data[1] = coord.y;
    topology_data[2] = coord.z;
    topology_data[3] = inner_rad;
    bfile.write((char *) topology_data, sizeof(float)*4 );
    bfile.write((char *) &outer_rad, sizeof(float) );
    bfile.write((char *) &flag, sizeof(int) );
    bfile.write((char *) &numcones, sizeof(int) );
    Cone *tmp = cone_head;
    for( int m=0; m<numcones; m++ )
       {
        topology_data[0] = tmp -> unit_axis.x;
        topology_data[1] = tmp -> unit_axis.y;
        topology_data[2] = tmp -> unit_axis.z;
        topology_data[3] = tmp -> cos_ang1_sq;
        bfile.write((char *) topology_data, sizeof(float)*4 );
        tmp = tmp -> next;
       }
   }
 delete topology_data;
 return bfile;
}

ostream& Shell::write_top_in_ascii(ostream& afile)
{
 if( !afile )
   {
    cerr<<"Shell::write_top_in_ascii : can not open the file "<<endl;
    return afile;
   }
 else
   {
    Cone *tmp = cone_head;
    coord.print(afile);
    afile << "  " << inner_rad << "  " << outer_rad << "  " << flag
          << "  " << numcones << endl;
    for( int m=0; m<numcones; m++ )
       {
        Coord cone_axis = tmp -> unit_axis;
        cone_axis.print(afile);
        afile << "  ";
        afile << tmp -> cos_ang1_sq << endl;
        tmp = tmp -> next;
       }
   }
 return afile;
}

int Shell::read_top_in_binary(istream& bfile)
{
 int num_cones;
 float top_data[4];
 bfile.read( (char*) top_data, sizeof(float)*4 );
 if( !bfile )
   {
    cerr << "ERROR : Shell::read_top_in_binary: " << endl;
    return 0;
   }
 coord.x = top_data[0];
 coord.y = top_data[1];
 coord.z = top_data[2];
 inner_rad = top_data[3];
 bfile.read( (char*) &outer_rad, sizeof(float) );
 if( !bfile )
   {
    cerr << "ERROR : Shell::read_top_in_binary: 1" << endl;
    return 0;
   }
 outer_rad_sq = outer_rad*outer_rad;
 bfile.read( (char*) &flag, sizeof(int) );
 if( !bfile )
   {
    cerr << "ERROR : Shell::read_top_in_binary: 2" << endl;
    return 0;
   }
 bfile.read( (char*) &num_cones, sizeof(int) );
 if( !bfile )
   {
    cerr << "ERROR : Shell::read_top_in_binary: 3" << endl;
    return 0;
   }
 for( int n=0; n<num_cones; ++n )
    {
     bfile.read( (char*) top_data, sizeof(float)*4 );
     if( !bfile )
       {
        cerr << "ERROR : Shell::read_top_in_binary: 4" << endl;
        return 0;
       }
     Cone new_cone;
     new_cone.unit_axis.x = top_data[0];
     new_cone.unit_axis.y = top_data[1];
     new_cone.unit_axis.z = top_data[2];
     new_cone.cos_ang1_sq = top_data[3];
     add_cone( new_cone );
    }
  return 1;
}

int Shell::read_top_in_ascii(istream& afile)
{
 string word;
 int shell_number;
 int flag_numb, num_cones;
 afile >> word;
 afile >> shell_number;
 coord.read(afile);
 afile >> inner_rad;
 afile >> outer_rad;
 outer_rad_sq = outer_rad*outer_rad;
 afile >> flag_numb;
 afile >> num_cones;
 if( flag_numb == 0 )
     flag = free;
 else
   if( flag_numb == 1 )
       flag = partially_buried;
   else
     if( flag_numb == 2 )
       flag = buried;
     else
       {
        cerr << "ERROR : Shell::read_top_in_ascii: "
             << "wrong flag type detected " << endl;
        return 0;
       }
 for( int n=0; n<num_cones; n++ )
    {
     Cone new_cone;
     new_cone.unit_axis.read(afile);
     afile >> new_cone.cos_ang1_sq;
     add_cone( new_cone );
    }
 return 1;
}

// Shell.cc ends here
