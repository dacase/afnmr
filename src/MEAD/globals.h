/* Definition of some global constants, variables and functions   -*- C++ -*-
    Copyright (c) 1993--1995 by Donald Bashford.

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

$Id: globals.h,v 2.11 2007/11/26 22:39:32 bashford Exp $
*/

#ifndef _globals_h
#define _globals_h 1
#include <string>
using std::string;
#include <iostream>
#include <fstream>

using std::ostream;
using std::istream;
using std::ofstream;
using std::ifstream;
using std::cout;
using std::cerr;
using std::cin;
using std::endl;

extern ofstream cnull;
extern ostream *blab1pt;
extern ostream *blab2pt;
extern ostream *blab3pt;
#define blab1 (*blab1pt)
#define blab2 (*blab2pt)
#define blab3 (*blab3pt)

extern const float pi;
extern const float ln10;

extern void error (const string& s1 = "", const string& s2 = "",
		   const string& s3 = "");

// open a new file named by the string for the ofstream,
// but call error rather than do the open if the named file already exists.
extern void safeopen(ofstream&, std::string);

extern std::string DOT;

#endif

// globals.h ends here
