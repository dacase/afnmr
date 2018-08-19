#include <iostream>
using std::ios;
#include <sstream>
#include <string>
#include "MEAD/Coord.h"
#include "MEAD/CubeLatSpec.h"
#include "MEAD/AvsScalarField.h"

// Some C standard libraries:
#include <string.h>
#include <ctype.h>
#include <time.h>

// FIXME. SGI's iostream library has ios::setstate as a non-public function,
// contrary to the Standard.  Therefore, calls to setstate(ios::failbit)
// have been replaced by clear(ios::failbit) as a work-around.
// This is not so nice since it blows away previous status bits.

// Function to look for something line "name =", store name in strng,
// and return istream positioned just after the "="
istream&
getname(istream& istr, string& strng)
{
  istr >> std::ws; // Skip leading whitespace
  char c;
  string ts;
  while (istr.get(c)) {  // Scan for the name
    if (isspace(c) || c == '=') break;
    ts += c;
  }
  if (!istr) return istr;
  if (c != '=') {
    // Skip any whitespace, then expect the '='
    istr >> std::ws;
    istr.get(c);
    if (c != '=' && istr) {
      istr.clear(ios::failbit);
      return istr;
    }
  }
  strng = ts;
  return istr;
}

/*
A big ugly function which parses an AVS header, checks whether data
in it are consistent with CubeLatSpec previously given, and decides
whether field type is uniform or rectilinear, and whether input mode
is native or data parsing.  In the data parsing case, store info in
file names, offsets, etc.
(It may be better to use sstream and string if generally available. FIXME?)
*/
istream&
AvsScalarField::read_header(istream& head)
{
  const int linebufsize = 500;
  char linebuf[linebufsize];

  int ndim = 0;  // Must appear and must be 3
  int dim1 = 0;  // Must appear and must be grid_dim.
  int dim2 = 0;  // ditto
  int dim3 = 0;  // ditto
  int veclen = 0; // Must appear and must be 1.
  int nspace = 0;  // If it appears, it must be three

  fieldtype = undefined;
  min_ext = Coord(0.0, 0.0, 0.0); // Extents needed for uniform fields
  min_ext_defined = false;
  max_ext = Coord(0.0, 0.0, 0.0);
  max_ext_defined = false;
  separator_seen = false;

  // Intialize some data parsing mode params to "undefined" state;
  var_stride = 0;
  coord_stride[0] = coord_stride[1] = coord_stride[2] = 0;

  for (int linenum=1; head; ++linenum) {
    linebuf[0] = '\0';
    head.get(linebuf, linebufsize, '\n');
    int linelength = head.gcount();
    if (!head) break;
    if (head.get() != '\n') {
      cerr << "WARNING: AvsScalarField::read_header encountered line > "
	<< (linebufsize-1) << " characters long\n"
	  << "Will try to interpret anyway" << endl;
    }

    if (linenum == 1) {  // The first line must begin with "# AVS"
      if (! (linelength >= 5) && (strncmp("# AVS", linebuf, 5) == 0)) {
	cerr << "ERROR:  AvsScalarField::read_header:  "
	  << "Input does not begin with \"# AVS\"." << endl;
	head.clear(ios::failbit);
	return head;
      }
      continue;
    }

    // Throw away comments and leading whitspaces, ignore empty lines
    char * comm = strchr(linebuf, '#');
    if (comm) {
      *comm = '\0';
      linelength = comm - linebuf;
    }
    char *pbuf = linebuf;
    for ( ; isspace(*pbuf); ++pbuf, --linelength) ;
    if (linelength==0) continue;

    // Look for data parsing input mode "variable" or "coordinate" lines
    enum {var, coord, neither} dptype = neither;
    if (strncmp(pbuf, "variable", 8) == 0) dptype = var;
    if (strncmp(pbuf, "coord", 5) == 0) dptype = coord;
    if (dptype == var || dptype == coord) {
      string filespec, filetype;
      int skip=0;
      int stride=1;
      int n;
      string strng;
      std::istringstream dpstr(pbuf);
      dpstr >> strng;  // Consume the "varable" or "coord" word
      dpstr >> n; // This is required.
      if (dpstr.fail()) {
	cerr << "ERROR:  AvsScalarField::read_header:  "
	  << "A \"varable\" or \"coord\" is missing \"n\" specifier." << endl;
	head.clear(ios::failbit);
	return head;
      }
      while (getname(dpstr,strng)) {
	if (strng == string("file"))
	  dpstr >> filespec;
	else if (strng == string("filetype"))
	  dpstr >> filetype;
	else if (strng == string("skip"))
	  dpstr >> skip;
	else if (strng == string("stride"))
	  dpstr >> stride;
	else if (strng == string("offset")) {
	  cerr << "WARNING:  AvsScalarField::read_header: offset specifier\n"
	    << "ignored since ASCII data parsing not supported" << endl;
	  int junk;
	  dpstr >> junk;
	}
	else {
	  cerr << "WARNING:  AvsScalarField::read_header:\n"
	    << "Data parsing mode specifier, " << strng
	      << ", not recognized" << endl;
	  // Consume something and hope it is OK to continue..
	  string junk;
	  dpstr >> junk;
	}
      }
      // Done reading a data-parsing line...
      // check and/or store values read from it.
      if (!(filetype == "" || filetype == string("binary"))) {
	cerr << "ERROR:  AvsScalarField::read_header: Sorry,\n"
	  << "data parsing filetypes other than \"binary\" not supported."
	    << endl;
	head.clear(ios::failbit);
	return head;
      }
      if (dptype == var) {
	if (n != 1) {
	  cerr << "ERROR:  AvsScalarField::read_header:  "
	    << "varible n specifier != 1."
	      << endl;
	  head.clear(ios::failbit);
	  return head;
	}
	var_filespec = filespec;
	var_skip = skip;
	var_stride = stride;
      }
      else if (dptype == coord) {
	if (n < 1 || n > 3) {
	  cerr << "ERROR:  AvsScalarField::read_header:  "
	    << "coord n specifier not in range 1..3"
	      << endl;
	  head.clear(ios::failbit);
	  return head;
	}
	coord_filespec[n-1] = filespec;
	coord_skip[n-1] = skip;
	coord_stride[n-1] = stride;
      }
      continue;
    }  // Finished processing a data-parsing line

    // Process lines with a single NAME=VALUE pair;
    std::istringstream nvstr(pbuf);
    string name;
    if (getname(nvstr, name)) {
      blab3 << "Processing name=value pair:" << pbuf << endl;
      if (name == string("ndim")) {
	nvstr >> ndim;
	blab3 << "READ: ndim = " << ndim << endl;
      }
      else if (name == string("dim1")) {
	nvstr >> dim1;
      }
      else if (name == string("dim2")) {
	nvstr >> dim2;
      }
      else if (name == string("dim3")) {
	nvstr >> dim3;
      }
      else if (name == string("nspace")) {
	nvstr >> nspace;
      }
      else if (name == string("veclen")) {
	nvstr >> veclen;
      }
      else if (name == string("data")) {  // Value must be "float"
	string data;
	nvstr >> data;
	if (data != string("float")) {
	  cerr << "WARNING AvsScalarField::read_header: "
	    << "data type must be float" << endl;
	  head.clear(ios::failbit);
	  return head;
	}
      }
      else if (name == string("field")) {
	string field;
	nvstr >> field;
	if (field == string("uniform")) fieldtype = uniform;
	else if (field == string("rectilinear")) fieldtype = rectilinear;
	else {
	  cerr << "WARNING AvsScalarField::read_header: "
	    << "field type, " << field << ", not supported." << endl;
	  head.clear(ios::failbit);
	  return head;
	}
      }
      else if (name == string("min_ext")) {
	nvstr >> min_ext;
	min_ext_defined = true;
      }
      else if (name == string("max_ext")) {
	nvstr >> max_ext;
	max_ext_defined = true;
      }
      else if (name == string("min_val")) ; // Ignore
      else if (name == string("max_val")) ; // Ignore
      else if (name == string("labels")) ; // Ignore
      else if (name == string("unit")) ; // Ignore
      else {
	cerr << "WARNING: AvsScalarField::read_header: Don't recognize \""
	  << name << "\" in input line,\n" << linebuf << endl;
      }
    } // end of "if this is a name=value pair"

    // Look out for the separator before the binary area;
    if (head) {
      int ch = head.get();
      if (ch == 12) { // The ASCII form-feed character.  There should be two.
	ch = head.get();
	if (ch != 12) {
	  cerr << "WARNING AvsScalarField::read_header: "
	    << "Expected a second form-feed in separator, but did not find it."
	      << endl;
	}
	separator_seen = true;
	break;
      }
      else
	head.putback((char) ch);
    }
  }

    // Done parsing ASCII header.  Check results so far.

  if (ndim != 3) {
    cerr << "ERROR AvsScalarField::read_header: ndim != 3." << endl;
    head.clear(ios::failbit);
    return head;
  }
  if (!(dim1 == grid_dim && dim2 == grid_dim && dim3 == grid_dim)) {
    cerr << "ERROR AvsScalarField::read_header: Not all dims == grid_dim."
      << endl;
    head.clear(ios::failbit);
    return head;
  }
  if (veclen != 1) {
    cerr << "ERROR AvsScalarField::read_header: veclen != 1." << endl;
    head.clear(ios::failbit);
    return head;
  }
  if (!(nspace == 0 || nspace == 3)) {
    cerr << "ERROR AvsScalarField::read_header: nspace defined, but not as 3."
      << endl;
    head.clear(ios::failbit);
    return head;
  }
  switch (fieldtype) {
  case undefined:
    cerr << "ERROR AvsScalarField::read_header: field type undefined."
      << endl;
    head.clear(ios::failbit);
    return head;
  case uniform:
    if ((!separator_seen) && !(min_ext_defined && max_ext_defined)) {
      cerr << "WARNING AvsScalarField::read_header: "
	<< "min and/or max extent data not in header\n"
	  << "of a uniform/data-formatting input\n"
	    << "This means no check of agreement with CubeLatSpec" << endl;
    }
    if ((min_ext_defined && max_ext_defined)
	&& !(min_ext == bottom_corner_in_space
	     && max_ext == top_corner_in_space)) {
      cerr << "ERROR AvsScalarField::read_header:\n"
	<< "min and/or max extent data does not agree with CubeLatSpec."
	  << endl;
      head.clear(ios::failbit);
      return head;
    }
    break;
  case rectilinear:
    // Nothing to be done.  min_ext and max_ext to be ignored.
    break;
  }

  return head;
}

// AvsSF_read_header.cc ends here
