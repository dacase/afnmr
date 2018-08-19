#include <iostream>
using std::ios;
#include <string>
#include "MEAD/Coord.h"
#include "MEAD/CubeLatSpec.h"
#include "MEAD/AvsScalarField.h"
#include "MEAD/Bigmem.h"

// Some C standard libraries:
#include <string.h>
#include <ctype.h>
#include <time.h>

// FIXME. SGI's iostream library has ios::setstate as a non-public function,
// contrary to the Standard.  Therefore, calls to setstate(ios::failbit)
// have been replaced by clear(ios::failbit) as a work-around.
// This is not so nice since it blows away previous status bits.

bool
AvsScalarField::read(const CubeLatSpec& cls, const string& fieldname)
{
  grid_dim = cls.get_grid_dim();
  spacing = cls.get_spacing();
  grid_center_in_space = cls.get_center();
  float grlen = (float) (grid_dim - 1);
  float halfgrlen = grlen/2;
  grid_center_in_grid = Coord(halfgrlen, halfgrlen, halfgrlen);
  float halfsplen = spacing * halfgrlen;
  bottom_corner_in_space = grid_center_in_space
    - Coord(halfsplen, halfsplen, halfsplen);
  top_corner_in_space = grid_center_in_space
    + Coord(halfsplen, halfsplen, halfsplen);

  // Open .fld (header info) and .avsdat (binary) files based on fieldname
  string headname_string = fieldname + ".fld";
  const char *headname = headname_string.c_str();
  ifstream head(headname, ios::in);
  if (!head.good()) {
    cerr << "ERROR: AvsScalarField::read: File, " << headname
      << " could not be opened" << endl;
    return 0;
  }
  read_header(head);
  blab3 << "In AvsScalarField::read: read_header has returned\n";
  if (separator_seen) { // This is an "native" style AVS field file
    read_binary_field(head);
    if (fieldtype == rectilinear) {
      read_rectilinear_coord(head, 1);
      read_rectilinear_coord(head, 2);
      read_rectilinear_coord(head, 3);
    }
    else
      read_uniform_coord(head);
  }
  else { // The header was from a "data parsing mode" field file
    blab3 << "opening var_filespec = \"" << var_filespec << "\"" << endl;
    ifstream data(var_filespec.c_str(), ios::in | ios::binary);
    if (!data) {
      cerr << "ERROR: AvsScalarField::read: Could not open file, "
	<< var_filespec << endl;
      return 0;
    }
    read_binary_field(data);
    blab3 << "In AvsScalarField::read: read_binary_field has returned"
      << endl;;
    if (fieldtype == rectilinear) {
      string prev_filespec = var_filespec;
      for (int i=0; i<3; ++i) {
	if (coord_filespec[i] != prev_filespec) {
	  data.close();
	  data.open(coord_filespec[i].c_str(), ios::in | ios::binary);
	}
	read_rectilinear_coord(data, i+1);
	blab3 << "In AvsScalarField::read: read_rectilinear(data, "
	  << (i+1) << ") has returned" << endl;;
	prev_filespec = coord_filespec[i];
      }
    }
  }
  return 1;
}

ifstream&
AvsScalarField::read_binary_field(ifstream& data)
{
  int ncube = grid_dim*grid_dim*grid_dim;
  float * data_array = (float*) big_rigid_malloc(ncube*(sizeof (float)));

  blab3 << "in read_binary_field with var_stride = "
    << var_stride << ", var_skip = " << var_skip
      << ", ncube = " << ncube << endl;

  // if stride is non-zero, we are under data parsing mode so use
  // and skip, stride to control reading.  Otherwise assume data is
  // positioned to the start of a sequence
  if (var_stride) {
    if (var_stride < 0 || var_skip < 0) {
      cerr << "ERROR: AvsScalarField::read_binary_field:\n"
	<< "either var_stride = " << var_stride << " or var_skip = "
	  << var_skip << " is nonsensical" << endl;
      data.clear(ios::failbit);
      return data;
    }
    data.seekg(var_skip, ios::beg);
    if (var_stride==1) {
      blab3 << "Attempting read of " << ncube << " floats" << endl;
      data.read((char *) data_array, (sizeof (float)) * ncube);
    }
    else {
      for (int i=0; i<ncube-1; ++i) {
	data.read((char *) &data_array[i], sizeof (float));
	data.seekg(var_stride * (sizeof (float)), ios::cur);
      }
      data.read((char *) &data_array[ncube-1], sizeof (float));
    }
  }
  else
    data.read((char *) data_array, (sizeof (float)) * ncube);
  if (data.fail()) {
    cerr << "ERROR: AvsScalarField::read_binary_field got istream failure\n"
      << "after extracting " << data.gcount() << " bytes. "<< endl;
  }

  put_val_array(data_array); // Install the data that was read.
  big_rigid_free(data_array);
  return data;
}

// Run through the data array and store the values from it via putval.
// data_array is in FORTRAN order.
// Many derived classes will want to override this since putval calls
// are inefficient.
bool AvsScalarField::put_val_array (const float* data_array)
{
  Coord rspace;
  Coord rgrid;
  int h=0;
  bool ok = true;
  for (int k=0; k<grid_dim; ++k) {
    for (int j=0; j<grid_dim; ++j) {
      for (int i=0; i<grid_dim; ++i) {
	Coord rgrid(static_cast<float>(i),
		static_cast<float>(j),
		static_cast<float>(k));
	rspace = (rgrid - grid_center_in_grid)*spacing + grid_center_in_space;
	ok = ok && putval(rspace, data_array[h]);
	++h;
      }
    }
  }
  if (!ok)
    cerr << "ERROR: AvsScalarField::put_val_array: putval failed"
	    << endl;
  return ok;
}

istream&
AvsScalarField::read_uniform_coord(istream& data) const
{
  float buf[6];
  data.read((char *) buf, (sizeof (float))*6);
  if (data.fail()) {
    cerr << "ERROR AvsScalarField::read_uniform_coord:\n"
      << "Input failure during binary read." << endl;
    return data;
  }
  Coord min_coor(buf[0], buf[2], buf[4]);
  Coord max_coor(buf[1], buf[3], buf[5]);
  if (!(min_coor == bottom_corner_in_space
	&& max_coor == top_corner_in_space)) {
    cerr << "ERROR AvsScalarField::read_uniform_coord:\n"
      << "min and/or max extent data does not agree with CubeLatSpec."
	<< endl;
    data.clear(ios::failbit);
  }
  return data;
}

ifstream&
AvsScalarField::read_rectilinear_coord(ifstream& data, int n) const
{
  if (n<1 || n>3) {
    cerr << "ERROR: AvsScalarField::read_rectilinear_coord called with n= "
      << n << ",\nrather than n in range 1..3" << endl;
    data.clear(ios::failbit);
    return data;
  }
  float *coord_array = new float[grid_dim];
  // if coord_stride is non-zero, we are under data parsing mode so use
  // and coord_skip, coord_stride to control reading.
  // Otherwise assume the istream
  // ,data, is positioned to the start of the desired input already.
  if (coord_stride[n-1]) {
    if (coord_stride[n-1] < 0 || coord_skip[n-1] < 0) {
      cerr << "ERROR: AvsScalarField::read_rectilinear_coord:\n"
	<< "either coord_stride = " << coord_stride[n-1] << " or coord_skip = "
	  << coord_skip[n-1] << " is nonsensical" << endl;
      data.clear(ios::failbit);
      return data;
    }
    data.seekg(coord_skip[n-1], ios::beg);
    if (coord_stride[n-1]==1)
      data.read((char *) coord_array, (sizeof (float)) * grid_dim);
    else {
      for (int i=0; i<grid_dim-1; ++i) {
	data.read((char *) &coord_array[i], sizeof (float));
	data.seekg(coord_stride[n-1] * (sizeof (float)), ios::cur);
      }
      data.read((char *) &coord_array[grid_dim-1], sizeof (float));
    }
  }
  else
    data.read((char *) coord_array, (sizeof (float)) * grid_dim);
  if (data.fail()) {
    cerr << "ERROR: AvsScalarField::read_rectilinear_coord got istream failure"
      << "\nafter extracting " << data.gcount() << " bytes. "<< endl;
    return data;
  }

  blab3 << "AvsScalarField::read_rectilinear_coord(data, " << n
    << ")\n successfuly extracted " << data.gcount() << " characters " << endl;
  // Comparison to expected values
  float xpct = 9999.9F;
  switch (n) {
  case 1:
    xpct = bottom_corner_in_space.x;
    break;
  case 2:
    xpct = bottom_corner_in_space.y;
    break;
  case 3:
    xpct = bottom_corner_in_space.z;
    break;
  }
  bool ok = true;
  for (int i=0; i<grid_dim; ++i, xpct += spacing) {
    ok = ok && xpct == coord_array[i];
  }
  if (!ok) {
    cerr << "ERROR: AvsScalarField::read_rectilinear_coord got discrepency\n"
      << "in values for coord " << n << endl;
    data.clear(ios::failbit);
  }
  delete coord_array;
  return data;
}

// AvsScalarField_input.cc ends here
