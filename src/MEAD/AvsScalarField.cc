#include <iostream>
#include <string>
#include "MEAD/Coord.h"
#include "MEAD/CubeLatSpec.h"
#include "MEAD/AvsScalarField.h"
#include "MEAD/Bigmem.h"
#include "MEAD/globals.h"

// Some C standard libraries:
#include <string.h>
#include <ctype.h>
#include <time.h>

void
AvsScalarField::write_uniform_native(const CubeLatSpec& cls,
				     const string& fieldname) const
{
  expand_cls(cls);
  // Open .fld (header info) file based on fieldname
  string headname = fieldname + ".fld";
  ofstream head;
  safeopen(head, headname);

  write_header(cls, head);

  head << "field=uniform"
    << "             # field type (uniform, rectilinear, irregular)" << endl;
  // Done with header
  // Two form-feeds separate header from binary data
  char sep[2];
  sep[0] = sep[1] = 12;
  head.write(sep, 2);
  // Write the binary data
  write_binary_field(head);
  // Write the binary coord corners
  float buf[6];
  buf[0] = bottom_corner_in_space.x;
  buf[1] = top_corner_in_space.x;
  buf[2] = bottom_corner_in_space.y;
  buf[3] = top_corner_in_space.y;
  buf[4] = bottom_corner_in_space.z;
  buf[5] = top_corner_in_space.z;

  head.write((char*)buf, 6 * sizeof (float));
  if (!head.good()) {
    cerr << "ERROR: AvsScalarField::write_uniform_native:\n"
      << "An error occured during writing to " << headname << endl;
  }
}

void
AvsScalarField::write_rectilinear_native (const CubeLatSpec& cls,
					  const string& fieldname) const
{
  expand_cls(cls);
  // Open .fld (header info) file based on fieldname
  string headname = fieldname + ".fld";
  ofstream head;
  safeopen(head, headname);

  write_header(cls, head);

  // Add rectilinear/data-formatting stuff to header
  head << "field=rectilinear"
    << "             # field type (uniform, rectilinear, irregular)" << endl;
  // Done with header
  // Two form-feeds separate header from binary data
  char sep[2];
  sep[0] = sep[1] = 12;
  head.write(sep, 2);

  write_binary_field(head);
  write_rectilinear_coord(head);
}

void
AvsScalarField::write_uniform_parsing(const CubeLatSpec& cls,
				      const string& fieldname) const
{
  expand_cls(cls);
  // Open .fld (header info) file based on fieldname
  string headname = fieldname + ".fld";
  ofstream head;
  safeopen(head, headname);

  write_header(cls, head);

  // Add rectilinear/data-formatting stuff to header
  head << "field=uniform"
    << "             # field type (uniform, rectilinear, irregular)" << endl;
  // The .avsdat (binary) file name is based on fieldname
  string dataname = fieldname + ".avsdat";
  head << "variable 1 file=" << dataname << " filetype=binary skip=0" << endl;
  // Done with header

  // Write the binary data
  ofstream data;
  safeopen(data, dataname);
  write_binary_field(data);
  // No Coord data, I guess.
}

void
AvsScalarField::write_rectilinear_parsing(const CubeLatSpec& cls,
					  const string& fieldname) const
{
  expand_cls(cls);
  // Open .fld (header info) file based on fieldname
  string headname = fieldname + ".fld";
  ofstream head;
  safeopen(head, headname);

  write_header(cls, head);

  // Add rectilinear/data-formatting stuff to header
  head << "field=rectilinear"
    << "             # field type (uniform, rectilinear, irregular)" << endl;
  // The .avsdat (binary) file name is based on fieldname
  string dataname = fieldname + ".avsdat";
  int ncube = grid_dim*grid_dim*grid_dim;
  int binfieldsize = ncube*sizeof(float);
  int bincoordsize = grid_dim*sizeof(float);
  head << "variable 1 file=" << dataname << " filetype=binary skip=0" << endl;
  for (int i=0; i<3; ++i) {
    head << "coord " << i+1 << " file=" << dataname
      << " filetype=binary skip=" << binfieldsize + i*bincoordsize << endl;
  }
  // Done with header

  // Write the binary data
  ofstream data;
  safeopen(data, headname);
  write_binary_field(data);
  write_rectilinear_coord(data);
}

// Perhaps this should be made const via casting away.
void
AvsScalarField::expand_cls(const CubeLatSpec& cls)
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

}

void
AvsScalarField::expand_cls(const CubeLatSpec& cls) const
{
  AvsScalarField* th = (AvsScalarField*) this; // Cast away const-ness!
  th->grid_dim = cls.get_grid_dim();
  th->spacing = cls.get_spacing();
  th->grid_center_in_space = cls.get_center();
  float grlen = (float) (grid_dim - 1);
  float halfgrlen = grlen/2;
  th->grid_center_in_grid = Coord(halfgrlen, halfgrlen, halfgrlen);
  float halfsplen = spacing * halfgrlen;
  th->bottom_corner_in_space = grid_center_in_space
    - Coord(halfsplen, halfsplen, halfsplen);
  th->top_corner_in_space = grid_center_in_space
    + Coord(halfsplen, halfsplen, halfsplen);
}

ostream&
AvsScalarField::write_header(const CubeLatSpec& cls, ostream& head) const
{
  head << "# AVS field file " << endl;
  head << "# creation date: ";
  time_t clock = time(0);
  head << ctime(&clock);
  head.flush();
  head << "# Scalar field for cubic lattice: " << cls << endl;
  head << "ndim=3         # number of dimensions in the field " << endl;
  int i=0;
  for (i=0; i<3; ++i) {
    head << "dim" << i+1 << "=" << cls.get_grid_dim()
      << "            # dimension of axis " << i+1 << endl;
  }
  head << "nspace=3       # number of physical coordinates per point " << endl;
  head << "veclen=1       # number of components per point " << endl;
  head << "data=float     # data type (byte,integer,float,double) " << endl;
  head << "min_ext=" << bottom_corner_in_space.x
    << " " << bottom_corner_in_space.y << " " << bottom_corner_in_space.z
      << endl;
  head << "max_ext=" << top_corner_in_space.x
    << " " << top_corner_in_space.y << " " << top_corner_in_space.z
      << endl;
  if (!head.good())
    cerr << "WARNING: AvsScalarField::write_header:\n "
      << "header file output stream is not in good condition after last write"
	<< endl;
  return head;
}


ostream&
AvsScalarField::write_binary_field(ostream& data) const
{
  int ncube = grid_dim*grid_dim*grid_dim;
  float * data_array = (float*) big_rigid_malloc(ncube*(sizeof (float)));
  if (get_val_array(data_array))
    data.write((char *) data_array, ncube * sizeof (float));
  else
    cerr << "ERROR: AvsScalarField::write_binary_field: call to\n"
      << "get_val_array failed.  Skipping the writing of binary field values."
	<< endl;
  big_rigid_free(data_array);
  return data;
}

bool
AvsScalarField::get_val_array(float* data_array) const
{
  Coord rspace;
  Coord rgrid;
  int h=0;
  for (int k=0; k<grid_dim; ++k) {
    for (int j=0; j<grid_dim; ++j) {
      for (int i=0; i<grid_dim; ++i) {
	Coord rgrid(static_cast<float>(i),
		static_cast<float>(j),
		static_cast<float>(k));
	rspace = (rgrid - grid_center_in_grid)*spacing + grid_center_in_space;
	data_array[h] = value(rspace);
	++h;
      }
    }
  }
  return true;  // No way to detect error from calling value.
}

ostream&
AvsScalarField::write_rectilinear_coord(ostream& data) const
{
  float *xcoor = new float[grid_dim];
  float *ycoor = new float[grid_dim];
  float *zcoor = new float[grid_dim];
  for (int i=0; i<grid_dim; ++i) {
    Coord rgrid(static_cast<float>(i),
		static_cast<float>(i),
		static_cast<float>(i));
    Coord rspace = (rgrid - grid_center_in_grid)*spacing
      + grid_center_in_space;
    xcoor[i] = rspace.x;
    ycoor[i] = rspace.y;
    zcoor[i] = rspace.z;
  }
  data.write((char *) zcoor, grid_dim*sizeof(float));
  data.write((char *) ycoor, grid_dim*sizeof(float));
  data.write((char *) xcoor, grid_dim*sizeof(float));
  delete xcoor;
  delete ycoor;
  delete zcoor;
  return data;
}

// FIXME?  Why is this private func. a no-op?  Is it used anywhere?
ostream&
AvsScalarField::write_uniform_coord(ostream& data) const
{
  return data;
}

// AvsScalarField.cc ends here
