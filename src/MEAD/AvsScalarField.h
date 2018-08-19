/* This is -*- C++ -*-

AvsScalarFeild is an abstract base class that provides its
derivatives with the ability to read and write CubeLatSpecs full of
scalar float data to AVS field files.  The derived class must supply the
in-memory storage space for the data along with concrete definitions
for the pure virtual functions, value and putval, which extract and
deposit data from that storage, respectively.  The format of AVS field
files is documented in the "read field" and "Field Data File Format"
AVS man pages, and the "AVS Data Types" chapter of the "AVS
Developer's Guide".  The read field man page seems to be the most
complete of these.  The goal of the AvsScalarField class is to support
the subset of AVS Field files that are consistent with a CubeLatSpec
of scalar float data.  AvsScalarField is intended mainly as a
"mix-in" class to add binary I/O capability to classes that contain
data fields of some kind.

In what follows, cls is a CubeLatSpec, name is a string, val is a float,
c is a Coord, and [io]datstream is an [io]stream or [io]fstream.

Public member functions:

  read(cls, name): Look for a file, name.cls and read AVS field data
from it.  The general format may be either AVS "native field" in which
the data data follows the header info in name.fld, or AVS
"data-parsing input mode" in which the data are in other files that
are specified name.fld.  The field must be consistent with a
cls of scalar values stored as floats which means, in AVS-ese: ndim=3;
dim1=dim2=dim3=cls.grid_dim; nspace=3; veclen=1; field = uniform or
rectilinear; data=float; min_ext and max_ext values exactly agreeing
with the corners of the cube defined by cls.  In the case of
data-parsing input mode, the variable and coord lines must
filetype=binary.  The cls function argument is used as the "authority"
for what is acceptable from name.fld.  If the AVS file(s) supply
coordinate data it is used purely for checking consistency with cls.

  write_* (cls, name): A set of write functions write the data in AVS
field file.  Can select AVS field type as uniform or rectilinear and
AVS I/O mode as native or data-parsing by selection the functions,
write_uniform_native, write_rectilinear_native, write_uniform_parsing
or write_rectilinear_parsing.  Coordinate limits or (rectilinear)
coordinate values are generated from the cls.  Name is used to
generate filenames for the AVS files: name.fld for the field file
containing header info and in the native case, data; and name.avsdat
to contain data in the data-parsing case.

  write (cls, name): The plain write function executes the default
from the above set, which for the base class is write_uniform_native.
Derived classes can easily change this.

Protected member functions:

  DERIVED CLASS MUST SUPPLY value and putval.

  value(c):  Return the value associated with the point, c.
This function should succeed if it is called for c exactly on a lattice
point.  What it does for off lattice points shouldn't matter for the
purposes of AvsScalarField.

  putval(c, val): Store the value, val, so that it is associated with
the lattice point, c.  Return true if operation succeeds, false
otherwise.  This function should succeed if it is called for c exactly
on a lattice point.  put_val_array is the only member function that
calls putval, so if put_val_array is redefined not to use putval,
putval can be a no-op.

  The following protected member functions are supplied so that
users can get more efficiency.  They extract or
and installing data _en masse_.  The implementations here operate
be calls to value or putval, which may result and many function
calls for large sets.  A derived class may overload these to transfer
the data in a more efficient way, bypassing the use of value or
putval (which can then be no-ops).

  put_val_array(arr): Install values in the array, arr, as the
values of the scalar field.  Return true or false on sucess or
failure, respectively.  The array, arr, is ordered like a FORTRAN 3D
array where the first (X) coordinate varies most rapidly and the last
(Z) is least rapid.  THIS IS THE OPPOSITE OF WHAT C PROGRAMMERS ARE
USED TO (and this convention is due to the AVS format).  The array
must be of dimension grid_dim**3 (where grid_dim is supplied by the
same cls that was given to the read function).  A successful call to
read will always result in (among other things) the allocation of an
array such as arr; reading of binary data from an AVS file
into arr; a call, put_val_array(arr), to install the data; and the
deallocation of the array, arr.  SO put_val_array MUST ARRANGE FOR
INDEPENDENT STORAGE OF THE VALUES IN arr SINCE arr WILL BE
DEALLOCATED.  The implementation provided by the base class makes
a call to putval for every element of arr, so it might be slow.

  get_val_array(arr): Load the array, arr, with values of the scalar
field.  Return true or false on sucess or failure, respectively.  The
array is is described above for put_val_array.  A successful call to
one of the public write functions will always result in (among other
things) the allocation of an array such as arr; a call,
get_val_array(arr), to load the data into arr; writing of arr to an
AVS file; and the deallocation of the array, arr.  The implementation
provided by the base class makes a call to value for every element of
arr, so it might be slow.

NOTES FOR DERIVED CLASS WRITERS:  Although Calling a public function does
result in changes to some internally stored values
of an AvsScalarField these change do not affect any subsequent calls
to PUBLIC functions.  One can say the class is "publicly stateless".
However, the public functions call protected functions which ARE
influenced by these changes.  Specifically, get_val_array and
put_val_array are affected by the values of grid_dim, spacing,
and grid_center_in_grid, which are set by the read and write functions
on the basis of their cls arguments.  I think that this will generally
be harmless.  Another possible "gotcha" for derived class writers
trying to override {get,put}_val_array is that grid_dim, et al. are
private to the base class.  The derived class will need to maintain
its own copy of the CubeLatSpec or equivalent.  I suspect most classes
that will want to use AvsScalarField as a "mix-in" class will already
do this anyway.

SHORTCOMINGS: The read function handles multiple input formats,
but the write function produces only one format.
Xdr, as in xdr_float, which would be nice for data
portability, is not supported (yet).  Filetype=ascii is not supported
(yet).  Filetype=unformatted (a FORTRAN thing) is not supported.  If
this were a template it could be made to support data=double, etc.
(FIXME?).  For the uniform field type in data-parsing mode, the AVS
documents do not seem to provide for any representation of coordinate
limits in the binary data, while the min_ext and and max_ext header
entries are optional.  So it can happen that there is no way to
check agreement between the AVS data and the CubeLatSpec.  I guess
this is an AVS shortcoming.
*/
#ifndef _AvsScalarField_h
#define _AvsScalarField_h 1

#include <string>
using std::string;
#include <iostream>

class CubeLatSpec;

#include "MEAD/Coord.h"

class AvsScalarField {
public:
  virtual void write(const CubeLatSpec& cls, const string& fieldname) const
    {write_uniform_native(cls, fieldname);}
  virtual void write_uniform_native(const CubeLatSpec&,
				    const string&) const ;
  virtual void write_rectilinear_native (const CubeLatSpec&,
					 const string&) const ;
  virtual void write_uniform_parsing(const CubeLatSpec&,
				     const string&) const ;
  virtual void write_rectilinear_parsing(const CubeLatSpec&,
					 const string&) const ;
  virtual bool read(const CubeLatSpec&, const string&);
protected:
  virtual float value (Coord) const = 0;
  virtual bool putval (Coord, float) = 0;
// Non-private to allow derived classes to use efficient methods
  virtual bool put_val_array (const float*); // Arg is a FORTRAN-ordered array!
  virtual bool get_val_array (float*) const ; // ditto

private:
  virtual ostream& write_header(const CubeLatSpec&, ostream&) const ;
  virtual ostream& write_binary_field(ostream&) const ;
  virtual ostream& write_rectilinear_coord(ostream&) const ;

  virtual ostream& write_uniform_coord(ostream&) const ;

  virtual istream& read_header(istream&);
  virtual ifstream& read_binary_field(ifstream&);
  virtual istream& read_uniform_coord(istream&) const ;
  virtual ifstream& read_rectilinear_coord(ifstream&, int n) const ;

  virtual void expand_cls(const CubeLatSpec&);
  virtual void expand_cls(const CubeLatSpec&) const ;

  int grid_dim;
  float spacing;
  Coord grid_center_in_space;
  Coord grid_center_in_grid;
  Coord bottom_corner_in_space;
  Coord top_corner_in_space;
  bool separator_seen;
  enum {undefined, uniform, rectilinear} fieldtype;
  Coord min_ext; // Extents needed for uniform fields
  bool min_ext_defined;
  Coord max_ext;
  bool max_ext_defined;
  string var_filespec;
  int var_skip, var_stride;
  string coord_filespec[3];
  int coord_skip[3], coord_stride[3];
};

#endif

// AvsScalarField.h ends here
