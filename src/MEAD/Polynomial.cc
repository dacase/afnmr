// polynomial.cpp Member function definitions for class
// Polynomial
#include <iostream>
using std::ostream;
using std::istream;
using std::endl;
#include <iomanip>

#include <assert.h>
#include "MEAD/Polynomial.h"

#include <algorithm>

#ifdef WIN32
// Because contrary to Windows documentation, algorithm does not
// provide max and min (WinBug)

template <class T> inline T max(const T& a, const T& b)
{
  return b > a ? b : a;
}

template <class T> inline T min(const T& a, const T& b)
{
  return b < a ? b : a;
}
#else
using std::max;
using std::min;
#endif

// Initialize static data member at file scope
int Polynomial::polynomialCount = 0;   // no objects yet

// Default constructor for class Polynomial (default size 1)
Polynomial::Polynomial( unsigned polynomialSize )
{
  ++polynomialCount;          // count one more object
  p.resize(polynomialSize);
  for (unsigned i = 0; i < polynomialSize; i++ )
    p[ i ] = 0;          // initialize polynomial
}
//  Explicit Constructor for class Polynomial
Polynomial::Polynomial( double c0, double c1 )
{
   p.resize(1);
  ++polynomialCount;          // count one more object
   p[ 0 ] = c0;          // initialize polynomial
   if (c1){
     p.resize(2);
     p[ 1 ] = c1;          // initialize polynomial
   }
}



// Copy constructor for class Polynomial
// must receive a reference to prevent infinite recursion
Polynomial::Polynomial( const Polynomial &init )
{
   ++polynomialCount;          // count one more object
   p.resize(init.p.size());
   for (unsigned i = 0; i < init.p.size(); i++ )
     p[ i ] = init.p[ i ];  // copy init into object
}

// Destructor for class Polynomial
Polynomial::~Polynomial()
{
   --polynomialCount;             // one fewer objects
}
// Get the degree  of the polynomial
int Polynomial::getDegree() const { return p.size() - 1; }
// Get the size  of the polynomial
int Polynomial::getSize() const { return p.size(); }
double Polynomial::getCoefficient(const int i) const
{
  assert (i < getSize());
  return p[i];
}
// take the derivative
Polynomial Polynomial::derivative() const
{
  Polynomial temp(unsigned(p.size() - 1));
  for (unsigned i = 1; i < p.size(); i++ )
    temp.p[ i - 1] = p[i] * i;
  return temp;
}
// Overloaded assignment operator
// const return avoids: ( a1 = a2 ) = a3
Polynomial &Polynomial::operator=( const Polynomial &right )
{
   if ( &right != this ) {  // check for self-assignment
      p.resize(right.p.size());
     for (unsigned i = 0; i < right.p.size(); i++ )
         p[ i ] = right.p[ i ];  // copy polynomial into object
   }

   return *this;   // enables x = y = z;
}
const Polynomial &Polynomial::operator+=( const Polynomial &right )
{
  int maxsize = max(p.size() , right.p.size());
  int minsize = min(p.size() , right.p.size());
  p.resize(maxsize);
  int i;
  for (i = 0; i < minsize; i++ )
    p[i] += right.p[ i ];
  for (i = int(p.size()); i < maxsize; i++ )
    p[i] = right.p[ i ];
  return *this;
}

const Polynomial &Polynomial::operator-=( const Polynomial &right )
{
  int maxsize = max(p.size(),right.p.size());
  int minsize = min(p.size(),right.p.size());
  p.resize(maxsize);
  int i;
  for (i = 0; i < minsize; i++ )
    p[i] -= right.p[ i ];
  for (i = int(p.size()); i < maxsize; i++ )
    p[ i ] = - right.p[ i ];
  return *this;
}

const Polynomial &Polynomial::operator*=( const Polynomial &right )
{
  vector<double> temp(p.size() + right.p.size() -1);
  p.resize(temp.size());
  for (unsigned i = 0; i < p.size(); i++ )
    for (unsigned j = 0; j < right.p.size(); j++ )
      temp[ i+j ] += p[i] * right.p[ j ];
  for (unsigned i = 0; i < temp.size(); i++ )
    p[ i ] = temp[ i ];  // copy polynomial into object
  return *this;
}

Polynomial Polynomial::operator+( const Polynomial &right ) const
{
  int maxsize = max(p.size(),right.p.size());
  int minsize = min(p.size(),right.p.size());
  Polynomial temp(maxsize);
  int i;
  for (i = 0; i < minsize; i++ )
    temp.p[ i ] = p[i] + right.p[ i ];
  for (i = int(p.size()); i < maxsize; i++ )
    temp.p[ i ] = right.p[ i ];
  for (i = int(right.p.size()); i < maxsize; i++ )
    temp.p[ i ] = p[ i ];
  return temp;

}
Polynomial Polynomial::operator-( const Polynomial &right ) const
{
  int maxsize = max(p.size(),right.p.size());
  int minsize = min(p.size(),right.p.size());
  Polynomial temp(maxsize);
  int i;
  for (i = 0; i < minsize; i++ )
    temp.p[ i ] = p[i] - right.p[ i ];
  for (i = int(p.size()); i < maxsize; i++ )
    temp.p[ i ] = - right.p[ i ];
  for (i = int(right.p.size()); i < maxsize; i++ )
    temp.p[ i ] = p[ i ];
  return temp;
}

Polynomial Polynomial::operator*( const Polynomial &right ) const
{
  Polynomial temp(unsigned(p.size() + right.p.size() -1));
  for (unsigned i = 0; i < p.size(); i++ )
      for (unsigned j = 0; j < right.p.size(); j++ )
	temp.p[ i+j ] += p[i] * right.p[ j ];
  return temp;
}
// multiplication by a double
Polynomial Polynomial::operator*( const double a ) const
{
  Polynomial temp(unsigned(p.size()));
  for (unsigned i = 0; i < p.size(); i++ )
    temp.p[ i ] = a*p[i];
  return temp;
}
// division by a double
Polynomial Polynomial::operator/( const double a ) const
{
  assert(a != 0);
  Polynomial temp(unsigned(p.size()));
  for (unsigned i = 0; i < p.size(); i++ )
    temp.p[ i ] = p[i]/a;
  return temp;
}
// multiplication by an integer
Polynomial Polynomial::operator*( const int a ) const
{
  Polynomial temp(unsigned(p.size()));
  for (unsigned i = 0; i < p.size(); i++ )
    temp.p[ i ] = p[i]*a;
  return temp;
}
// Return the number of Polynomial objects instantiated
// static functions cannot be const
int Polynomial::getPolynomialCount() { return polynomialCount; }

// friend Polynomial operator*(double, const Polynomial & );
// Overloaded multiplication operator for class Polynomial;
// allows scalar double * polynomial
Polynomial operator*( double a, const Polynomial & poly)
{
  Polynomial temp(unsigned(poly.p.size()));
  for (unsigned i = 0; i < poly.p.size(); i++ )
    temp.p[i] = a * poly.p[ i ];
  return temp;
}

// Overloaded multiplication operator for class Polynomial;
// allows scalar int * polynomial
Polynomial operator*( int a, const Polynomial & poly)
{
  Polynomial temp(unsigned(poly.p.size()));
  for (unsigned i = 0; i < poly.p.size(); i++ )
    temp.p[i] = a * poly.p[ i ];
  return temp;
}

// Overloaded input operator for class Polynomial;
// inputs values for entire polynomial.
istream &operator>>( istream &input, Polynomial &a )
{
   for (unsigned i = 0; i < a.p.size(); i++ )
      input >> a.p[ i ];

   return input;   // enables cin >> x >> y;
}

// Overloaded output operator for class Polynomial
ostream &operator<<( ostream &output, const Polynomial &a )
{
   unsigned i;

   for ( i = 0; i < a.p.size(); i++ ) {
      output << std::setw( 12 ) << a.p  [ i ];

      if ( ( i + 1 ) % 4 == 0 ) // 4 numbers per row of output
         output << endl;
   }

   if ( i % 4 != 0 )
      output << endl;

   return output;   // enables cout << x << y;
}

// Polynomial.cc ends here
