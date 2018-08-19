// Polynomial.h
// Simple class Polynomial (with integer coefficients)
#ifndef Polynomial_h
#define Polynomial_h

#include <iostream>
#include <vector>
using std::vector;

class Polynomial;
// These will be friends of Polynomial
ostream &operator<<( ostream &, const Polynomial & );
istream &operator>>( istream &, Polynomial & );
Polynomial operator*(double, const Polynomial & );
Polynomial operator*(int, const Polynomial & );

class Polynomial {
   friend ostream &operator<<( ostream &, const Polynomial & );
   friend istream &operator>>( istream &, Polynomial & );
   friend Polynomial operator*(double, const Polynomial & );
   friend Polynomial operator*(int, const Polynomial & );
public:
   Polynomial( unsigned = 1 );                   // default constructor
   Polynomial( double, double);                   // explicit constructor
   Polynomial(const vector<double>& parg) : p(parg) {}
   Polynomial( const Polynomial & );              // copy constructor
   ~Polynomial();                            // destructor
   int getDegree() const;                 // return degree
   int getSize() const;                 // return degree
   double getCoefficient(const int) const;
   Polynomial &operator=( const Polynomial & ); // assign Polynomials
   Polynomial operator+( const Polynomial & ) const; // add Polynomials
   Polynomial operator-( const Polynomial & ) const; // subtract Polynomials
   Polynomial operator*( const Polynomial & ) const; // multiply Polynomials
   Polynomial operator*( const double ) const; // multiply Polynomial by double
   Polynomial operator*( const int ) const; // multiply Polynomial by int
   Polynomial operator/( const double ) const; // divide Polynomial by double
   const Polynomial &operator+=( const Polynomial & ); // addition assignment
   const Polynomial &operator-=( const Polynomial & );// subtraction assignment
   const Polynomial &operator*=( const Polynomial & );// multiplication assignment
   Polynomial derivative() const;// derivative
   static int getPolynomialCount();          // Return count of 
   // Polynomials instantiated.
  private:
   vector<double> p;
   static int polynomialCount;  // # of Arrays instantiated
};

#endif
