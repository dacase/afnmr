// This is -*- C++ -*-
#ifndef _Angle_h
#define _Angle_h 1

/* A class for Angles in 3-space
*/

#include "MEAD/Coord.h"
#include <iomanip>
#include <assert.h>
// On msvc the define below pulls in macro defs line M_PI
#define _USE_MATH_DEFINES 1
#include <math.h>
class Angle {
public:
  Angle ();
  Angle (Coord,Coord);
  Angle (const Angle& a);
  ostream& print(ostream& ost) {
    return ost << "(" << v1 << ", " << v2 << ", " << cos <<  ", " << sin<<")";
  }
  Coord v1,v2;
  double cos;
  double sin;
private:
static const double epsilon;
};

inline Angle::Angle ()
{v1.x = 0; v1.y = 0; v1.z = 1;
 v2.x = 0; v2.y = 0; v2.z = 1;
 cos = 1;
 sin = 0;}

//inline Angle::Angle ()
//{v1 = (0,0,1); v2 = (0,0,1); cos = 0;}

inline Angle::Angle (Coord v1i, Coord v2i)
{v1=v1i; v2=v2i;
//cout  << setprecision(16)  << "v1 = " << v1 << endl;
//cout  << setprecision(16)  << "v2 = " << v2 << endl;
//double dot_pdk = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
double dot_pdk = v1 * v2;
 double norm = sqrt(v1*v1) * sqrt(v2*v2);
 //double norm = sqrt(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z)
 //* sqrt( v2.x * v2.x + v2.y * v2.y + v2.z * v2.z);
 //cout  << setprecision(16)  << "dot_pdk = " << dot_pdk << endl;
 //cout  << setprecision(16)  << "norm = " << norm << endl;
 if (norm == 0){
   //   std::cerr << "Warning:  you are attempting to take the angle between v1 and v2, and at least one of them has zero length." << std::endl;
   cos = 1;
   sin = 0;
   //cout   << setprecision(16) << "cos = " << cos << endl;
   //cout   << setprecision(16) << "sin = " << sin << endl;
 }
 else
   {
     cos = dot_pdk/norm;
     sin =  sqrt(1 - cos*cos);
     //  cout  << setprecision(16)  << "cos*cos = " << cos*cos << endl;
     // cout   << setprecision(16) << "1 - cos*cos = " <<1 - cos*cos << endl;
     //cout   << setprecision(16) << " sqrt(1 - cos*cos) = "
     //    <<  sqrt(1 - cos*cos) << endl;
     //cout   << setprecision(16) << "cos = " << cos << endl;
     //cout   << setprecision(16) << "sin = " << sin << endl;
     //cout  << setprecision(16)  << "fabs(sin) = "
     //   << fabs(sin) << endl;

     if (fabs(fabs(cos) - 1) < epsilon)
       {
	 sin = 0;
       }
   }

}

inline Angle::Angle (const Angle& a)
{  v1 = a.v1;  v2 = a.v2;  cos = a.cos; sin = a.sin;}

inline ostream& operator<<(ostream& s, Angle a)
{return a.print(s);}

#endif

// Angle.h ends here
