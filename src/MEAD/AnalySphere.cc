#include "MEAD/AnalySphere.h"
#include "MEAD/DielectricSphere.h"
#include "MEAD/ElectrolyteEnvironment.h"
#include "MEAD/ChargeDist.h"
#include "MEAD/PhysCond.h"
#include "MEAD/Angle.h"
#include <math.h>

#define EPSILON 1.0e-16
#define SIZE_INCR 2

const double AnalySphere::epsilon = EPSILON;

AnalySphere::AnalySphere(DielectricSphere* e, ChargeDist_lett* r,
			 ElySphere* ely, int maxterm)
: AnalyticEP(e, r, ely), eps(e), rho(r), electrolyte(ely), l(maxterm)
{
  blab3 << "Entering AnalySphere ctor." << endl;
  size = l + SIZE_INCR;
  epsin = e->epsin_value();
  epsext = e->epsext_value();
  rad_diel = e->radius_value();
  rad_ely = ely->get_radius();
  center = e->get_center();
  ionic_strength = ely->ionic_strength();
  if (ionic_strength) {
    double kappasq = 8 * pi * PhysCond::get_conconv() * ionic_strength
      / (PhysCond::get_kBolt() * PhysCond::get_T() * epsext);
    kappa = float(sqrt(kappasq));
  }
  else
    kappa = 0.0;
  //solved = false;
  solved = true;
  blab3 << "Leaving AnalySphere ctor." << endl;
}

void
AnalySphere::solve() {
  double b = rad_diel;
  double a = rad_ely;
  cout << "l = " << l << endl;
  cout << "kappa = " << kappa << endl;
  cout << "a = " << a << endl;
  cout << "kappa * a = " << kappa * a << endl;
  // Make an array of Legendre Polynomials
  legendre.resize(size);
  Polynomial L0(1,0);
  Polynomial L1(0,1);
  legendre[0] = L0;
  if (size > 1)
    legendre[1] = L1;
  int i;
  for (i = 2; i < size; i++ )
    legendre[ i ] = (L1*legendre[i-1]*(2*i - 1) - legendre[i-2]*(i - 1))
      /double(i) ;

  // Make an array of Kirkwood Polynomials
   //Polynomial kirkwood[size];
   kirkwood = new Polynomial[size + 1];
   //kirkwood.resize(size + 10);
  Polynomial K0(1,0);
  Polynomial K1(1,1);
  Polynomial x(0,1);
  kirkwood[0] = K0;
  if (size > 1)
    kirkwood[1] = K1;
  for (i = 2; i < size; i++ )
    kirkwood[ i ] = kirkwood[i-1] + x*x*kirkwood[i-2]
                                    / double((2*i - 1)*(2*i-3)) ;
  //cout << endl;
  Ka.resize(size + 2);
  //D.resize(size); // this line produces warnings
  //cout << "size = " << size;
  A = new double[size];
  B = new double[size ];
  C = new double[size ];
  D = new double[size ];
  for (int ll = 0; ll<= (l+1); ll++){
    // for (int ll = 0; ll<= l; ll++){
   //cout << endl;
   Ka[ll] = 0;
   for (int i = 0; i<= ll; i++){
      //cout << "in AnalySphere.cc , i = " << i
     //<< " kirkwood[ll].getSize() = " << kirkwood[ll].getSize()  << endl;
     Ka[ll] += kirkwood[ll].getCoefficient(i)*pow(kappa*a,i);
     //cout  << "ll = " << ll << "    i= " << i
     // << "  kirkwood[i].getCoefficient(i)*pow(kappa*a,i) = "
     // << kirkwood[i].getCoefficient(i)*pow(kappa*a,i) << endl;
      //cout  << "ll = " << ll << "    i= " << i
     //  << "  kirkwood[i].getCoefficient(i) = "
     //  << kirkwood[ll].getCoefficient(i) << endl;
   }

  }

}

float
AnalySphere::value(Coord c) const
{
  double pot = 0.0;
  double b = rad_diel;
  double a = rad_ely;
  const double eps_in = epsin;
  const double eps_ex = epsext;
  Coord disp = c - center;
  double r = sqrt(disp*disp);
  for (ChargeDist::const_iterator p = rho->pc_begin(); p!=rho->pc_end(); ++p)
    {
      PointCharge pc = *p;
      if (pc.charge == 0)
	continue;
      Coord disp = pc.coord - center;
      Angle theta(disp,c);
    double cosine = theta.cos;
    double r1 = sqrt(disp*disp);
    double q = pc.charge;

    if (r1 >= rad_diel) {
      ::error("ERROR: AnalySphere::value: ",
	      "Charge outside the dielectric radius not supported.",
	      "  Sorry.");
    }

    /* At first it seems very odd to be assembling these series
       coeffs.  here in value!  But Dave L's shortcut here was to
       avoid the m summations of Kirkwood (and associated Legendres,
       etc.) by regarding each field point (or charge) as a new Z-axis.
    */
    //Compute the series coefficients, A_{0}, B_{0}, C_{0},  D_{0}.
    // But he's not following Kirkwood's notation here! (see below)
    C[0] = q/(eps_ex*b);
    //cout  << "C[0] = "  << C[0] << endl;
    D[0] = q*exp(kappa*a)/(a*eps_ex*(1 + kappa*a));
    //cout  << "D[0] = "  << D[0] << endl;
    B[0] = -kappa*q/(eps_ex*(1 + kappa*a));
    //cout  << "B[0] = "  << B[0] << endl;
    A[0] = B[0] + (q/b)*(1/eps_ex - 1/eps_in); // Kirkwood's B_{00}  (below 11)
    //cout  << "A[0] = "  << A[0] << endl;
    //Compute the series coefficient, D_{l}, l > 0
    //const double exp = exp(-kappa*a);
    int i;
    for (i = 1; i <= l; i++ )
      {
	double AA = (2*i+1)*q*pow(r1/a,i)*exp(kappa*a)/a;
	double BB = (i*eps_in + (i+1)*eps_ex)* Ka[i+1];
	double CC = pow(b/a,2*i+1)*i*(eps_ex - eps_in)*(Ka[i+1] - Ka[i]);
	D[i] = AA/(BB+CC);
      }

     //Compute the coefficient B_{l}, l > 0

    for (i = 1; i <= l; i++ )
      {
	double AA = (i+1)*(2*i+1)*q*pow(r1/a,i)/
	  ((i*eps_in + (i+1)*eps_ex)*i*a);
	double BB = exp(-kappa*a)*(i*Ka[i] - (2*i + 1)*Ka[i+1])/double(i);
	double CC = 1 - pow(b/a,2*i+1)*(i+1)*(eps_ex - eps_in)
	  /(i*eps_in + (i+1)*eps_ex);
	B[i] = (AA + BB*D[i])/CC;
      }

      //Compute the coefficient C_{l}, l > 0
    for (i = 1; i <= l; i++ )
      {
	double AA = (2*i+1)*q*pow(r1/b,i)/b;
	double BB = (eps_ex - eps_in)*i*pow(b/a,i);
	double CC = i*eps_in + (i+1)*eps_ex;
	C[i] = (AA+BB*B[i])/CC;
      }

    //Compute the coefficient A_{l}, l > 0
    for (i = 1; i <= l; i++ )
      A[i] = B[i]*pow(b/a,i) + C[i] - q*pow(r1/b,i)/(eps_in * b);

    // BEGIN compute the value of the Kirkwood Polynomial at x = kappa*r
    double* Kr;
    Kr  = new double[size];
    double x = kappa*r;
    int ll;
    for (ll = 0; ll<= l; ll++){
      Kr[ll] = 0;
      for (int i = 0; i<= ll; i++)
	Kr[ll] += kirkwood[ll].getCoefficient(i)*pow(x,i);
    }
    // END compute the value of the Kirkwood Polynomial at x = kappa*r

    // BEGIN compute the value of the Legendre Polynomial at x = cosine(theta)
    double* Px;
    Px  = new double[size];
    //cout << "cosine =  = " <<  cosine << endl;
    for (ll = 0; ll<= l; ll++){
      Px[ll] = 0;
      for (int i = 0; i<= ll; i++)
	Px[ll] += legendre[ll].getCoefficient(i)*pow(cosine,i);
    }
    // END compute the value of the Legendre Polynomial at x = cosine(theta)
    if (r >= a)
      for (int i = 0; i <= l; i++){
	pot += D[i]* pow(a/r,i + 1)*exp(-x)*Kr[i]*Px[i];
	// so from above, D[n] is Kirkwood's A_{nO} / a^{n+1} ?

	/*cout  << "    i= " << i
	      << "    D[i]= " << D[i]
	      << " D[i]* pow(r/a,-i)*pow(r,-1)*pow(r,-i-1)*exp(-x)*Kr[i]*Px[i]    = "
	      <<  D[i]* pow(r/a,-i)*pow(r,-1)*pow(r,-i-1)*exp(-x)*Kr[i]*Px[i]
	      << endl; */
      }
    else if ((r >= b) && (r<a))
      for (int i = 0; i <= l; i++){
	pot += Px[i]*(B[i]*pow(r/a,i) + C[i]*pow(b/r,i+1));
	// so  C[n] is Kirkwood's C_{n0} / b^{n+1}
	// and B[n] is Kirkwood's G_{n0} a^n

	/*cout  << "    i= " << i
	<< "  C[i]*pow(r,-i-1)) = "
	<< C[i]*pow(r,-i-1) << endl;
	cout  << "    i= " << i
	      << " B[i]*pow(r,i) = "
	      << B[i]*pow(r,i) << endl; */
      }
    else if (r < b){
      for (int i = 0; i <= l; i++)
	pot += Px[i]*A[i]*pow(r/b,i);
      // so A[n] is Kirkwood's B_{n0}*b^n
      Coord disp = c - pc.coord;
      double dist = sqrt(disp*disp);
      if (q != 0)
	{
	  assert (dist > 0);
	  pot += q/(eps_in*dist);
	}
    }
    else
      cout << "Logical flaw in AnalySphere" << endl;
  }
  //cout << "pot at " << c << " = " << pot << endl;
return float(pot);
}

//  field is first solved in a coord system in which a charge lies on the
//  positive z-axis and the point at which the field is measured lies in
// the x-z plane.  The transformation from the physical coordinates to the
// convenient coordinates can be accomplished by two succesive rotations,
// R1 and R2.  After the field is calculated in the convenient (primed)
// coordinate system, we apply R2_inverse and then R1_inverse to obtain
// the true field.
Coord
AnalySphere::field(Coord c) const
{
  Coord fld(0,0,0);
  double b = rad_diel;
  double a = b + rad_ely;
  const double eps_in = epsin;
  const double eps_ex = epsext;
  Coord rdisp = c - center;
  double r = sqrt(rdisp*rdisp);

  for (ChargeDist::const_iterator p = rho->pc_begin(); p!=rho->pc_end(); ++p)
    {
      PointCharge pc = *p;
      if (pc.charge == 0)
	continue;
      Coord r1disp = pc.coord - center;
      Angle theta(r1disp,c);
    double cosine = theta.cos;
    //cout << "cosine = " << cosine << endl;
    double sine = theta.sin;
    //cout << "sine = " << sine << endl;
    double r1 = sqrt(r1disp*r1disp);
    double q = pc.charge;
    C[0] = q/(eps_ex*b);
    D[0] = q*exp(kappa*a)/(a*eps_ex*(1 + kappa*a));
    B[0] = -kappa*q/(eps_ex*(1 + kappa*a));
    A[0] = B[0] + (q/b)*(1/eps_ex - 1/eps_in);
    //cout  << "A[0] = "  << A[0] << endl;
    int i;
    for (i = 1; i <= l; i++ )
      {
	double AA = (2*i+1)*q*pow(r1/a,i)*exp(kappa*a)/a;
	double BB = (i*eps_in + (i+1)*eps_ex)* Ka[i+1];
	double CC = pow(b/a,2*i+1)*i*(eps_ex - eps_in)*(Ka[i+1] - Ka[i]);
	D[i] = AA/(BB+CC);
      }
    for (i = 1; i <= l; i++ )
      {
	double AA = (i+1)*(2*i+1)*q*pow(r1/a,i)/
	  ((i*eps_in + (i+1)*eps_ex)*i*a);
	double BB = exp(-kappa*a)*(i*Ka[i] - (2*i + 1)*Ka[i+1])/double(i);
	double CC = 1 - pow(b/a,2*i+1)*(i+1)*(eps_ex - eps_in)
	  /(i*eps_in + (i+1)*eps_ex);
	B[i] = (AA + BB*D[i])/CC;
	//cout << "B["<<i<<"] = " << B[i] << endl;
      }
    for (i = 1; i <= l; i++ )
      {
	double AA = (2*i+1)*q*pow(r1/b,i)/b;
	double BB = (eps_ex - eps_in)*i*pow(b/a,i);
	double CC = i*eps_in + (i+1)*eps_ex;
	C[i] = (AA+BB*B[i])/CC;
	//cout << "C["<<i<<"] = " << C[i] << endl;
      }
    for (i = 1; i <= l; i++ )
      {
	A[i] = B[i]*pow(b/a,i) + C[i] - q*pow(r1/b,i)/(eps_in * b);
	//cout << "A["<<i<<"] = " << A[i] << endl;
      }
    // BEGIN compute the value of the Kirkwood Polynomial at x = kappa*r
    double* Kr;
    Kr  = new double[size];
    double x = kappa*r;
    int ll;
    for (ll = 0; ll<= l; ll++){
      Kr[ll] = 0;
      for (int i = 0; i<= ll; i++)
	Kr[ll] += kirkwood[ll].getCoefficient(i)*pow(x,i);
    }
    // END compute the value of the Kirkwood Polynomial at x = kappa*r
    //
    // BEGIN compute the derivative of the Kirkwood Polynomial at x = kappa*r
    double* dKr;
    dKr  = new double[size];
    if (x != 0)
      for (int ll = 0; ll<= l; ll++) {
	dKr[ll] = 0;
	for (int i = 0; i<= ll; i++)
	  dKr[ll] += ((2*ll + 1 + x)*kirkwood[ll].getCoefficient(i)*pow(x,i)
		      - (2*ll + 1)*kirkwood[ll+1].getCoefficient(i)*pow(x,i))
	    /x;
	dKr[ll] -= (2*ll + 1)*kirkwood[ll+1].getCoefficient(ll+1)*pow(x,ll+1)
	  /x;
      }
    else
      {
	dKr[0] = 0;
	for (int ll = 1; ll<= l; ll++)
	  dKr[ll] = 1;
      }
    // END compute the derivative of the Kirkwood Polynomial at x = kappa*r
    //
    // BEGIN compute the value of the Legendre Polynomial at x = cosine(theta)
    double* Px;
    Px  = new double[size];
    for (ll = 0; ll<= l; ll++){
      Px[ll] = 0;
      for (int i = 0; i<= ll; i++)
	Px[ll] += legendre[ll].getCoefficient(i)*pow(cosine,i);
    }
    // END compute the value of the Legendre Polynomial at x = cosine(theta)
    //
    // BEGIN compute the value of the derivative of the Legendre Polynomial
    // with respect to theta (chain rule)
    double* dPx;
    dPx  = new double[size];
    dPx[0] = 0;
    for (ll = 0; ll< l; ll++)
      dPx[ll+1] = -sine*((ll+1)*Px[ll] + cosine*dPx[ll]);

    //for (int ll = 0; ll <= l; ll++)
    //cout << "dPx["<<ll<<"] = "<< dPx[ll] << endl;
    // END compute the value of the derivative of the Legendre Polynomial

    // BEGIN Transform to Primed Coordinates
    //R1 Rotate r1disp into z
    Coord z(0,0,1);
    Coord n = cross(z,r1disp);
    if (n*n < epsilon)
      n = z;
    n /= sqrt(n*n);
    double cosine_R1 = r1disp.cos_th();
    double sine_R1 = r1disp.sin_th();
    Coord r1p( r1disp * cosine_R1 + n*(n*r1disp)*(1.0 - cosine_R1)
	       + (cross(r1disp,n))*sine_R1);
    //cout << "r1p = " << r1p << endl;
    //R1 Same rotation on rdisp
    Coord rp( rdisp * cosine_R1 + n*(n*rdisp)*(1.0 - cosine_R1)
	      + (cross(rdisp,n))*sine_R1);

    //R2  Rotate rp into x-z plane
    double cosine_R2 = rp.cos_phi();
    double sine_R2 = rp.sin_phi();
    n = z;
    rp = rp * cosine_R2 + n*(n*rp)*(1 - cosine_R2)
      + (cross(rp,n))*sine_R2;
    double E_rp = 0; // radial component of E-field in primed coord syst
    double E_tp = 0; //  theta component of E-field in primed coord syst
    if (r >= a)
      for (int i = 0; i <= l; i++){
	E_tp -= D[i]* pow(a/r,i + 1)*exp(-x)*Kr[i]*dPx[i]/r;
	E_rp += D[i]* pow(a/r,i + 1)*exp(-x)*Px[i]
	  *((i+1)*Kr[i]/r + kappa*Kr[i] - kappa*dKr[i]);
      }
    else if ((r >= b) && (r<a))
      for (int i = 0; i <= l; i++){
	E_rp += (-i*B[i]*pow(r/a,i) + (i+1)*C[i]*pow(b/r,i+1))*Px[i]/r;
	E_tp -= (B[i]*pow(r/a,i) + C[i]*pow(b/r,i+1))*dPx[i]/r;
	//cout << "(B[i]*pow(r/a,i) + C[i]*pow(b/r,i+1))*dPx[i]/r = "
	//   <<  (B[i]*pow(r/a,i) + C[i]*pow(b/r,i+1))*dPx[i]/r << endl;
	//cout << "B[i]*pow(r/a,i) = "
	//   <<  B[i]*pow(r/a,i)  << endl;
	//cout << " C[i]*pow(b/r,i+1))*dPx[i]/r = "
	//   <<   C[i]*pow(b/r,i+1)*dPx[i]/r << endl;
	//cout << "E_tp  = " << E_tp << endl;
      }
    else if (r < b){
      int i;
      for (i = 1; i <= l; i++)
	if (fabs(A[i]) > epsilon)
	  E_rp -= i*Px[i]*A[i]*pow(r/b,i)/r;
      for (i = 0; i <= l; i++)
	if (fabs(A[i]) > epsilon)
	  E_tp -= dPx[i]*A[i]*pow(r/b,i)/r;
      Coord disp = c - pc.coord;
      double dist = sqrt(disp*disp);
      if (q != 0)
	{
	  assert (dist > 0);
	  E_rp += q*(r - r1*cosine)/(eps_in*dist*dist*dist);
	  E_tp += q*r1*sine/(eps_in*dist*dist*dist);
	}
    }
    else
      cout << "Logical flaw in AnalySphere" << endl;
    //cout << "E_rp  = " << E_rp << endl;
    //cout << "E_tp  = " << E_tp << endl;
    Coord E_p(E_rp * sine +  E_tp * cosine, 0, E_rp * cosine -  E_tp * sine);
    //cout << "E_p = " << E_p << endl;
    // END Transform to primed Coordinates


    // Undo R2
    n = z;
    Coord E( E_p * cosine_R2 + n*(n*E_p)*(1 - cosine_R2)
	     - (cross(E_p,n))*sine_R2);
     //cout << "Undo R2, E_p -> " << E << endl;
    // Undo R1
    n = cross(z,r1disp);
 if (n*n < epsilon)
      n = z;
    n /= sqrt(n*n);
    //cout << "n = " << n << endl;
    //cout << "cosine_R1 = " << cosine_R1 << endl;
     //cout << "sine_R1 = " << sine_R1 << endl;
    E = E * cosine_R1 + n*(n*E)*(1 - cosine_R1) - (cross(E,n))*sine_R1;
    fld += E;

    //cout << "E = " << E << endl;
    }
  cout << "fld at " << c << " = " << fld << endl;
  return fld;
}

Coord
AnalySphere::displacement(Coord c) const
{
  Coord Disp(0,0,0);
  double b = rad_diel;
  const float eps_in = epsin;
  const float eps_ex = epsext;
  Coord rdisp = c - center;
  double r = sqrt(rdisp*rdisp);
  if (r >= b)
    Disp = eps_ex * field(c);
  else if (r < b)
    Disp = eps_in * field(c);
   cout << "D at " << c << " = " << Disp << endl;
  return Disp;
}

/*
double AnalySphere::operator* (const ChargeDist& c) const
{
  return 0.0;
}
*/

// AnalySphere.cc ends here
