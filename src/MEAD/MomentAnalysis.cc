#include "MEAD/MomentAnalysis.h"
#include "MEAD/MomentAnalysis_tmplts.h"

#include <functional>
using std::bind1st;
using std::mem_fun;

Moments::Moments(unsigned maxell)
{
  // Allocate all space needed for moments up to ell=maxell
  // and intiallize to zero
  _momvec.reserve(maxell+1);
  for (unsigned ell = 0; ell <= maxell; ++ell)
    _momvec.push_back(std::vector<momtype> (2*ell+1, momtype(0,0)));
}

Moments Moments::operator-() const
{
  return (*this)*(-1.0);
}

Moments Moments::operator+(const Moments o) const
{
  int mx = ellmax();
  Moments retval(mx);
  for (int ell=0; ell <= mx; ++ell)
    for (int m = -ell; m <= ell; ++m)
      retval(ell,m) = (*this)(ell,m) + o(ell,m);
  return retval;
}

Moments Moments::operator-(const Moments o) const
{
  int mx = ellmax();
  Moments retval(mx);
  for (int ell=0; ell <= mx; ++ell)
    for (int m = -ell; m <= ell; ++m)
      retval(ell,m) = (*this)(ell,m) - o(ell,m);
  return retval;
}

Moments Moments::operator*(momtype f) const
{
  int mx = ellmax();
  Moments retval(mx);
  for (int ell=0; ell <= mx; ++ell)
    for (int m = -ell; m <= ell; ++m)
      retval(ell,m) = (*this)(ell,m)*f;
  return retval;
}

Moments Moments::operator*(double f) const
{
  return (*this) * momtype(f);
}


Moments
momentsOfElstatPot(ElstatPot esp, unsigned maxell,
		   float radius, Coord center)
{
  return momentsOfPotfun(bind1st(mem_fun(&ElstatPot::value), &esp),
			 maxell, radius, center);
}  

inline double trigtrim(double cosorsin)
{
  // Arg is supposed to be cos or sin and therefore in (-1.0 ... 1.0)
  // Return the arg unchanged if this condition is met.
  // Silently trim it if the error's not too bad. Puke if it's bad.
  const double tol = 0.001;
  if (cosorsin > 1.0) {
    if (cosorsin > 1.0 + tol)
      ::error("trigtrim in module MomentAnalysis: cos or sin > 1.0");
    return 1.0;
  }
  else if (cosorsin < -1.0) {
    if (cosorsin < -1.0 - tol)
      ::error("trigtrim in module MomentAnalysis: cos or sin < -1.0");
    return -1.0;
  }
  return cosorsin;
}

#include <iomanip>
using std::setprecision;

Moments
momentsOfChargeDist(ChargeDist::const_iterator chbegin,
		    ChargeDist::const_iterator chend,
		    unsigned maxell, Coord center)
  // Return the spherical moments, up to ell=maxell, for rho.
{
  Moments mom(maxell);
  SphericalHarmonic_series Y(maxell);
  for (ChargeDist::const_iterator ich=chbegin; ich != chend; ++ich) {
    const Coord rvec = (*ich).coord - center;
    const double r = sqrt(rvec*rvec);
    if (r==0) {
      mom(0,0) += double((*ich).charge) * Y(0,0)(0.0, 0.0);
      continue;
    }
    const double costheta = trigtrim(rvec.z/r);
    const double theta = acos(costheta); // wasteful, but Y needs it.
    const double sintheta = trigtrim(costheta*costheta < 1.0 ?
				     sqrt(1.0 - costheta*costheta) : 0.0);
    const double sinphi = trigtrim(sintheta != 0.0 ?
				   rvec.y/r/sintheta : 0.0);
    const double maybe_phi = asin(sinphi);
    const double phi = ( rvec.x > 0.0 ? maybe_phi : M_PI - maybe_phi );

    for (unsigned ell=0; ell <= maxell; ++ell) {
      const double rpow = pow(r, static_cast<int>(ell));
      for (unsigned m=0; m <= ell; ++m) {
	mom(ell,m) += conj(Y(ell,m)(theta, phi)) * rpow
	  * double((*ich).charge);
      }
      for (int m=1; m <= static_cast<int>(ell); ++m) {
	Moments::momtype cmom = conj(mom(ell,m));
	mom(ell,-m) = ( (m % 2 == 0) ? cmom : -cmom );
      }
    }
  }
  return mom;
}  

void print_moments(Moments debmoms)
{
  cerr << "Non-zero moments up to " << debmoms.ellmax() << endl;
  for (int ell=0; ell <= static_cast<int>(debmoms.ellmax()); ++ell) {
    for (int m=-ell; m <= ell; ++m) {
      Moments::momtype mom = debmoms(ell, m);
      if (norm(mom) < 1.0e-15)
	continue;
      cout << ell << "   " << m << "   " << mom << endl;
    }
  }
}

void compare_moments(const Moments& moma, const Moments& momb,
		     string aname, string bname)
{
  const double zero_if_smaller = 1.0e-15;
  const int ellmax = moma.ellmax();
  if (momb.ellmax() != ellmax) {
    cerr << "WARNING attempted comparing moment sets of different sizes"
	 << endl;
    return;
  }
  cout << "Ratios of moments of " << aname << " to " << bname
       << " for non-zero moments up to " << ellmax << endl;
  for (int ell=0; ell <= ellmax; ++ell) {
    Moments::momtype sum(0.0);
    bool suminf = false;
    for (int m=-ell; m <= ell; ++m) {
      const Moments::momtype ma = moma(ell, m);
      const Moments::momtype mb = momb(ell, m);
      const bool azero = norm(ma) < zero_if_smaller;
      const bool bzero = norm(mb) < zero_if_smaller;
      if (azero && bzero) continue;
      cout << ell << "   " << m << "   ";
      if (azero) {
	cout << "0";
	sum += ma; // which is zero, so does nothing really
      }
      else if (bzero) {
	cout << "1/0";
	suminf = true; // so averaging is impossilbe
      }
      else {
	cout << ma/mb;
	sum += ma/mb;
      }
      cout << endl;
    }
    if ( ! suminf ) {
      cout << "For ell = " << ell << " ave. ratio = "
	   << (sum/double(2*ell+1)) << endl;
    }
  }
}
