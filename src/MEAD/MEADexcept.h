#ifndef _MEADexcept_h
#define _MEADexcept_h

#include <string>

using std::string;

class MEADexcept {
public:
  MEADexcept(	const string& s1 = "",
		const string& s2 = "",
		const string& s3 = "");
  string get_error1();
  string get_error2();
  string get_error3();
private:
  string error1;
  string error2;
  string error3;
};

inline MEADexcept::MEADexcept(	const string& s1,
				const string& s2,
				const string& s3) :
  error1(s1),
  error2(s2),
  error3(s3) {}

inline string MEADexcept::get_error1() {return error1;}
inline string MEADexcept::get_error2() {return error2;}
inline string MEADexcept::get_error3() {return error3;}

#endif //  _MEADexcept_h
