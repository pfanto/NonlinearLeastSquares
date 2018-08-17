#include "Sensitivity.hpp"

static const double SENSDIFF = 1.0e-4;

double Sensitivity(const double xi, const int k, const int npars, double * a,
	double (*model) (const double, const int, double * )) {

	//cout.setf(ios::scientific,ios::floatfield);
	//cout << setprecision(8);
	//cout << "xi, k, npars = " << xi << " " << k << " " << npars << endl;
	//PrintPointerArray("a = ",a,npars);

	double ak = a[k];
	double da = ak*SENSDIFF;
	//cout << "ak, da = " << ak << " " << da << endl;

	a[k] += da;
	//PrintPointerArray("a plus = ",a,npars);
	double fplus = (*model)(xi,npars,a);
	//cout << "fplus = " << fplus << endl;
	a[k] -= 2.0*da;
	double fminus = (*model)(xi,npars,a);
	a[k] = ak*1.0;
	double df = 0.5*(fplus-fminus)/da;
	return df;

}