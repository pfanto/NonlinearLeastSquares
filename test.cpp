#include "test.hpp"

int main(void) {
	cout.setf(ios::scientific,ios::floatfield);
	cout << setprecision(8);
	cout << "TEST" << endl;

	const int n = 100;
	const double xlow = 0.0;
	const double dx = 0.1;
	const double alpha = 0.5;

	double * x = new double [n];
	double * y = new double [n];
	for (int i = 0; i < n; i++) {
		double xval = xlow + i*dx;
		x[i] = xval;
		y[i] = exp(-alpha*xval);
		cout << x[i] << " " << y[i] << endl;
	}
	PrintPointerArray("x = ",x,n);
	PrintPointerArray("y = ",y,n);
	const int npars = 1;
	const int m = 1;
	double * a = new double [npars];
	a[0] = 10.0;

	double (*func) (const double, const int, double * ) = expdecay;

	double s00 = Sensitivity(1.0,0,npars,a,(*func));
	cout << "s00 = " << s00 << endl;

	LeastSquaresSolver(n,m,npars,x,y,a,(*func), 500, 1.0e-5, 0.001);
	PrintPointerArray("a = ",a,npars);



	delete [] x;
	delete [] y;
	delete [] a;

}