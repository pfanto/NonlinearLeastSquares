#include "test.hpp"

int main(void) {
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
	a[0] = 1.0;

	double (*func) (const double, const int, double * ) = expdecay;



	delete [] x;
	delete [] y;
	delete [] a;

}