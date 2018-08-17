#include "LeastSquaresSolver.hpp"

void LeastSquaresSolver(const int n, const int m, const int npars, 
	double * x, double * y, double * a, double (*model)(const double, const int, double *)) {

	cout.setf(ios::scientific,ios::floatfield);
	cout << setprecision(8);

	cout << "Least Squares Solver" << endl;

	double ** S = new double * [n]; // sensitivity matrix
	cout << "S = " << endl;
	for (int i = 0; i < n; i++) {
		S[i] = new double [m];
		double xi = x[i];
		for (int k = 0; k < m; k++) {
			S[i][k] = Sensitivity(xi,k,npars,a,(*model));
		}
		PrintPointerArray("",S[i],m);
	}

	// memory deallocation
	for (int i = 0; i < n; i++) delete [] S;
}