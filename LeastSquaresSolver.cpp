#include "LeastSquaresSolver.hpp"

static const int NSTEP = 100;
static const double STAT = 1.0e-2;
static const double lambda_initial = 0.001;
void LeastSquaresSolver(const int n, const int m, const int npars, 
	double * x, double * y, double * a, double (*model)(const double, const int, double *)) {

	cout.setf(ios::scientific,ios::floatfield);
	cout << setprecision(8);

	cout << "Least Squares Solver" << endl;
	// memory allocation
	double ** S = new double * [n]; // sensitivity matrix
	for (int i = 0; i < n; i++) S[i] = new double [m];

	double ** H = new double * [m];
	for (int i = 0; i < m; i++) H[i] = new double [m];

	double lambda = lambda_initial*1.0;

	for (int it = 0; it < NSTEP; it++) {
		for (int i = 0; i < n; i++) {
			double xi = x[i];
			for (int k = 0; k < m; k++) {
				S[i][k] = Sensitivity(xi,k,npars,a,(*model));
			}
			//PrintPointerArray("",S[i],m);
		}
		// set up Hessian matrix
		for (int k = 0; k < m; k++) {
			for (int l = k; l < m; l++) {
				for (int i = 0; i < n; i++) {
					H[k][l] += S[i][k]*S[i][l];
				}
				H[l][k] = H[k][l]*1.0;
				if (k==l) H[k][l] *= (1.0 + lambda);
			}
		}
		
		return;



	}



	// memory deallocation
	for (int i = 0; i < n; i++) delete [] S;
	for (int i = 0; i < m; i++) delete [] H;
}