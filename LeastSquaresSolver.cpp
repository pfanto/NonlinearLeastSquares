#include "LeastSquaresSolver.hpp"

static const int NSTEP = 100;
static const double STAT = 1.0e-5;
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

	// memory deallocation
	for (int i = 0; i < n; i++) delete [] S[i];
	for (int i = 0; i < m; i++) delete [] H[i];

	
}









