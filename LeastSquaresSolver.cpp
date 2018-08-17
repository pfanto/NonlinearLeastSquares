#include "LeastSquaresSolver.hpp"

static const int NSTEP = 500;
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

	double * f = new double [n]; // function evaluation
	double * b = new double [m]; // rhs Hdx = b
	double * da = new double [m];

	bool convergence = false; // convergence criterion

	for (int it = 0; it < NSTEP; it++) {
		// set up S and f
		for (int k = 0; k < m; k++) {
			for (int i = 0; i < n; i++) {
				S[i][k] = Sensitivity(x[i],k,npars,a,(*model));
				f[i] = (*model)(x[i],npars,a);
				//cout << f[i] <<  endl;
			}
		}

		// set up H and b
		for (int k = 0; k < m; k++) {
			b[k] = 0.0;
			for (int i = 0; i < n; i++) {
				b[k] += S[i][k]*(y[i]-f[i]);
			}

			for (int l = 0; l < m; l++) {
				H[k][l] = 0.0;
				for (int i = 0; i < n; i++) {
					H[k][l] += S[i][k]*S[i][l];
				}
			}
		}
		//cout << "b = " << b[0] << endl;
		//cout << "H = " << H[0][0] << endl;

		// solve for da
		symsolve(m,da,H,b);
		PrintPointerArray("da = ",da,npars);
		// convergence criterion
		double absda = 0.0;
		
		// update a
		for (int k = 0; k < m; k++) {
			absda += da[k]*da[k];
			a[k] += da[k];
		}
		absda = sqrt(absda);
		if (absda <= STAT) {
			convergence = true;
			break;
		}
		// continue

	}
	

	// memory deallocation
	for (int i = 0; i < n; i++) delete [] S[i];
	for (int i = 0; i < m; i++) delete [] H[i];
	delete [] f; 
	delete [] b;
	delete [] da;

	if (convergence) {
		cout << "converged!" << endl;
		return;
	}
	cerr << "could not converge!" << endl;
	throw exception();
}









