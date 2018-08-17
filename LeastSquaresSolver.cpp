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


	double * da = new double [m]; // update of parameters
	double * b = new double [m]; // lhs of da equation = S^t (y-f)
	double * f = new double [n]; // model prediction
	double * dloss  = new double [m];

	bool converged = false;

	// presets for first iteration
	double loss = 0.0; 
	for (int i = 0; i < n; i++) {
		double xi = x[i];
		f[i] = (*model)(xi,npars,a);
		loss += (y[i]-f[i])*(y[i]-f[i]);
		for (int k = 0; k < m; k++) {
			S[i][k] = Sensitivity(xi,k,npars,a,(*model));
		}
	}
	for (int k = 0; k < m; k++) {
		b[k] = 0.0;
		for (int l = k; l < m; l++) {
			H[k][l] = 0.0;
			for (int i = 0; i < n; i++) {
				H[k][l] += S[i][k]*S[i][l];
			}
			H[l][k] = H[k][l]*1.0;
		}
		
		for (int i = 0; i < n; i++) {
			b[k] += S[i][k]*(y[i]-f[i]);
		}
	}


	for (int it = 0; it < NSTEP; it++) {
		symsolve(m,da,H,b);

		// approximate loss update
		/*double loss_next = loss*1.0;
		for (int k = 0; k < m; k++) {
			dloss[k] = 0.0;
			for (int i = 0; i < n; i++) dloss[k] -= 2.0*(y[i]-f[i])*S[i][k];
			loss_next += da[k]*dloss[k];
		}
		*/
		for (int k = 0; k < m; k++) a[k] += da[k];
		double loss_next = 0.0;
		for (int i = 0; i < n; i++) {
			double xi = x[i];
			f[i] = (*model)(xi,npars,a);
			loss_next += (y[i]-f[i])*(y[i]-f[i]);
		}

		if (abs(loss_next-loss) <= STAT) {
			converged = true;
			break;
		}
		
		loss = loss_next*1.0;

		for (int i = 0; i < n; i++) {
			double xi = x[i];
			for (int k = 0; k < m; k++) {
				S[i][k] = Sensitivity(xi,k,npars,a,(*model));
			}
		}

		for (int k = 0; k < m; k++) {
			b[k] = 0.0;
			for (int l = 0; l < m; l++) {
				H[k][l] = 0.0;
				for (int i = 0; i < n; i++) {
					H[k][l] += S[i][k]*S[i][l];
				}
				H[l][k] = H[k][l]*1.0;
			}

			for (int i = 0; i < n; i++) b[k] += S[i][k]*(y[i] - f[i]);
		}
	}

	// memory deallocation
	for (int i = 0; i < n; i++) delete [] S[i];
	for (int i = 0; i < m; i++) delete [] H[i];
	delete [] da;
	delete [] b;
	delete [] f;

	if (converged) {
		cout << "converged!" << endl;
		return;
	}
	else {
		cout << "could not converge!" << endl;
		return;

	}
}









