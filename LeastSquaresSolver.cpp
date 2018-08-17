#include "LeastSquaresSolver.hpp"


void LeastSquaresSolver(const int n, const int m, const int npars, 
	double * x, double * y, double * a, double (*model)(const double, const int, double *),
	const int NSTEP, const double STATIONARY_TOL, const double lambda_initial) {

	cout.setf(ios::scientific,ios::floatfield);
	cout << setprecision(8);
	//cout << NSTEP << " " << STATIONARY_TOL << " " << lambda_initial << endl;
	//cout << "Least Squares Solver" << endl;
	// memory allocation
	double ** S = new double * [n]; // sensitivity matrix
	for (int i = 0; i < n; i++) S[i] = new double [m];


	double ** H = new double * [m];
	for (int i = 0; i < m; i++) H[i] = new double [m];

	double * f = new double [n]; // function evaluation
	double * b = new double [m]; // rhs Hdx = b
	double * da = new double [m];

	bool convergence = false; // convergence criterion
	// presets

	double loss = 0.0;
	// set up S and f

	for (int i = 0; i < n; i++) {
		f[i] = (*model)(x[i],npars,a);
		loss += (y[i]-f[i])*(y[i]-f[i]);
		for (int k = 0; k < m; k++) S[i][k] = Sensitivity(x[i],k,npars,a,(*model));
		//cout << f[i] <<  endl;
	}
	//cout << "loss = " << loss << endl;

	double lambda = lambda_initial*1.0;

	for (int it = 0; it < NSTEP; it++) {
		//cout << "lambda = " << lambda << endl;
		// set up H and b
		for (int k = 0; k < m; k++) {
			b[k] = 0.0;
			for (int i = 0; i < n; i++) {
				b[k] += S[i][k]*(y[i]-f[i]);
			}

			for (int l = k; l < m; l++) {
				H[k][l] = 0.0;
				for (int i = 0; i < n; i++) {
					H[k][l] += S[i][k]*S[i][l];
				}
				H[l][k] = H[k][l]*1.0;
			}
			H[k][k] = H[k][k] + lambda*H[k][k];
		}
		//cout << "b = " << b[0] << endl;
		//cout << "H = " << H[0][0] << endl;

		// solve for da
		symsolve(m,da,H,b);
		//PrintPointerArray("da = ",da,npars);
		// convergence criterion
		
		// update a
		for (int k = 0; k < m; k++) {
			a[k] += da[k];
		}
		//PrintPointerArray("updated a = ",a,npars);
		/*
		absda = sqrt(absda);
		if (absda <= STAT) {
			convergence = true;
			break;
		}
		*/
		
		// calculate new loss
		// update S and f as well
		double loss_new = 0.0;

		for (int i = 0; i < n; i++) {
			f[i] = (*model)(x[i],npars,a);
			loss_new += (y[i]-f[i])*(y[i]-f[i]);
			for (int k = 0; k < m; k++) S[i][k] = Sensitivity(x[i],k,npars,a,(*model));
		}
		//cout << "loss_new = " << loss_new << endl;

		// convergence criterion

		if(abs(loss_new-loss) <= STATIONARY_TOL) {
			convergence = true;
			break;
		}
		
		if (loss_new >= loss) {
			lambda *= 10.0;
			for (int k = 0; k < m; k++) a[k] -= da[k];
			// re-update f[i] and S
			for (int i = 0; i < n; i++) {
				f[i] = (*model)(x[i],npars,a);
				for (int k = 0; k < m; k++) S[i][k] = Sensitivity(x[i],k,npars,a,(*model));
			}

		}
		else {
			lambda /= 10.0;
		}
		

		// update loss
		loss = loss_new*1.0;

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









