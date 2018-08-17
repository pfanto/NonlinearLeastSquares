#include "lapackwrapper.hpp"
void diag(const int n, double** A, complex<double> * e, complex<double> ** u) {
	char jobvl = 'N';
	char jobvr = 'V';
	int ncopy = n*1;

	double * Acopy = new double[ncopy*ncopy];

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Acopy[j*n + i] = A[i][j]*1.0;
		}
	}

	int lda = n*1;

	double * wr = new double [ncopy]; 
	double * wi  = new double [ncopy];

	int ldvl = 1;
	double *vl = new double[ldvl*n];

	int ldvr = n*1;
	double * vr = new double [ncopy*ncopy];

	int lwork = 6*n;
	double * work = new double[lwork];
	
	int info = 0;

	dgeev(& jobvl, & jobvr, & ncopy, Acopy, & lda, wr, wi, vl, & ldvl, vr, & ldvr, work, & lwork, & info);


	int j = 0;
	while (j < n) {
		//cout << "j = " << j << endl;
		bool iscomp = false;
		complex<double> z = complex<double>(wr[j],wi[j]);
		if (abs(imag(z)) > 1e-10) iscomp = true;
		e[j] = z;
		if (iscomp) e[j+1] = conj(z);
		for (int i = 0; i < n; i++) {
			if (!iscomp) {
				u[i][j] = complex<double>(vr[j*n+i],0.0);
			}
			else {
				u[i][j] = complex<double>(vr[j*n+i],vr[(j+1)*n+i]);
				u[i][j+1] = complex<double>(vr[j*n+i],-vr[(j+1)*n+i]);
			}

		}
		if (!iscomp) j += 1;
		else j += 2;
	}

	delete [] Acopy;
	delete [] wr;
	delete [] wi;
	delete [] vl;
	delete [] vr;
	delete [] work;

}

void symsolve(const int n,  double * x, double ** A, double * y) {
	if (n==1) {
		x[0] = y[0]/A[0][0];
		return;
	}
	int ncopy = n*1;
	double * Acopy = new double[ncopy*ncopy];
	double * b = new double[ncopy];
	int * ipiv = new int [ncopy];
	for (int i = 0; i < ncopy; i++) {
		b[i] = y[i]*1.0;
		ipiv[i] = 0;
		for (int j = 0; j < ncopy; j++) {
			Acopy[j*n + i] = A[i][j]*1.0;
		}
	}
	int lda = n*1;
	int nrhs = 1;
	int ldb = ncopy*1;
	int lwork = 6*n;
	double * work = new double[lwork];
	int info = 0;

	dsysv("Lower",&ncopy,&nrhs,Acopy,&lda,ipiv,b,&ldb,work,&lwork,&info);

	for (int i = 0; i < n; i++) x[i] = b[i];	

}

void symdiag(const int n, double ** A, double * l, double ** O) {
	int ncopy = n*1;
	double * Acopy = new double [ncopy*ncopy];
	double * w = new double[ncopy];
	for (int i = 0; i < ncopy; i++) {
		for (int j = 0; j < ncopy; j++) {
			Acopy[j*ncopy+i] = A[i][j];
		}
	}
	int lda = n*1;
	int lwork = 6*n;
	double * work = new double[lwork];
	int info = 0;

	dsyev("V","L",&ncopy,Acopy,&lda,w,work,&lwork,&info);

	for (int i = 0; i < n; i++) {
		l[i] = w[i]*1.0;
		for (int j = 0; j < n; j++) {
			O[i][j] = Acopy[j*n+i];
		}
	}

	delete [] Acopy;
	delete [] w;
	delete [] work;
}

void syminv(const int n, double ** A, double ** Ainv) {
	int ncopy = n*1;
	double * Acopy = new double [ncopy*ncopy];
	for (int i = 0; i < ncopy; i++) {
		for (int j = 0; j < ncopy; j++) {
			Acopy[j*ncopy+i] = A[i][j];
		}
	}
	int lda = ncopy*1;
	int * ipiv = new int [ncopy];
	double * work = new double [ncopy];
	int info = 0;
	dsytri("Upper",&ncopy,Acopy,&lda,ipiv,work,&info);
	cout << "info = " << info << endl;
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			Ainv[i][j] = Acopy[j*n+i];
			Ainv[j][i] = Ainv[i][j]*1.0;
		}
	}
	
	delete [] Acopy;
	delete [] ipiv;
	delete [] work;
}
/*
void diag(vector< vector<double> > & A, const int n, vector< complex<double> > & eigvals, vector< vector< complex<double> > > & eigvecs) {
	//cout << "diag" << endl;
	char jobvl = 'N';
	char jobvr = 'V';
	int ncopy = n*1;
	//cout << "n = " << ncopy << endl;
	double * Acopy = new double [ncopy*ncopy];
	for (int i  = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Acopy[j*n+i] = A.at(i).at(j)*1.0;
		}
	}
	/*cout << "Acopy = " << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) cout << setw(12) << Acopy[j*n+i];
		cout << endl;
	}

	int lda = n*1;

	double * wr = new double [ncopy]; 
	double * wi  = new double [ncopy];

	int ldvl = 1;
	double *vl = new double[ldvl*n];

	int ldvr = n*1;
	double * vr = new double [ncopy*ncopy];

	int lwork = 6*n;
	double * work = new double[lwork];
	
	int info = 0;

	dgeev(& jobvl, & jobvr, & ncopy, Acopy, & lda, wr, wi, vl, & ldvl, vr, & ldvr, work, & lwork, & info);

	for (int i = 0; i < n; i++) {
		complex<double> z = complex<double>(wr[i],wi[i]);
		cout << z << endl;
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << setw(12) << vr[j*n+i];
		}
		cout << endl;
	}

	eigvals.resize(n,complex<double>(0.0,0.0));
	for (int i = 0; i < n; i++)	{
		vector<complex<double> > temp (n,complex<double>(0.0,0.0));
		eigvecs.push_back(temp);
	}
	int j = 0;
	while (j < n) {
		//cout << "j = " << j << endl;
		bool iscomp = false;
		complex<double> z = complex<double>(wr[j],wi[j]);
		if (abs(imag(z)) > 1e-10) iscomp = true;
		eigvals.at(j) = z;
		if (iscomp) eigvals.at(j+1) = conj(z);
		for (int i = 0; i < n; i++) {
			if (!iscomp) {
				eigvecs.at(i).at(j) = complex<double>(vr[j*n+i],0.0);
			}
			else {
				eigvecs.at(i).at(j) = complex<double>(vr[j*n+i],vr[(j+1)*n+i]);
				eigvecs.at(i).at(j+1) = complex<double>(vr[j*n+i],-vr[(j+1)*n+i]);
			}

		}
		if (!iscomp) j += 1;
		else j += 2;
	}
	/*cout << "eigvals = " << endl;
	for (int i = 0; i < n; i++) {
		cout << eigvals.at(i) << endl;
	}

	cout << "eigvecs = " << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) cout << setw(20) << eigvecs.at(i).at(j);
		cout << endl;
	}
	

	return;
}*/
/*
void zdiag(const vector<vector<complex<double> > > & A, const int n, vector<complex<double> > & eigvals, vector< vector< complex<double> > > & eigvecs) {
	//cout << "diag" << endl;
	char jobvl = 'N';
	char jobvr = 'V';
	int ncopy = n*1;
	//cout << "n = " << ncopy << endl;
	complex<double> * Acopy = new complex<double> [n*n];
	for (int i  = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Acopy[j*n+i] = A.at(i).at(j)*1.0;
		}
	}
	/*cout << "Acopy = " << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) cout << setw(12) << Acopy[j*n+i];
		cout << endl;
	}
	int lda = n*1;

	complex<double> * w = new complex<double> [n]; // eigenvalues
	complex<double> * vl;
	int ldvl = 1;
	complex<double> * vr  = new complex<double> [n*n];
	int ldvr = n;

	int lwork = 6*n;
	complex<double> * work = new complex<double> [lwork];

	double * rwork = new double [2*n];

	int info = 0;

	zgeev(& jobvl, & jobvr, & ncopy, Acopy, & lda, w, vl, & ldvl, vr, & ldvr, work , & lwork, rwork, & info);

	eigvals.resize(n,complex<double>(0.0,0.0));
	for (int i = 0; i < n; i++)	{
		vector<complex<double> > temp (n,complex<double>(0.0,0.0));
		eigvecs.push_back(temp);
	}
	
	for (int i = 0; i < n; i++) {
		eigvals.at(i) = w[i];
		for (int j = 0; j < n; j++) {
			eigvecs.at(i).at(j) = vr[j*n+i];
		}
	}

	/*

	cout << "eigvals = " << endl;
	for (int i = 0; i < n; i++) {
		cout << eigvals.at(i) << endl;
	}

	cout << "eigvecs = " << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) cout << setw(20) << eigvecs.at(i).at(j);
		cout << endl;
	}
	
}



int main() {

	int n = 2;
	vector<vector<complex<double> > > A (n,vector<complex<double> > (n));
	A.at(0).at(0) = complex<double>(0.0,0.0); A.at(0).at(1) = complex<double>(0.0,-1.0);
	A.at(1).at(0) = complex<double>(0.0,1.0); A.at(1).at(1) = complex<double>(0.0,0.0);
	vector<complex<double > > eigvals;
	vector<vector<complex<double> > > eigvecs;

	cout << "A = " << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) cout << setw(12) << A.at(i).at(j);
		cout << endl;
	}

	zdiag(A,n,eigvals,eigvecs);

	cout << "eigvals = " << endl;
	for (int i = 0; i < n; i++) {
		cout << eigvals.at(i) << endl;
	}

	cout << "eigvecs = " << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) cout << setw(20) << eigvecs.at(i).at(j);
		cout << endl;
	}

}
*/