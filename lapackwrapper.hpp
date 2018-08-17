#ifndef LAPACKWRAPPER_H
#define LAPACKWRAPPER_H

#include <iostream>
#include <vector>
#include <complex>
#include <iomanip>
//#include "writeout.h"
using namespace std;

extern "C" {
	extern void dgeev(char* jobvl, char* jobvr, int* n , double * a, int* lda,
		double * wr, double * wi, double * vl , int * ldvl, double * vr, int * ldvr, double * work, 
		int * lwork, int * info);
	extern void zgeev(char * jobvl,char * jobvr, int * n, complex<double> * a, int * lda, complex<double> * w,
		complex<double> * vl, int * ldvl, complex<double> * vr, int * ldvr, complex<double> * work, int * lwork,
		double * rwork, int * info);

	extern void dsysv(char* uplo, int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, 
		int* ldb, double* work, int* lwork, int* info);

	extern void dsyev(char* jobz, char* uplo, int* n, double* a, int* lda, double * w,
		double* work, int* lwork, int* info);

	extern void dsytri(char* uplo, int* n, double * a, int* lda, int* ipiv, double* work, int* info);
}

//void diag(vector< vector<double> > & A, const int n, vector< complex<double> > & eigvals, vector< vector< complex<double> > > & eigvecs);
//void zdiag(const vector<vector<complex<double> > > & A, const int n, vector<complex<double> > & eigvals, vector< vector< complex<double> > > & eigvecs);

void diag(const int n, double** A, complex<double> * e, complex<double> ** u);
void symsolve(const int n,  double * x, double ** A, double * y);
void symdiag(const int n, double ** A, double * l, double ** O);
void syminv(const int n, double ** A, double ** Ainv);
//int main();

#endif