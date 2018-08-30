#ifndef LSQSOLV_H
#define LSQSOLV_H

#include <iostream>
#include <iomanip>
#include <cstdlib>

using namespace std;

#include "Sensitivity.hpp"
#include "writeout.hpp"
#include "lapackwrapper.hpp"

// Levenberg-Marquardt minimization
void LeastSquaresSolver(const int n, const int m, const int npars, 
	double * x, double * y, double * a, double (*model)(const double, const int, double *),
	const int NSTEP, const double STATIONARY_TOL, const double lambda_initial);

// nonlinear least-squares minimization using inverse Hessian with regularization
void LeastSquaresSolverRegularized (const int n, const int m, const int npars, 
	double * x, double * y, double * a, double (*model)(const double, const int, double *),
	const int NSTEP, const double STATIONARY_TOL, const double mu);

#endif
