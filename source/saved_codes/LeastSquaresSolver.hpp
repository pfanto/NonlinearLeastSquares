#ifndef LSQSOLV_H
#define LSQSOLV_H

#include <iostream>
#include <iomanip>
#include <cstdlib>

using namespace std;

#include "Sensitivity.hpp"
#include "writeout.hpp"
#include "lapackwrapper.hpp"

void LeastSquaresSolver(const int n, const int m, const int npars, 
	double * x, double * y, double * a, double (*model)(const double, const int, double *));


#endif