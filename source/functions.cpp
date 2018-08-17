#include "functions.hpp"

double expdecay(const double x, const int np, double * a) {
	double alpha = a[0];
	return exp(-alpha*x);
}