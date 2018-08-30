#ifndef SENS_H
#define SENS_H

#include <iostream>
#include <iomanip>
#include <cstdlib>

using namespace std;

#include "writeout.hpp"

// derivative of model function with respect to parameters

double Sensitivity(const double xi, const int k, const int npars, double * a, double (*model) (const double, const int, double * ));

#endif
