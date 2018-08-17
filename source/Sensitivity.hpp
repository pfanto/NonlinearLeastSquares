#ifndef SENS_H
#define SENS_H

#include <iostream>
#include <iomanip>
#include <cstdlib>

using namespace std;

#include "writeout.hpp"

double Sensitivity(const double xi, const int k, const int npars, double * a, double (*model) (const double, const int, double * ));

#endif