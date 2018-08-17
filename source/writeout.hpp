#ifndef WRITEOUT_H
#define WRITEOUT_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <array>
#include <vector>
#include <complex>

using namespace std;


// automatically prints out vectors and double pointers.
void PrintVector(const string & s, const vector< double > & vec);
void PrintComplexVector(const string & s, const vector<complex<double> > & vec);
void PrintPointerArray(const string & s, double * x, const int n);
void PrintComplexPointerArray(const string & s, complex<double> * x, const int n);

#endif