#include "writeout.hpp"

void PrintVector(const string & s, const vector<double> & vec) {
	cout << s << endl;
	int size = vec.size();
	for (int i = 0; i < size; i++) {
		cout << setw(20) << vec.at(i);
		if ((i+1)%5 == 0) cout << endl;
	}
	cout << endl;
}

void PrintComplexVector(const string & s, const vector<complex<double> > & vec) {
	cout << s << endl;
	int size = vec.size();
	for (int i = 0; i < size; i++) {
		cout << setw(20) << vec.at(i);
		if ((i+1)%5 == 0) cout << endl;
	}
	cout << endl;
}

void PrintPointerArray(const string & s, double * x, const int n) {
	cout << s << endl;
	for (int i = 0; i < n; i++) {
		cout << " " << x[i];
		if((i+1) % 5 == 0) cout << endl;
	}
	cout << endl;
}

void PrintComplexPointerArray(const string & s, complex<double> * x, const int n) {
	cout << s << endl;
	for (int i = 0; i < n; i++) {
		cout << setw(12) << x[i];
		if((i+1) % 5 == 0) cout << endl;
	}
	cout << endl;
}