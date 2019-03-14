#include <iostream>
#include <iomanip>
#include <complex>

#include "qd/qd_real.h"

using namespace std;

int main(void)
{
	complex<qd_real> a, b, c;

	// a = sqrt(2) + sqrt(3) * I
	// b = sqrt(3) - sqrt(5) * I

	a = complex<qd_real>(sqrt((qd_real)2), sqrt((qd_real)3));
	b = complex<qd_real>(sqrt((qd_real)3), -sqrt((qd_real)5));

	cout << "a = " << a << endl;
	cout << "b = " << b << endl;

	cout << "a + b = " << a + b << endl;
	cout << "a - b = " << a - b << endl;
	cout << "a * b = " << a * b << endl;
	cout << "a / b = " << a / b << endl;
	cout << "sqrt(a) = " << sqrt(a) << endl;

	return 0;
}
