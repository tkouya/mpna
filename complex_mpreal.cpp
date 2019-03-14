#include <iostream>
#include <iomanip>
#include <complex>

#include "mpreal.h"

using namespace std;
using namespace mpfr;

int main(void)
{

	mpreal::set_default_prec(128);
	complex<mpreal> a, b, c;

	// a = sqrt(2) + sqrt(3) * I
	// b = sqrt(3) - sqrt(5) * I

	a = complex<mpreal>(sqrt((mpreal)2), sqrt((mpreal)3));
	b = complex<mpreal>(sqrt((mpreal)3), -sqrt((mpreal)5));

	cout << "a = " << a << endl;
	cout << "b = " << b << endl;

	cout << "a + b = " << a + b << endl;
	cout << "a - b = " << a - b << endl;
	cout << "a * b = " << a * b << endl;
	cout << "a / b = " << a / b << endl;
	cout << "sqrt(a) = " << sqrt(a) << endl;

	return 0;
}