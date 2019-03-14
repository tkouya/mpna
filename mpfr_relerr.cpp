#include <iostream>
#include <iomanip>

#include "mpreal.h"

using namespace std;
using namespace mpfr;

void relerr(mpreal &c, mpreal a, mpreal ad)
{
	c = abs((a - ad) / ad);

	cout << "c = " << c << endl;

	return;
}

int main()
{
	unsigned long dprec;

	cout << "Input decimal prec: "; cin >> dprec;	

	mpreal::set_default_prec(mpfr::digits2bits(dprec));
	mpreal a, b, c;
	mpreal ad(0, digits2bits(dprec * 2), MPFR_RNDN), bd(0, digits2bits(dprec * 2)), cd(0, digits2bits(dprec * 2));

	a = sqrt((mpreal)2UL);
	b = sqrt((mpreal)3UL);

	cout << "a.prec = " << a.get_prec() << endl;
	cout << "ad.prec = " << ad.get_prec() << endl;
	//ad = sqrt((mpreal)2UL);
	//bd = sqrt((mpreal)3UL);
	mpfr_sqrt_ui(ad.mpfr_ptr(), 2UL, MPFR_RNDN);
	mpfr_sqrt_ui(bd.mpfr_ptr(), 3UL, MPFR_RNDN);

	relerr(c, a, ad);
	cout << setprecision(3) << "relerr(a) = " << c << endl;
	relerr(c, b, bd);
	cout << setprecision(3) << "relerr(b) = " << c << endl;

	// working setprecision for mpreal variables very well !
	relerr(c, a + b, ad + bd);
	cout << setprecision(3) << "relerr(a + b) = " << c << endl;
	relerr(c, a - b, ad - bd);
	cout << setprecision(3) << "relerr(a - b) = " << c << endl;
	relerr(c, a * b, ad * bd);
	cout << setprecision(3) << "relerr(a * b) = " << c << endl;
	relerr(c, a / b, ad / bd);
	cout << setprecision(3) << "relerr(a / b) = " << c << endl;


	return 0;
}
