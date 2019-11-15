//******************************************************************************
// mpfr_relerr.cpp : Computation of Relative errors
// Copyright (C) 2019 Tomonori Kouya
// 
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation, either version 3 of the License or any later
// version.
// 
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//******************************************************************************
#include <iostream>
#include <iomanip>

#include "mpreal.h"

using namespace std;
using namespace mpfr;

void relerr(mpreal &c, mpreal a, mpreal ad)
{
	if(ad != 0)
		c = abs((a - ad) / ad);
	else
		c = abs(a - ad);

//	cout << "c = " << c << endl;

	return;
}

int main()
{
	unsigned long prec;

	cout << "Input prec in bits: "; cin >> prec;	

	mpreal::set_default_prec(prec);
	mpreal a, b, c;
	mpreal ad(0, prec * 2, MPFR_RNDN), bd(0, prec * 2), cd(0, prec * 2);

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
