//******************************************************************************
// mpf_relerr.cpp : Print relative errors of basic arithmetic
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

#include "gmpxx.h"

using namespace std;

// relative error
void relerr(mpf_class &c, mpf_class a, mpf_class ad)
{
	if(ad == 0)
		c = abs(a - ad);
	else
		c = abs((a - ad) / ad);
}

int main(void)
{
	unsigned long prec; // precision

	cout << " Input default prec in bits: "; cin >> prec;

	mpf_set_default_prec(prec); // in bits

	mpf_class a, b, c, d;
	mpf_class ad(0, prec * 2), bd(0, prec * 2), cd(0, prec * 2);

	//a = sqrt(2);
	//b = sqrt(3);
	//mpf_sqrt_ui(a.get_mpf_t(), 2UL);
	//mpf_sqrt_ui(b.get_mpf_t(), 3UL);
	mpf_sqrt_ui(a.get_mpf_t(), 5UL);
	mpf_sqrt_ui(b.get_mpf_t(), 7UL);

	mpf_sqrt_ui(ad.get_mpf_t(), 5UL);
	mpf_sqrt_ui(bd.get_mpf_t(), 7UL);

	//gmp_printf("ad = %97.90Fe\n", ad.get_mpf_t());
	relerr(c, a, ad);
	gmp_printf("relerr(a) = %10.3Fe\n", c.get_mpf_t());
	relerr(c, b, bd);
	gmp_printf("relerr(b) = %10.3Fe\n", c.get_mpf_t());

// unworking setprecision(xx) for mpf_class variables !
//	cout << setprecision(50) << a << " + " << b << " = " << a + b << endl;

	c = a + b; // short
	cd = ad + bd; // long
	relerr(c, c, cd);
	gmp_printf("relerr(a + b) = %10.3Fe\n", c.get_mpf_t());

	c = a - b; // short 
	cd = ad - bd; // long
	relerr(c, c, cd);
	gmp_printf("relerr(a - b) = %10.3Fe\n", c.get_mpf_t());

	c = a * b; // short
	cd = ad * bd; // long
	relerr(c, c, cd);
	gmp_printf("relerr(a * b) = %10.3Fe\n", c.get_mpf_t());

	c = a / b; // short
	cd = ad / bd; // long
	relerr(c, c, cd);
	gmp_printf("relerr(a / b) = %10.3Fe\n", c.get_mpf_t());

	return 0;
}
