//******************************************************************************
// mpfr_printf.c : Test program to show how to use mpfr_printf
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
#include <stdio.h>
#include "mpfr.h"

int main()
{
	mpfr_t a;

	mpfr_init2(a, 299); // 299 bits approx decimal 100 digits
	mpfr_const_pi(a, MPFR_RNDN); // a = 3.1415....

	mpfr_printf("(1) a = %Re\n", a); // default rounding mode
	mpfr_printf("(1) a = %100.90Re\n", a); // default rounding mode
	mpfr_printf("(2) a = %100.90RNe\n", a);// MPFR_RNDN (1)
	mpfr_printf("(3) a = %100.90R*e\n", MPFR_RNDN, a); // MPFR_RNDN(2)

	mpfr_printf("(b) a    = %RNb\n", a); // print as binary number
	mpfr_mul_ui(a, a, 10UL, MPFR_RNDN); // a *= 10
	mpfr_printf("(b) a*10 = %RNb\n", a); // print as binary number

	mpfr_clear(a); // remove variable

	return 0;
}
