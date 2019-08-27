//******************************************************************************
// mpfr_pi_simple.c : PI calculation program with MPFR/GMP
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
// [usage] $ gcc mpfr_pi_simple.c -lmpfr -lgmp -lm
//         $ ./a.out
//           Input decimal digits of pi: 10000
//           Decimal digits = 10000
//           Binary digits  = 33220
//           pi(33220 bits) = 3.141592653589793238462643383279502884197169399375
//******************************************************************************
#include <stdio.h>
#include <math.h>
#include "mpfr.h" // MPFR/GMP

int main()
{
	mpfr_prec_t dec_prec, prec;
	mpfr_t mp_pi;

	printf("Input decimal digits of pi: "); while(scanf("%ld", &dec_prec) < 1);

	prec = (mpfr_prec_t)ceil(dec_prec / log10(2.0)); // get binary digits

	printf("Decimal digits = %ld\n", dec_prec);
	printf("Binary digits  = %ld\n", prec);

	mpfr_init2(mp_pi, prec); // initialize mpfr_t variable

	mpfr_const_pi(mp_pi, MPFR_RNDN); // get the value of pi
	mpfr_printf("pi(%ld bits) = %RNe\n", prec, mp_pi); // output pi

	mpfr_clear(mp_pi); // delete variable

	return 0;
}
