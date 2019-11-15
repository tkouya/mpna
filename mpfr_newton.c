//******************************************************************************
// mpfr_newton.c : Newton method with MPFR
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
#include <stdlib.h>

#include "mpfr.h"

// func(x) = 0

// func(x) 
void mpfr_func(mpfr_t ret, mpfr_t x, mpfr_rnd_t rmode)
{
	// ret := x - cos(x)
	//mpfr_cos(ret, x, rmode); mpfr_sub(ret, x, ret, rmode);

	// ret := 1 - a * x
	mpfr_mul_ui(ret, x, 7UL, rmode); mpfr_ui_sub(ret, 1UL, ret, rmode);

	// ret := x^2 - a
	//mpfr_sqr(ret, x, rmode); mpfr_sub_ui(ret, ret, 2, rmode);
}

// dfunc(x) = func(x)' = 
void mpfr_dfunc(mpfr_t ret, mpfr_t x, mpfr_rnd_t rmode)
{
	// ret := 1 + sin(x)
	//mpfr_sin(ret, x, rmode); mpfr_add_ui(ret, ret, 1UL, rmode);

	// ret := -a
	mpfr_set_ui(ret, 7UL, rmode); mpfr_neg(ret, ret, rmode);

	// ret := 2 * x
	//mpfr_mul_ui(ret, x, 2UL, rmode);
}

int main(int argc, char *argv[])
{
	mp_bitcnt_t old_prec, new_prec, max_prec;
	mpfr_t old_x, new_x, func_x, dfunc_x, relerr;
	int itimes;

	if(argc <= 1)
	{
		printf("Usage: %s [max_prec in bits] \n", argv[0]);
		return 0;
	}

	max_prec = (mp_bitcnt_t)atoi(argv[1]);
	if(max_prec <= 2)
	{
		fprintf(stderr, "Error: max_prec = %ld is too small!\n", max_prec);
		return -1;
	}

	mpfr_set_default_prec(max_prec);

	mpfr_init(old_x);
	mpfr_init(new_x);
	mpfr_init(func_x);
	mpfr_init(dfunc_x);
	mpfr_init(relerr);

	//mpfr_set_ui(old_x, 1UL, MPFR_RNDN);
	mpfr_set_d(old_x, 7.0, MPFR_RNDN);
	//mpfr_set_d(old_x, 0.3, MPFR_RNDN);

	// newton iteration
	for(itimes = 1; itimes <= 20; itimes++)
	{
		// new_x := old_x - func(old_x) / dfunc(old_x)
		mpfr_func(func_x, old_x, MPFR_RNDN); mpfr_dfunc(dfunc_x, old_x, MPFR_RNDN); mpfr_div(new_x, func_x, dfunc_x, MPFR_RNDN);
		mpfr_sub(new_x, old_x, new_x, MPFR_RNDN);
		
		// new_x := x * (2 - a * x)
		//mpfr_mul_ui(new_x, old_x, 3UL, MPFR_RNDN); mpfr_ui_sub(new_x, 2UL, new_x, MPFR_RNDN); mpfr_mul(new_x, old_x, new_x, MPFR_RNDN);

		// detect precision
		
		mpfr_reldiff(relerr, new_x, old_x, MPFR_RNDN);

		mpfr_printf("%3d %10.3RNe %50.43RNe\n", itimes, relerr, new_x);

		// old_x := new_x
		mpfr_set(old_x, new_x, MPFR_RNDN);
	}

	mpfr_clear(old_x);
	mpfr_clear(new_x);
	mpfr_clear(func_x);
	mpfr_clear(dfunc_x);
	mpfr_clear(relerr);

	return 0;
}
