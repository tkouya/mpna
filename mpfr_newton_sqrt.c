//******************************************************************************
// mpfr_newton.c : Newton method to compute square root with MPFR
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

int main(int argc, char *argv[])
{
	mp_bitcnt_t max_prec;
	mpfr_t old_x, new_x, func_x, const_a, relerr;
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
	mpfr_init(const_a);
	mpfr_init(relerr);

	// a := 3
	mpfr_set_ui(const_a, 3UL, MPFR_RNDN);

	mpfr_set(old_x, const_a, MPFR_RNDN);

	// newton iteration
	for(itimes = 1; itimes <= 20; itimes++)
	{
		// new_x := 1/2 * (old_x +  a / old_x)
		mpfr_div(new_x, const_a, old_x, MPFR_RNDN);
		mpfr_add(new_x, old_x, new_x, MPFR_RNDN);
		mpfr_div_ui(new_x, new_x, 2UL, MPFR_RNDN);		

		mpfr_reldiff(relerr, new_x, old_x, MPFR_RNDN);
		mpfr_printf("%3d %10.3RNe %50.43RNe\n", itimes, relerr, new_x);

		// old_x := new_x
		mpfr_set(old_x, new_x, MPFR_RNDN);
	}

	mpfr_clear(old_x);
	mpfr_clear(new_x);
	mpfr_clear(const_a);
	mpfr_clear(relerr);

	return 0;
}
