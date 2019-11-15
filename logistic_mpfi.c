//******************************************************************************
// logistic_mpfi.c : logistic map with MPFI
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
#include "mpfi.h"

#define MAX_NUM 128

int main()
{
	int i;
	unsigned long prec;
	mpfi_t x[MAX_NUM];
	mpfr_t relerr;

	printf("prec(bits) = "); while(scanf("%ld", &prec) < 1);
	mpfr_set_default_prec(prec);

	// initialize variables
	mpfr_init(relerr);
	for(i = 0; i < MAX_NUM; i++)
		mpfi_init(x[i]);

	// set a initial interval
	mpfi_set_str(x[0], "0.7501", 10);

	for(i = 0; i <= 100; i++)
	{
		if((i % 10) == 0)
		{
			printf("%5d, ", i);
			mpfi_out_str(stdout, 10, 17, x[i]);
			mpfi_diam(relerr, x[i]);
			mpfr_printf("%10.3RNe\n", relerr);
		}

		//x[i + 1] = 4 * x[i] * (1 - x[i]);
		mpfi_ui_sub(x[i + 1], 1UL, x[i]);
		mpfi_mul(x[i + 1], x[i + 1], x[i]);
		mpfi_mul_ui(x[i + 1], x[i + 1], 4UL);

	}

	// delete variables
	mpfr_clear(relerr);
	for(i = 0; i < MAX_NUM; i++)
		mpfi_clear(x[i]);

	return 0;
}
