//******************************************************************************
// logistic_mpfr.c : Logistic function with MPFR
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

#define MAX_NUM 128

int main()
{
	int i;
	unsigned long prec;
//	mpfr_t x[102];
	mpfr_t x[MAX_NUM];
	mpfr_t lx[MAX_NUM];
	mpfr_t relerr;

	printf("prec(bits) = "); scanf("%ld", &prec);
	mpfr_set_default_prec(prec);

	// Initialize mpfr array
	mpfr_init(relerr);
	for(i = 0; i < MAX_NUM; i++)
	{
		mpfr_init(x[i]);
		mpfr_init2(lx[i], prec * 2);
	}

	//x[0] = 0.7501;
	mpfr_set_str(x[0], "0.7501", 10, MPFR_RNDN);
	mpfr_set_str(lx[0], "0.7501", 10, MPFR_RNDN);

//	for(i = 0; i <= 100; i++)
	for(i = 0; i <= 101; i++)
	{
		//if((i % 10) == 0)
		{
			printf("%5d, ", i);
			mpfr_printf("%25.17RNe, ", x[i]);
			mpfr_reldiff(relerr, lx[i], x[i], MPFR_RNDN);
			mpfr_printf("%10.3RNe\n", relerr);
		}

		//x[i + 1] = 4 * x[i] * (1 - x[i]);
		mpfr_ui_sub(x[i + 1], 1UL, x[i], MPFR_RNDN);
		mpfr_mul(x[i + 1], x[i + 1], x[i], MPFR_RNDN);
		mpfr_mul_ui(x[i + 1], x[i + 1], 4UL, MPFR_RNDN);

		//x[i + 1] = 4 * x[i] * (1 - x[i]);
		mpfr_ui_sub(lx[i + 1], 1UL, lx[i], MPFR_RNDN);
		mpfr_mul(lx[i + 1], lx[i + 1], lx[i], MPFR_RNDN);
		mpfr_mul_ui(lx[i + 1], lx[i + 1], 4UL, MPFR_RNDN);

	}

	// delete variables
	mpfr_clear(relerr);
	for(i = 0; i < MAX_NUM; i++)
	{
		mpfr_clear(x[i]);
		mpfr_clear(lx[i]);
	}

	return 0;
}
