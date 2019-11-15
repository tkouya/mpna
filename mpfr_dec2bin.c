//******************************************************************************
// mpfr_dec2bin.c : Transform inputted decimal string to binary mpfr_t number
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
#include <math.h>

#include "mpfr.h"

int main(int argc, char *argv[])
{
	int i, int_val, digits, reminder, base;
	double val, frac_val;
	mpfr_t mpfr_val;
	char *int_array, *frac_array;

	if(argc <= 1)
	{
		fprintf(stderr, "[usage] %s [dec_format] [base]\n", argv[0]);
		return 0;
	}
	//printf("argv[1] = %s\n", argv[1]);

	mpfr_init2(mpfr_val, 1024);

	base = 2;
	if(argc >= 3)
		base = atoi(argv[2]);

	//val = atof(argv[1]);
	mpfr_set_str(mpfr_val, argv[1], 10, MPFR_RNDN);

	mpfr_printf("%RNf = (", mpfr_val);
	mpfr_out_str(stdout, base, 0, mpfr_val, MPFR_RNDN);
	printf(")_%d\n", base);

	mpfr_clear(mpfr_val);

	return 0;
}

