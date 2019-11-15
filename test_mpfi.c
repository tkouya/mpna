//******************************************************************************
// test_mpfi.c : Test program of MPFI
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

//#include "mpfr.h"
#include "mpfi.h"

int main()
{
	unsigned long prec;
	mpfi_t a, b, c;

	printf("prec, a, b, c\n");
	for(prec = 16; prec <= 1024; prec *= 2)
	{
		// Initialize all variables
		mpfi_init2(a, prec);
		mpfi_init2(b, prec);
		mpfi_init2(c, prec);

		// a = sqrt(2), b = sqrt(3), c = sqrt(5)
		mpfi_set_ui(a, 2UL); mpfi_sqrt(a, a);
		mpfi_set_ui(b, 3UL); mpfi_sqrt(b, b);
		mpfi_set_ui(c, 5UL); mpfi_sqrt(c, c);

		// Print intervals
		//printf("a = ");
		mpfi_out_str(stdout, 10, 10, a); printf(", "); //printf("\n");
		//printf("b = ");
		mpfi_out_str(stdout, 10, 10, b); printf(", "); //printf("\n");
		//printf("c = ");
		mpfi_out_str(stdout, 10, 10, c); printf("\n");

		// Clear all variables
		mpfi_clear(a);
		mpfi_clear(b);
		mpfi_clear(c);
	}

	return 0;
}
