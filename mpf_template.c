//******************************************************************************
// mpf_template.c : Sample code of mpf_t variables
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
#include "gmp.h"

int main(void)
{
	unsigned long prec;
	mpf_t a, b, c; // default precision
	mpf_t ad, bd, cd; // larger precsion of a, b and c

	printf("Input default prec in bits: "); while(scanf("%ld", &prec) < 1);

	// set default precision
	mpf_set_default_prec((mp_bitcnt_t)prec); // in bits

	// print default precision
	printf("default prec in bits: %ld\n", prec);

	mpf_init(a); // initialize with default precision
	mpf_inits(b, c, NULL); // multiple initialization with default precision

	// initialize with twice precsion
	mpf_init2(ad, (mp_bitcnt_t)(prec * 2));
	mpf_init2(bd, (mp_bitcnt_t)(prec * 2));
	mpf_init2(cd, (mp_bitcnt_t)(prec * 2));

	// print precision of ad
	printf("prec(ad) = %ld\n", mpf_get_prec(ad));

	// ---------------------
	// variable calculations
	// ---------------------

	// clear one variable
	mpf_clear(a);

	// clear more than one variables
	mpf_clears(b, c, NULL);       
	mpf_clears(ad, bd, cd, NULL); 

	return 0;
}
