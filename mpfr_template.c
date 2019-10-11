//******************************************************************************
// mpfr_template.c : Sample code of mpfr_t variables
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

int main(void)
{
	unsigned long prec; // precision in bits
	mpfr_t a, b, c; // default precision
	mpfr_t ad, bd, cd; // longer precision tha default

	printf("Input default prec in bits: "); while(scanf("%ld", &prec) < 1);

	// set default precision of MPFR
	mpfr_set_default_prec((mp_prec_t)prec); // in bits

	// print default precision
	printf("default prec in bits: %ld\n", mpfr_get_default_prec());

	mpfr_init(a); // initialize variable with default precision
	mpfr_inits(b, c, NULL); // initialize variables with default precision

	// initialize variables with twice default precision
	mpfr_init2(ad, (mp_prec_t)(prec * 2));
	mpfr_init2(bd, (mp_prec_t)(prec * 2));
	mpfr_init2(cd, (mp_prec_t)(prec * 2));

	// print each precision of variable
	printf("prec(ad) = %ld\n", mpfr_get_prec(ad));

	// --------------------
	//  computations 
	// --------------------

	// delete variable
	mpfr_clear(a); // single

	// delete variables simultaneously
	mpfr_clears(b, c, NULL);       
	mpfr_clears(ad, bd, cd, NULL); 

	return 0;
}
