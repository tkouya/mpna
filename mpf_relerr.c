//******************************************************************************
// mpf_relerr.c : Print relative errors of basic arithmetic
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

/*
[tkouya@cs-muse mpna]$ gcc mpf_relerr.c -lgmp
[tkouya@cs-muse mpna]$ ./a.out
*/
// Relative error
// a ... Approximation
// ad ... Longer approximation using as True value
// void relerr(mpf_t c, mpf_t a, mpf_t ad)
void relerr(mpf_t c, mpf_t a, mpf_t ad)
{
	// c := |(a - ad) / ad|
	mpf_sub(c, a, ad);
	if(mpf_cmp_ui(ad, 0UL) != 0)
		mpf_div(c, c, ad);
	mpf_abs(c, c);
}

int main(void)
{
	unsigned long prec;
	mpf_t a, b, c;
	mpf_t ad, bd, cd; // longer mantissa than a, b and c

	printf("Input default prec in bits: "); while(scanf("%ld", &prec) < 1);

	mpf_set_default_prec(prec); // in bits

	//a = sqrt(2);
	//b = sqrt(3);
//	mpf_init_set_ui(a, 2UL); mpf_sqrt(a, a);
//	mpf_init_set_ui(b, 3UL); mpf_sqrt(b, b);
	mpf_init_set_ui(a, 5UL); mpf_sqrt(a, a);
	mpf_init_set_ui(b, 7UL); mpf_sqrt(b, b);
	mpf_init(c);

	// 256 bits
	mpf_init2(ad, prec * 2); mpf_sqrt_ui(ad, 5UL);
	mpf_init2(bd, prec * 2); mpf_sqrt_ui(bd, 7UL);
	mpf_init2(cd, prec * 2);

//	gmp_printf("ad = %97.90Fe\n", ad);
//	gmp_printf("bd = %97.90Fe\n", bd);
	
	//mpf_sub(c, a, ad); mpf_div(c, c, ad); mpf_abs(c, c);
	relerr(c, a, ad);
	gmp_printf("relerr(a) = %10.3Fe\n", c);
	relerr(c, b, bd);
	gmp_printf("relerr(b) = %10.3Fe\n", c);

	mpf_add(c, a, b); // [short] c := a + b;
	mpf_add(cd, ad, bd); // [long] cd := ad + bd;
	relerr(c, c, cd);
	gmp_printf("relerr(a + b) = %10.3Fe\n", c);

	mpf_sub(c, a, b); // [short] c = a - b;
	mpf_sub(cd, ad, bd); // [long] cd = ad - bd;
	relerr(c, c, cd);
	gmp_printf("relerr(a - b) = %10.3Fe\n", c);

	mpf_mul(c, a, b); // [short] c = a * b;
	mpf_mul(cd, ad, bd); // [long] cd = ad * bd;
	relerr(c, c, cd);
	gmp_printf("relerr(a * b) = %10.3Fe\n", c);

	mpf_div(c, a, b); // [short] c = a / b;
	mpf_div(cd, ad, bd); // [long] c = a / b;
	relerr(c, c, cd);
	gmp_printf("relerr(a / b) = %10.3Fe\n", c);

	mpf_clear(a);
	mpf_clear(b);
	mpf_clear(c);

	mpf_clear(ad);
	mpf_clear(bd);
	mpf_clear(cd);

	return 0;
}
