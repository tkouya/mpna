//******************************************************************************
// mpq_input_convert.c : Covert cyclic fixed point number to rational number
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
#include <string.h>

#include "gmp.h"

/*
a = 0.123
a = 0.123 - 3 charaters
c = 41/333
  = 0.123123
*/

int main(void)
{
	mpz_t a, b;
	mpq_t c;
	char str_a[1024];
	int len_str_a;

	mpz_init(a);
	mpz_init(b);
	mpq_init(c);

	// Input
	printf("a = 0."); while(!scanf("%s", str_a));
	len_str_a = (int)strlen(str_a);
	printf("a = 0.%s - %d charaters\n", str_a, len_str_a);

	// set numerator
	mpz_set_str(a, str_a, 10);

	// set denominator
	mpz_set_ui(b, 0UL);
	while(1)
	{
		//mpz_add_ui(mpq_denref(c), mpq_denref(c), 9UL);
		mpz_add_ui(b, b, 9UL);
		if(--len_str_a <= 0)
			break;

		//mpz_mul_ui(mpq_denref(c), mpq_denref(c), 10UL);
		mpz_mul_ui(b, b, 10UL);
	}
	mpz_set(mpq_denref(c), b);
	mpz_set(mpq_numref(c), a);

	// to canonicalize
	mpq_canonicalize(c);

	gmp_printf("c = %Qd\n", c);
	gmp_printf("  = %g\n", mpq_get_d(c));

	mpz_clear(a);
	mpz_clear(b);
	mpq_clear(c);

	return 0;
}
