//******************************************************************************
// mpz_factorization.c : Factorization of multiple precision integer
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
$ gcc mpz_prime_factorization.c -lgmp
$ ./a.out
a = 246857968
a = 246857968
246857968 = 2 * 2 * 2 * 2 * 7 * 73 * 109 * 277
*/

int main(void)
{
	int num_prime;
	mpz_t a, prime, c;

	mpz_init(a);
	mpz_init(c);
	mpz_init(prime);

	// input
	printf("a = ");
	gmp_scanf("%Zd", a);
	gmp_printf("a = %Zd\n", a);

	// prime := 2
	mpz_set_ui(prime, 2UL);
	gmp_printf("%Zd = ", a);
	num_prime = 0;
	while(mpz_cmp(a, prime) >= 0)
	{
		while(1)
		{
			mpz_mod(c, a, prime);
			if(mpz_cmp_ui(c, 0UL) != 0)
				break;

			mpz_div(a, a, prime);
			num_prime++;
			if(num_prime >= 2)
				printf(" * ");

			gmp_printf("%Zd", prime);
		}

		mpz_nextprime(prime, prime);
	}
	printf("\n");

	mpz_clear(a);
	mpz_clear(c);
	mpz_clear(prime);

	return 0;
}
