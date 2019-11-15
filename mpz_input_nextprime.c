//******************************************************************************
// mpz_input_nextprime.c : Print all prime numbers under inputed integer
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
#include "gmp.h"

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

    if(mpz_cmp_ui(a, 1UL) < 0)
    {
        gmp_printf("a = %Zd is less than zero!\n", a);
        exit(EXIT_FAILURE);
    }

    // prime := 2
    mpz_set_ui(prime, 2UL);
    num_prime = 0;
    while(mpz_cmp(a, prime) >= 0)
    {
        gmp_printf("%d: %Zd\n", ++num_prime, prime);
        mpz_nextprime(prime, prime);
    }

    gmp_printf("<= %Zd\n", a);

	mpz_clear(a);
	mpz_clear(c);
	mpz_clear(prime);

	return 0;
}
