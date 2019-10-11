//******************************************************************************
// mpz_mersenne.c : Mersenne Prime Number Calculator by GNU MP
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
// Compile: cc mpz_mersenne.c -lgmp
// 
//******************************************************************************
#include <stdio.h>
#include <stdlib.h>

#include "gmp.h"

/* Mersenne Prime Number */
/* 2^p - 1

/* cf.) http://primes.utm.edu/mersenne/ */

//int max_index = 50; // 2018-01-05
int max_index = 51; // 2019-02-05
unsigned long known_mersenne_prime_p[] = {
	       2,        3,        5,        7,       13,       17,       19,       31,       61,       89, \
	     107,      127,      521,      607,     1279,     2203,     2281,     3217,     4253,     4423, \
	    9689,     9941,    11213,    19937,    21701,    23209,    44497,    86243,   110503,   132049, \
	  216091,   756839,   859433,  1257787,  1398269,  2976221,  3021377,  6972593, 13466917, 20996011, \
	24036583, 25964951, 30402457, 32582657, 37156667, 42643801, 43112609, 57885161, 74207281, 77232917, \
	82589933 \
};

int main(void)
{
	mpz_t a;
	int index;

	mpz_init(a);

	// input
	fprintf(stderr, "2^p[index] - 1 = ? ");
	fprintf(stderr, "index(1 - %d) = ", max_index);
	while(scanf("%d", &index) < 1);
	if((index > max_index) || (index <= 0))
	{
		fprintf(stderr, "ERROR: index = %d is illegal!\n", index);
		return EXIT_FAILURE;
	}

	index--;
	printf("p[%d] = %ld\n", index + 1, known_mersenne_prime_p[index]);

	// a := 2^p - 1
	mpz_ui_pow_ui(a, 2UL, known_mersenne_prime_p[index]);
	//gmp_printf("2^%ld     = %Zd\n", known_mersenne_prime_p[index], a);
	mpz_sub_ui(a, a, 1UL);
	gmp_printf("2^%ld - 1 = %Zd\n", known_mersenne_prime_p[index], a);

	// is prime number ?
	if(mpz_probab_prime_p(a, 15))
		printf(" is a prime number!\n");

	mpz_clear(a);

	return 0;
}
