//******************************************************************************
// mpz_factorial.c : Factorial calculation with multiple precision integer
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
[tkouya@cs-muse mpna]$ gcc mpz_factorial.c -lgmp
tkouya@cs-athena:~/na/mpna$ ./mpz_factorial
n = 10
10! = 3628800
tkouya@cs-athena:~/na/mpna$ ./mpz_factorial
n = 100
100! = 93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000
tkouya@cs-athena:~/na/mpna$ ./mpz_factorial
n = 1000
1000! = 4023872600...939410970027753472000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
*/

// unsigned long factorial
unsigned long int ul_factorial(unsigned long int n)
{
	unsigned long ret = n;

	if(n <= 1)
		ret = 1;

	while(n-- >= 2)
	{
		ret *= n;
		printf("n, ret = %lu, %lu\n", n, ret);
	}

	return ret;
}

int main(void)
{
	unsigned long int n;
	mpz_t factorial;

	mpz_init(factorial);

	// input
	printf("n = ");
	while(scanf("%ld", &n) < 1);

	// unsigned long int factorial
	printf("[unsigned long] %ld! = %ld\n", n, ul_factorial(n));

	// factorial
	mpz_fac_ui(factorial, n);

	gmp_printf("[GNU MP mpz] %d! = %Zd\n", n, factorial);

	mpz_clear(factorial);

	return 0;
}
