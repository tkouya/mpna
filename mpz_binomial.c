//******************************************************************************
// mpz_binomial.c : Binomial Calculator by GNU MP
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
// Compile: cc mpz_binomial.c -lgmp
// 
//******************************************************************************
#include <stdio.h>
#include "gmp.h"

#define MIN(i, j) (((i)<(j)) ? (i) : (j))

int main(void)
{
	unsigned long int i, j, n, k;
	mpz_t binomial;

	mpz_init(binomial);

	// input
	printf("n = ");
	while(scanf("%ld", &n) < 1);
	printf("k = ");
	while(scanf("%ld", &k) < 1);

	// binomial = combination(n, k)
	mpz_bin_uiui(binomial, n, k);

	gmp_printf("[GNU MP mpz] %ld_C_%ld = %Zd\n", n, k, binomial);

	// Pascal's triagle
/*	for(i = 1; i <= n; i++)
	{
		for(j = 0; j <= MIN(k, i); j++)
		{
			mpz_bin_uiui(binomial, i, j);
			gmp_printf("%Zd ", binomial);
		}
		printf("\n");
	}
*/
	mpz_clear(binomial);

	return 0;
}
