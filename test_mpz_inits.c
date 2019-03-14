#include <stdio.h>
#include <stdlib.h>

#include "gmp.h"

int main(void)
{
	int i;
	mpz_t val[3];

	mpz_inits(val[0], val[1], val[2], NULL);

	mpz_set_str(val[0],  "123456789876543212345678987654321", 10);
	mpz_set_str(val[1],  "987654321234567898765432123456789", 10);
	mpz_set_str(val[2], "-987654321098765432109876543210981", 10);

	for(i = 0; i < 3; i++)
		gmp_printf("val[%d] = %Zd\n", i, val[i]);

	mpz_clears(val[0], val[1], val[2], NULL);

	return EXIT_SUCCESS;
}
