#include <stdio.h>
#include <stdlib.h>

#include "gmp.h"

/* watch returns of mpz_get_str */

int main(int argc, char *argv[])
{
	int output_base;
	mpz_t input_val;

	if(argc <= 2)
	{
		fprintf(stderr, "[usage] %s input_decimal_integer_string(<=1023 chars) output_base(<=62)\n", argv[0]);
		return EXIT_SUCCESS;
	}

	output_base = atoi(argv[2]);
	if((output_base < -36) || (output_base < 62) && (abs(output_base) < 2))
	{
		fprintf(stderr, "ERROR!: Illegal base = %d!\n", output_base);
		return EXIT_FAILURE;
	}

	// input as decimal
	mpz_init_set_str(input_val, argv[1], 10);

	// output
	printf("input(decimal)   : %s\n", argv[1]);
	printf("output(base = %2d): %s\n", output_base, mpz_get_str(NULL, output_base, input_val));
	printf("check(input == ?): %s\n", mpz_get_str(NULL, 10, input_val));

	// free
	mpz_clear(input_val);

	return EXIT_SUCCESS;
}
