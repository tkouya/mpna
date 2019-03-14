/* MPFI test program */
#include <stdio.h>

#include "mpfr.h"
#include "mpfi.h"

int main()
{
	unsigned long prec;
	mpfi_t a, b, c;

	//prec = 256;

	printf("prec, a, b, c\n");
	for(prec = 16; prec <= 1024; prec *= 2)
	{
		// ‰Šú‰»
		mpfi_init2(a, prec);
		mpfi_init2(b, prec);
		mpfi_init2(c, prec);

		// a = sqrt(2), b = sqrt(3), c = sqrt(5)
		mpfi_set_ui(a, 2UL); mpfi_sqrt(a, a);
		mpfi_set_ui(b, 3UL); mpfi_sqrt(b, b);
		mpfi_set_ui(c, 5UL); mpfi_sqrt(c, c);

		// o—Í
		//printf("a = ");
		 mpfi_out_str(stdout, 10, 10, a); printf(", "); //printf("\n");
		//printf("b = ");
		 mpfi_out_str(stdout, 10, 10, b); printf(", "); //printf("\n");
		//printf("c = ");
		 mpfi_out_str(stdout, 10, 10, c); printf("\n");

		// Á‹Ž
		mpfi_clear(a);
		mpfi_clear(b);
		mpfi_clear(c);
	}

	return 0;
}
