/* ARB test program */
#include <stdio.h>
#include "arb.h"

int main()
{
	slong prec;
	arb_t a, b, c;

	//prec = 256;

	// ‰Šú‰»
	arb_init(a);
	arb_init(b);
	arb_init(c);

	printf("prec, a, b, c\n");
	for(prec = 16; prec <= 1024; prec *= 2)
	{

		// a = sqrt(2), b = sqrt(3), c = sqrt(5)
		arb_sqrt_ui(a, (ulong)2, prec);
		arb_sqrt_ui(b, (ulong)3, prec);
		arb_sqrt_ui(c, (ulong)5, prec);

		// o—Í
		//printf("a = ");
		 arb_printd(a, 10); printf(", "); //printf("\n");
		//printf("b = ");
		 arb_printd(b, 10); printf(", "); //printf("\n");
		//printf("c = ");
		 arb_printd(c, 10); printf("\n");
	}
	
	// Á‹
	arb_clear(a);
	arb_clear(b);
	arb_clear(c);

	return 0;
}
