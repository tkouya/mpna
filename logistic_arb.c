// logistic Ê‘œ
// ARB”Å
#include <stdio.h>
#include "arb.h"

int main()
{
	int i;
	slong prec;
	arb_t x[102];

	// ‰Šú‰»
	for(i = 0; i < 102; i++)
		arb_init(x[i]);

	printf("prec(bits) = "); scanf("%ld", &prec);

	// ‰Šú’l
	arb_set_str(x[0], "0.7501", prec);

	for(i = 0; i <= 100; i++)
	{
		if((i % 10) == 0)
		{
			printf("%5d, ", i);
			arb_printd(x[i], 17);
			printf("\n");
		}

		//x[i + 1] = 4 * x[i] * (1 - x[i]);
		arb_sub_ui(x[i + 1], x[i], 1UL, prec);
		arb_neg(x[i + 1], x[i + 1]);
		arb_mul(x[i + 1], x[i + 1], x[i], prec);
		arb_mul_ui(x[i + 1], x[i + 1], 4UL, prec);

	}

	// Á‹
	for(i = 0; i < 102; i++)
		arb_clear(x[i]);

	return 0;
}
