// logistic é ëú
// MPFR & MPFIî≈
#include <stdio.h>
#include "mpfr.h"
#include "mpfi.h"

#define MAX_NUM 128

int main()
{
	int i;
	unsigned long prec;
	mpfi_t x[MAX_NUM];
	mpfr_t relerr;

	printf("prec(bits) = "); scanf("%ld", &prec);
	mpfr_set_default_prec(prec);

	// èâä˙âª
	mpfr_init(relerr);
	for(i = 0; i < MAX_NUM; i++)
		mpfi_init(x[i]);

	// èâä˙íl
	mpfi_set_str(x[0], "0.7501", 10);

	for(i = 0; i <= 100; i++)
	{
		if((i % 10) == 0)
		{
			printf("%5d, ", i);
			mpfi_out_str(stdout, 10, 17, x[i]);
			mpfi_diam(relerr, x[i]);
			mpfr_printf("%10.3RNe\n", relerr);
		}

		//x[i + 1] = 4 * x[i] * (1 - x[i]);
		mpfi_ui_sub(x[i + 1], 1UL, x[i]);
		mpfi_mul(x[i + 1], x[i + 1], x[i]);
		mpfi_mul_ui(x[i + 1], x[i + 1], 4UL);

	}

	// è¡ãé
	mpfr_clear(relerr);
	for(i = 0; i < MAX_NUM; i++)
		mpfi_clear(x[i]);

	return 0;
}
