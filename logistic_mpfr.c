// logistic é ëú
// MPFRî≈
#include <stdio.h>
#include "mpfr.h"

#define MAX_NUM 128

int main()
{
	int i;
	unsigned long prec;
//	mpfr_t x[102];
	mpfr_t x[MAX_NUM];
	mpfr_t lx[MAX_NUM];
	mpfr_t relerr;

	printf("prec(bits) = "); scanf("%ld", &prec);
	mpfr_set_default_prec(prec);

	// èâä˙âª
	mpfr_init(relerr);
	for(i = 0; i < MAX_NUM; i++)
	{
		mpfr_init(x[i]);
		mpfr_init2(lx[i], prec * 2);
	}

	// èâä˙íl
	//x[0] = 0.7501;
	mpfr_set_str(x[0], "0.7501", 10, MPFR_RNDN);
	mpfr_set_str(lx[0], "0.7501", 10, MPFR_RNDN);

//	for(i = 0; i <= 100; i++)
	for(i = 0; i <= 101; i++)
	{
		//if((i % 10) == 0)
		{
			printf("%5d, ", i);
			mpfr_printf("%25.17RNe, ", x[i]);
			mpfr_reldiff(relerr, lx[i], x[i], MPFR_RNDN);
			mpfr_printf("%10.3RNe\n", relerr);
		}

		//x[i + 1] = 4 * x[i] * (1 - x[i]);
		mpfr_ui_sub(x[i + 1], 1UL, x[i], MPFR_RNDN);
		mpfr_mul(x[i + 1], x[i + 1], x[i], MPFR_RNDN);
		mpfr_mul_ui(x[i + 1], x[i + 1], 4UL, MPFR_RNDN);

		//x[i + 1] = 4 * x[i] * (1 - x[i]);
		mpfr_ui_sub(lx[i + 1], 1UL, lx[i], MPFR_RNDN);
		mpfr_mul(lx[i + 1], lx[i + 1], lx[i], MPFR_RNDN);
		mpfr_mul_ui(lx[i + 1], lx[i + 1], 4UL, MPFR_RNDN);

	}

	// è¡ãé
	mpfr_clear(relerr);
	for(i = 0; i < MAX_NUM; i++)
	{
		mpfr_clear(x[i]);
		mpfr_clear(lx[i]);
	}

	return 0;
}
