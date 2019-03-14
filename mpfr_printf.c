#include <stdio.h>

#include "mpfr.h"

int main()
{
	mpfr_t a;

	mpfr_init2(a, 299); // 299ビット ≒ 10進100桁
	mpfr_const_pi(a, MPFR_RNDN); // a = 3.1415....

	mpfr_printf("(1) a = %Re\n", a); // デフォルト丸め
	mpfr_printf("(1) a = %100.90Re\n", a); // デフォルト丸め
	mpfr_printf("(2) a = %100.90RNe\n", a);// MPFR_RNDN丸め指定(1)
	mpfr_printf("(3) a = %100.90R*e\n", MPFR_RNDN, a); // MPFR_RNDN丸め指定(2)

	mpfr_printf("(b) a    = %RNb\n", a); // 2進出力
	mpfr_mul_ui(a, a, 10UL, MPFR_RNDN); // a *= 10
	mpfr_printf("(b) a*10 = %RNb\n", a); // 2進出力

	mpfr_clear(a); // クリア

	return 0;
}
