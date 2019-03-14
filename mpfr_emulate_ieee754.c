#include <stdio.h>

#include "mpfr.h"

// Emulation program of IEEE754 double precision arithmetic
// ref. MPFR 4.0.2 Mamual

int main(void)
{
	mpfr_t xa, xb; int i;
	float fa, fb;
	double a, b;

// 単精度エミュレート

	// 仮数部 24 bits
	mpfr_set_default_prec(24);

	// アンダーフローの境界指数部: 2^(-125-23)
	mpfr_set_emin(-148);

	// オーバーフローの境界指数部: 2^(+126)
	mpfr_set_emax(126);

	// デフォルト精度桁で初期化
	mpfr_init(xa); mpfr_init(xb);

	// 単精度値をセット
	fb = 34.3f; mpfr_set_d(xb, (double)fb, MPFR_RNDN);

	// a = hex(1.1235)×2^(-125)
 	fa = 0x1.1235P-125; mpfr_set_d(xa, (double)fa, MPFR_RNDN);

	// 単精度計算
	fa /= fb;
	printf("float: a / b = %15.7e\n", fa); 

	// MPFRエミュレート倍精度計算
	i = mpfr_div(xa, xa, xb, MPFR_RNDN);
	mpfr_printf("MPFR emulated float: xa / xb = %15.7RNe\n", xa);
	mpfr_printf("                             = %15.6RNa\n", xa);
	i = mpfr_subnormalize (xa, i, MPFR_RNDN); /* new ternary value */
	mpfr_printf("              subnormalized -> %15.7RNe\n", xa);
	mpfr_printf("                             = %15.6RNa\n", xa);

	// 消去
	mpfr_clear(xa); mpfr_clear(xb);


// 倍精度エミュレート

	// 仮数部 53 bits
	mpfr_set_default_prec(53);

	// アンダーフローの境界指数部: 2^(-1021-52)
	mpfr_set_emin(-1073);

	// オーバーフローの境界指数部: 2^(+1024)
	mpfr_set_emax(1024);

	// デフォルト精度桁で初期化
	mpfr_init(xa); mpfr_init(xb);

	// 倍精度値をセット
	b = 34.3; mpfr_set_d(xb, b, MPFR_RNDN);

	// a = hex(1.1235)×2^(-1021)
 	a = 0x1.1235P-1021; mpfr_set_d(xa, a, MPFR_RNDN);

	// 倍精度計算
	a /= b;
	printf("double: a / b = %25.17e\n",a); 

	// MPFRエミュレート倍精度計算
	i = mpfr_div(xa, xa, xb, MPFR_RNDN);
	mpfr_printf("MPFR emulated double: xa / xb  = %25.17RNe\n", xa);
	mpfr_printf("                               = %25.12RNa\n", xa);
	i = mpfr_subnormalize (xa, i, MPFR_RNDN); /* new ternary value */
	mpfr_printf("                subnormalized -> %25.17RNe\n", xa);
	mpfr_printf("                               = %25.12RNa\n", xa);

	// 消去
	mpfr_clear(xa); mpfr_clear(xb);

	return 0;
}
