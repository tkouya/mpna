// mpf_template.c: mpf_t型の使い方

#include <stdio.h>
#include <math.h>
#include "gmp.h"

int main(void)
{
	unsigned long prec, dec_prec;
	char format[128] = "";
	mpf_t a, b, c; // デフォルトの精度
	mpf_t ad, bd, cd; // a, b, cより長い桁

	printf("Input default prec in bits: "); scanf("%ld", &prec);

	// デフォルトの仮数部ビット数をセット
	mpf_set_default_prec((mp_bitcnt_t)prec); // in bits

	// デフォルトの精度を表示
	printf("default prec in bits: %ld\n", prec);

	mpf_init(a); // デフォルト精度
	mpf_inits(b, c, NULL); // デフォルト精度でまとめて初期化

	// デフォルト仮数部の2倍に設定
	mpf_init2(ad, (mp_bitcnt_t)(prec * 2));
	mpf_init2(bd, (mp_bitcnt_t)(prec * 2));
	mpf_init2(cd, (mp_bitcnt_t)(prec * 2));

	// 変数ごとの仮数部ビット数を表示
	printf("prec(ad) = %ld\n", mpf_get_prec(ad));

	// 値をセットして表示
	mpf_set_ui(a, 5UL);
	mpf_set_d(b, sqrt(2.0));
	mpf_set_str(c, "3.1415926535897932384626433832795e+1568759", 10);
	gmp_printf("a = %Fe, b = %Fe\n", a, b);
	gmp_printf("c = %50.43Fe\n", c);

	// 標準型に変換して表示
	printf("a = %ld\n", mpf_get_ui(a));
	printf("b = %f\n", mpf_get_d(b));
	printf("c = %25.17e\n", mpf_get_d(c));

	// フル表示
	dec_prec = (unsigned long)ceil(log10(2.0) * prec);
	sprintf(format, "%%%ld.%ldFe", dec_prec + 8, dec_prec);

	mpf_set_ui(ad, 2UL); mpf_sqrt(ad, ad);
	printf("ad = "); gmp_printf(format, ad); printf("\n");

	// 変数消去
	mpf_clear(a); // 単独

	// まとめて消去
	mpf_clears(b, c, NULL);       
	mpf_clears(ad, bd, cd, NULL); 

	return 0;
}
