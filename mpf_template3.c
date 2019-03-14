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

	//a = sqrt(2);
	//b = sqrt(3);
	mpf_set_ui(a, 2UL); mpf_sqrt(a, a);
	mpf_set_ui(b, 3UL); mpf_sqrt(b, b);
	
	mpf_add(c, a, b); // c = a + b;
	gmp_printf("%50.43Fe + %50.43Fe = %50.43Fe\n", a, b, c);
	
	mpf_sub(c, a, b); // c = a - b;
	gmp_printf("%50.43Fe - %50.43Fe = %50.43Fe\n", a, b, c);
	
	mpf_mul(c, a, b); // c = a * b;
	gmp_printf("%50.43Fe * %50.43Fe = %50.43Fe\n", a, b, c);
	
	mpf_div(c, a, b); // c = a / b;
	gmp_printf("%50.43Fe / %50.43Fe = %50.43Fe\n", a, b, c);

	mpf_sqrt(c, a); // c = sqrt(a);
	gmp_printf("sqrt(%50.43Fe) = %50.43Fe\n", a, c);

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
