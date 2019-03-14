// mpf_template.c: mpf_t型の使い方

#include <stdio.h>
#include "gmp.h"

int main(void)
{
	unsigned long prec;
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

	// ・・・・・・・・・・
	// ・・・演算処理・・・
	// ・・・・・・・・・・

	// 変数消去
	mpf_clear(a); // 単独

	// まとめて消去
	mpf_clears(b, c, NULL);       
	mpf_clears(ad, bd, cd, NULL); 

	return 0;
}
