// mpf_template.cpp: mpf_class型の使い方

#include <iostream>
#include <iomanip>
#include "gmpxx.h"

using namespace std;

int main(void)
{
	unsigned long prec;

	cout << "Input default prec in bits: "; cin >> prec;

	// デフォルトの仮数部ビット数をセット
	mpf_set_default_prec((mp_bitcnt_t)prec); // in bits

	mpf_class a, b, c; // デフォルトの精度
	mpf_class ad, bd, cd; // a, b, cより長い桁

	// デフォルトの精度を表示
	cout << "default prec in bits: " << prec << endl;

	// デフォルト仮数部の2倍に設定
	ad.set_prec((mp_bitcnt_t)(prec * 2));
	bd.set_prec((mp_bitcnt_t)(prec * 2));
	cd.set_prec((mp_bitcnt_t)(prec * 2));

	// 変数ごとの仮数部ビット数を表示
	cout << "prec(ad) = " << ad.get_prec() << endl;

	//a = sqrt(2);
	//b = sqrt(3);
	a = sqrt((mpf_class)2UL);
	b = sqrt((mpf_class)3UL);

	c = a + b;
	gmp_printf("%50.43Fe + %50.43Fe = %50.43Fe\n", a.get_mpf_t(), b.get_mpf_t(), c.get_mpf_t());

	c = a - b;
	gmp_printf("%50.43Fe - %50.43Fe = %50.43Fe\n", a.get_mpf_t(), b.get_mpf_t(), c.get_mpf_t());

	c = a * b;
	gmp_printf("%50.43Fe * %50.43Fe = %50.43Fe\n", a.get_mpf_t(), b.get_mpf_t(), c.get_mpf_t());

	c = a / b;
	gmp_printf("%50.43Fe / %50.43Fe = %50.43Fe\n", a.get_mpf_t(), b.get_mpf_t(), c.get_mpf_t());

	c = sqrt(a);
	gmp_printf("sqrt(%50.43Fe) = %50.43Fe\n", a.get_mpf_t(), c.get_mpf_t());

	// 変数消去
	~a; // 単独

	// まとめて消去
	~b, ~c;
	~ad, ~bd, ~cd; 

	return 0;
}
