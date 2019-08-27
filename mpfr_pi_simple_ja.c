//******************************************************************************
// mpfr_pi_simple.c : PI calculation program with MPFR/GMP
// Copyright (C) 2019 Tomonori Kouya
// 
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation, either version 3 of the License or any later
// version.
// 
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// [usage] $ gcc mpfr_pi_simple.c -lmpfr -lgmp -lm
//         $ ./a.out
//           10進精度桁数を入力: 1000
//           10進精度桁数 = 1000
//            2進精度桁数 = 3321
// pi(3321 bits) = 3.14159265358979323846264338327950288419716939937510582097494
//******************************************************************************
#include <stdio.h>
#include <math.h>
#include "mpfr.h" // MPFR/GMP

int main()
{
	mpfr_prec_t dec_prec, prec;
	mpfr_t mp_pi;

	printf("10進精度桁数を入力: "); while(scanf("%ld", &dec_prec) < 1);

	prec = (mpfr_prec_t)ceil(dec_prec / log10(2.0)); // 計算ビット数算出

	printf("10進精度桁数 = %ld\n", dec_prec);
	printf(" 2進精度桁数 = %ld\n", prec);

	mpfr_init2(mp_pi, prec); // 変数初期化

	mpfr_const_pi(mp_pi, MPFR_RNDN); // piを計算
	mpfr_printf("pi(%ld bits) = %RNe\n", prec, mp_pi); // 出力

	mpfr_clear(mp_pi); // 変数消去

	return 0;
}
