//******************************************************************************
// mpz_input_gcd_lcm.c : Sample code to get GCD and LCM of multiple precision 
// integers
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
//******************************************************************************
#include <stdio.h>

#include "gmp.h"

/*
[tkouya@cs-muse mpna]$ gcc mpz_test.c -lgmp -lgmp
[tkouya@cs-muse mpna]$ ./a.out
123456789012345678901234567890 + 9876543210 = 123456789012345678911111111100
123456789012345678901234567890 - 9876543210 = 123456789012345678891358024680
123456789012345678901234567890 * 9876543210 = 1219326311248285321124828532111263526900
123456789012345678901234567890 / 9876543210 = 12499999887343749990
123456789012345678901234567890 % 9876543210 = 1562499990
*/

int main(void)
{
	mpz_t a, b, c, d, e, gcd_ab, lcm_ab;
	unsigned char str_a[1024], str_b[1024];

//	mpz_init_set_str(a, "123456789012345678901234567890", 10);
//	mpz_init_set_ui(b, 9876543210);
	mpz_init(a);
	mpz_init(b);
	mpz_init(c);
	mpz_init(d);
	mpz_init(e);
	mpz_init(gcd_ab);
	mpz_init(lcm_ab);

	// input
	printf("a = ");
#ifdef USE_STR_INPUT
	// string to mpz
	while(scanf("%s", str_a) < 1);
	printf("str_a = %s\n", str_a);
	mpz_set_str(a, str_a, 10);
#else // USE_STR_INPUT
	gmp_scanf("%Zd", a);
#endif // USE_STR_INPUT

	gmp_printf("a = %Zd\n", a);

	printf("b = ");
#ifdef USE_STR_INPUT
	// string to mpz
	while(scanf("%s", str_b) < 1);
	printf("str_b = %s\n", str_a);
	mpz_set_str(b, str_b, 10);
#else // USE_STR_INPUT
	gmp_scanf("%Zd", b);
#endif // USE_STR_INPUT

	mpz_gcd(gcd_ab, a, b);
	gmp_printf("GCD(%Zd, %Zd) = %Zd\n", a, b, gcd_ab);

	mpz_lcm(lcm_ab, a, b);
	gmp_printf("LCM(%Zd, %Zd) = %Zd\n", a, b, lcm_ab);

// validation(1)
	mpz_mod(c, a, gcd_ab);
	gmp_printf("%Zd mod GCD(%Zd, %Zd) = %Zd\n", a, a, b, c);
	mpz_mod(c, a, gcd_ab);
	gmp_printf("%Zd mod GCD(%Zd, %Zd) = %Zd\n", a, a, b, c);
	mpz_mod(c, lcm_ab, a);
	gmp_printf("LCM(%Zd, %Zd) mod %Zd = %Zd\n", a, b, a, c);
	mpz_mod(c, lcm_ab, b);
	gmp_printf("LCM(%Zd, %Zd) mod %Zd = %Zd\n", a, b, b, c);

// validation(2)
// (a / GCD(a, b)) * (b / GCD(a, b)) * GCD(a, b) = LCM(a, b)
	mpz_div(c, a, gcd_ab);
	mpz_div(d, b, gcd_ab);
	mpz_mul(e, c, d);
	mpz_mul(e, e, gcd_ab);
	gmp_printf("(a / GCD(a, b)) * (b / GCD(a, b)) * GCD(a, b) = %Zd\n", e);
	gmp_printf("LCM(a, b)                                     = %Zd\n", lcm_ab);

	mpz_clear(a);
	mpz_clear(b);
	mpz_clear(c);
	mpz_clear(d);
	mpz_clear(e);
	mpz_clear(gcd_ab);
	mpz_clear(lcm_ab);

	return 0;
}
