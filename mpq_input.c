//******************************************************************************
// mpq_test.c : Input multiple precision rational numbers
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
$ ./mpq_input
a = 3/5
b = 78/47
3/5 + 78/47 = 531/235
3/5 - 78/47 = -249/235
3/5 * 78/47 = 234/235
3/5 / 78/47 = 47/130
*/

int main(void)
{
	mpq_t a, b, c;
	char str_a[1024], str_b[1024];

	mpq_init(a); //mpq_set_str(a, "1 / 3", 10);
	mpq_init(b); //mpq_set_str(b, "2 / 5", 10);
	mpq_init(c);

	// Input
	printf("a = "); gmp_scanf("%Qd", a);
	printf("b = "); gmp_scanf("%Qd", b);

	mpq_add(c, a, b);
	gmp_printf("%Qd + %Qd = %Qd\n", a, b, c);

	mpq_sub(c, a, b);
	gmp_printf("%Qd - %Qd = %Qd\n", a, b, c);

	mpq_mul(c, a, b);
	gmp_printf("%Qd * %Qd = %Qd\n", a, b, c);

	mpq_div(c, a, b);
	gmp_printf("%Qd / %Qd = %Qd\n", a, b, c);

	mpq_clear(a);
	mpq_clear(b);
	mpq_clear(c);

	return 0;
}
