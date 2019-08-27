//******************************************************************************
// mpn_sample_full.c : Sample code of multiple precision natural number
//               arithmetic (addition, subtraction, multiplication and division)
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
#include <stdlib.h>
#include <string.h>

#include "gmp.h"

#define MAX_LIMB_SIZE 1024
#define MAX_LIMB_SIZE2 2048
#define MAX_STR_LEN 1024
#define MAX_STR_LEN2 2048

int main(void)
{
	int i;
	mp_limb_t a[MAX_LIMB_SIZE], b[MAX_LIMB_SIZE], c[MAX_LIMB_SIZE2];
	mp_limb_t carry, quotient[MAX_LIMB_SIZE], reminder[MAX_LIMB_SIZE];
	mp_size_t size_a, size_b, size_c, size_quotient, size_reminder;
	unsigned char str_a[MAX_STR_LEN] = "", str_b[MAX_STR_LEN] = "", str_c[MAX_STR_LEN2] = "";
	unsigned char str_quotient[MAX_STR_LEN] = "", str_reminder[MAX_STR_LEN] = "";
	size_t strlen_a, strlen_b, strlen_c, strlen_quotient, strlen_reminder;

	// set zeros
	mpn_zero(a, MAX_LIMB_SIZE);
	mpn_zero(b, MAX_LIMB_SIZE);
	mpn_zero(c, MAX_LIMB_SIZE);

	// Input decimal numbers
	printf("a = "); while(scanf("%s", str_a) < 1); printf("str_a = %s\n", str_a);
	printf("b = "); while(scanf("%s", str_b) < 1); printf("str_b = %s\n", str_b);
	//strcpy(str_a, "987654321012345678909876543210123456789098765432101234567890");
	//strcpy(str_b, "123456789098765432101234567890987654321012345678909876543210");

	// convert strings to decimal numbers
	strlen_a = strlen(str_a);
	for(i = 0; i < strlen_a; i++)
		str_a[i] -= '0';

	strlen_b = strlen(str_b);
	for(i = 0; i < strlen_b; i++)
		str_b[i] -= '0';

	// convert to binary numbers
	size_a = mpn_set_str(a, str_a, strlen_a, 10);
	size_b = mpn_set_str(b, str_b, strlen_b, 10);

	// output 
	gmp_printf("a = %Nd\n", a, size_a);
	gmp_printf("b = %Nd\n", b, size_b);

/* Addition */

	// c := a + b
	carry = mpn_add(c, a, size_a, b, size_b);

	// set c's size
	size_c = (size_a > size_b) ? size_a : size_b;
	if(carry >= 1)
	{
		size_c += 1;
		c[size_c - 1] = carry;
	}

	strlen_c = mpn_get_str(str_c, 10, c, size_c);
	mpn_get_str(str_c, 10, c, size_c);
	strlen_c = mpn_sizeinbase(c, size_c, 10);
	printf("c[0] = %lu, carry = %lu, size_c = %d, strlen_c = %d\n", c[0], carry, (int)size_c, (int)strlen_c);

	for(i = 0; i < strlen_c; i++)
		str_c[i] += '0';
	str_c[strlen_c] = '\0';

	gmp_printf("c(%d) = a(%d) + b(%d) = \n[str_c] %s\n [limb] %Nd\n", (int)size_c, (int)size_a, (int)size_b, str_c, c, size_c);

/* Subtraction */

	if(size_a < size_b)
		printf("cannot get a(%d) - b(%d) if a < b!\n", (int)size_a, (int)size_b);
	else
	{
		mpn_zero(c, MAX_LIMB_SIZE2);

		// c := a - b
		carry = mpn_sub(c, a, size_a, b, size_b);

		// set c's size
		size_c = (size_a > size_b) ? size_a : size_b;
		if((carry >= 1) && (size_c > 1))
			size_c -= 1;

	//	printf("c[0] = %lu, carry = %lu, size_c = %d\n", carry, c[0], (int)size_c);

		strlen_c = mpn_get_str(str_c, 10, c, size_c);
		printf("c[0] = %lu, carry = %lu, size_c = %d, strlen_c = %d\n", c[0], carry, (int)size_c, (int)strlen_c);

		for(i = 0; i < strlen_c; i++)
			str_c[i] += '0';
		str_c[strlen_c] = '\0';

		gmp_printf("c(%d) = a(%d) - b(%d) = \n[str_c] %s\n [limb] %Nd\n", (int)size_c, (int)size_a, (int)size_b, str_c, c, size_c);
	}

/* Multiplication */

	// c := a * b
	mpn_zero(c, MAX_LIMB_SIZE2);

	if(size_a < size_b)
		printf("cannot get a(%d) * b(%d) if a < b!\n", (int)size_a, (int)size_b);
	else
	{
		mpn_mul(c, a, size_a, b, size_b);

		// set c's size
		size_c = size_a + size_b;

		strlen_c = mpn_get_str(str_c, 10, c, size_c);
		for(i = 0; i < strlen_c; i++)
			str_c[i] += '0';
		str_c[strlen_c] = '\0';

		gmp_printf("c(%d) = a(%d) * b(%d) = \n[str_c] %s\n [limb] %Nd\n", (int)size_c, (int)size_a, (int)size_b, str_c, c, size_c);
	}

/* Division with Reminder */

	// a / b = quotient ... reminder
	mpn_zero(quotient, MAX_LIMB_SIZE);
	mpn_zero(reminder, MAX_LIMB_SIZE);

	if(size_a < size_b)
		printf("cannot get a(%d) / b(%d) if a < b!\n", (int)size_a, (int)size_b);
	else
	{
		size_quotient = size_a - size_b + 1;
		size_reminder = size_b;
		printf("size_quotient, size_reminder = %d, %d\n", (int)size_quotient, (int)size_reminder);
		mpn_tdiv_qr(quotient, reminder, 0, a, size_a, b, size_b);
		gmp_printf("quotient = %Nd, size_quotient = %d\n", quotient, size_quotient, (int)size_quotient);
		gmp_printf("reminder = %Nd, size_reminder = %d\n", reminder, size_reminder, (int)size_reminder);

		// quotient
		strlen_quotient = mpn_get_str(str_quotient, 10, quotient, size_quotient);
		for(i = 0; i < strlen_quotient; i++)
			str_quotient[i] += '0';
		str_quotient[strlen_quotient] = '\0';

		// reminder
		strlen_reminder = mpn_get_str(str_reminder, 10, reminder, size_reminder);
		for(i = 0; i < strlen_reminder; i++)
			str_reminder[i] += '0';
		str_reminder[strlen_reminder] = '\0';

		printf("a(%d) / b(%d) = %s ... %s\n", (int)size_a, (int)size_b, str_quotient, str_reminder);
		gmp_printf("a(%d) + b(%d) = \n[str_c] %s ... %s\n [limb] %Nd ... %Nd\n", (int)size_a, (int)size_b, str_quotient, str_reminder, quotient, size_quotient, reminder, size_reminder);

		// check : q * b + r == a ?
		mpn_zero(c, MAX_LIMB_SIZE2);

		if(size_quotient > size_b)
			mpn_mul(c, quotient, size_quotient, b, size_b);
		else
			mpn_mul(c, b, size_b, quotient, size_quotient);

		size_c = size_quotient + size_b;
		gmp_printf("q = %Nd\nb = %Nd\nq * b = %Nd\n", quotient, size_quotient, b, size_b, c, size_c);

		carry = mpn_add(c, c, size_c, reminder, size_reminder);
		size_c = (size_c > size_reminder) ? size_c : size_reminder;
		if(carry >= 1)
		{
			size_c += 1;
			c[size_c - 1] = carry;
		}
		gmp_printf("q * b + r == a ? -> %Nd\n", c, size_c);
		if(mpn_cmp(c, a, size_a) == 0)
			printf("OK!\n");
		else
			printf("Not Equal!\n");

	}

	return 0;

}
