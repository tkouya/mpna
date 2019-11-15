//******************************************************************************
// float128.c : Sample program to use __float128 with GCC
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
#include <quadmath.h> // GCC Quad-Math Library

int main(void)
{
	char input_str[128], output_str[128];

  	__float128 a, b, c;

	// input a and b from standard input
	printf("Input a ->"); while(scanf("%s", input_str) < 1); a = strtoflt128(input_str, NULL);
	printf("Input b ->"); while(scanf("%s", input_str) < 1); b = strtoflt128(input_str, NULL);

	// print a and b
	quadmath_snprintf(output_str, sizeof(output_str), "%+Qg", a);
	printf("a = %s\n", output_str);
	quadmath_snprintf(output_str, sizeof(output_str), "%+Qa", a);
	printf("a = %s\n", output_str);

	quadmath_snprintf(output_str, sizeof(output_str), "%+Qg", b);
	printf("b = %s\n", output_str);

	// print a + b, a - b, a * b, a / b, and sqrt(a)
	c = a + b;
	quadmath_snprintf(output_str, sizeof(output_str), "%+.34Qg", c);
	printf("a + b = %s\n", output_str);

	c = a - b;
	quadmath_snprintf(output_str, sizeof(output_str), "%+.34Qg", c);
	printf("a - b = %s\n", output_str);

	c = a * b;
	quadmath_snprintf(output_str, sizeof(output_str), "%+.34Qg", c);
	printf("a * b = %s\n", output_str);

	c = a / b;
	quadmath_snprintf(output_str, sizeof(output_str), "%+.34Qg", c);
	printf("a / b = %s\n", output_str);

	c = sqrtq(a);
	quadmath_snprintf(output_str, sizeof(output_str), "%+.34Qg", c);
	printf("sqrt(a) = %s\n", output_str);

	c = c * c;
	quadmath_snprintf(output_str, sizeof(output_str), "%+.34Qg", c);
	printf("(sqrt(a))^2 = %s\n", output_str);

	// end
	return 0;
}
