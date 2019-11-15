//******************************************************************************
// complex128.c : Sample program to use __complex128 with GCC
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
	char input_str[128], output_str[2][128];

	// both real and imaginary parts are __float128
  	__complex128 a, b, c;
  	__float128 a_real, a_imag, b_real, b_imag;

	// input a and b from standard input
	printf("Input Re(a) ->"); while(scanf("%s", input_str) < 1); a_real = strtoflt128(input_str, NULL);
	printf("Input Im(a) ->"); while(scanf("%s", input_str) < 1); a_imag = strtoflt128(input_str, NULL);

	// a := (a_real) + (a_imag) * I
	__real__ a = a_real; __imag__ a = a_imag; // GCC dependant instruction

	printf("Input Re(b) ->"); while(scanf("%s", input_str) < 1); b_real = strtoflt128(input_str, NULL);
	printf("Input Im(b) ->"); while(scanf("%s", input_str) < 1); b_imag = strtoflt128(input_str, NULL);

	// b := (b_real) + (b_imag) * I
	__real__ b = b_real; __imag__ b = b_imag;

	// print a and b separately 
	quadmath_snprintf(output_str[0], sizeof(output_str[0]), "%+Qg", crealq(a));
	quadmath_snprintf(output_str[1], sizeof(output_str[1]), "%+Qg", cimagq(a));
	printf("a = (%s) + (%s) * I\n", output_str[0], output_str[1]);

	quadmath_snprintf(output_str[0], sizeof(output_str[0]), "%+Qg", crealq(b));
	quadmath_snprintf(output_str[1], sizeof(output_str[1]), "%+Qg", cimagq(b));
	printf("b = (%s) + (%s) * I\n", output_str[0], output_str[1]);

	// print a + b, a - b, a * b, a / b, sqrt(a), and sqrt(a)^2
	c = a + b;
	quadmath_snprintf(output_str[0], sizeof(output_str[0]), "%+.34Qg", crealq(c));
	quadmath_snprintf(output_str[1], sizeof(output_str[1]), "%+.34Qg", cimagq(c));
	printf("a + b = (%s) + (%s) * I\n", output_str[0], output_str[1]);

	c = a - b;
	quadmath_snprintf(output_str[0], sizeof(output_str[0]), "%+.34Qg", crealq(c));
	quadmath_snprintf(output_str[1], sizeof(output_str[1]), "%+.34Qg", cimagq(c));
	printf("a - b = (%s) + (%s) * I\n", output_str[0], output_str[1]);

	c = a * b;
	quadmath_snprintf(output_str[0], sizeof(output_str[0]), "%+.34Qg", crealq(c));
	quadmath_snprintf(output_str[1], sizeof(output_str[1]), "%+.34Qg", cimagq(c));
	printf("a * b = (%s) + (%s) * I\n", output_str[0], output_str[1]);

	c = a / b;
	quadmath_snprintf(output_str[0], sizeof(output_str[0]), "%+.34Qg", crealq(c));
	quadmath_snprintf(output_str[1], sizeof(output_str[1]), "%+.34Qg", cimagq(c));
	printf("a / b = (%s) + (%s) * I\n", output_str[0], output_str[1]);

	c = csqrtq(a);
	quadmath_snprintf(output_str[0], sizeof(output_str[0]), "%+.34Qg", crealq(c));
	quadmath_snprintf(output_str[1], sizeof(output_str[1]), "%+.34Qg", cimagq(c));
	printf("csqrt(a) = (%s) + (%s) * I\n", output_str[0], output_str[1]);

	c = c * c;
	quadmath_snprintf(output_str[0], sizeof(output_str[0]), "%+.34Qg", crealq(c));
	quadmath_snprintf(output_str[1], sizeof(output_str[1]), "%+.34Qg", cimagq(c));
	printf("(csqrt(a))^2 = (%s) + (%s) * I\n", output_str[0], output_str[1]);

	return 0;
}
