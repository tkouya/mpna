//******************************************************************************
// complex_d.c : Test program of double complex
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
#include <complex.h> // C99 complex type

int main(int argc, char* argv[])
{
	// all variable are double precision
  	double complex a, b, c;
  	double a_real, a_imag, b_real, b_imag;

	// input a and b b from standard input
	printf("Input Re(a) ->"); while(scanf("%lf", &a_real) < 1);
	printf("Input Im(a) ->"); while(scanf("%lf", &a_imag) < 1);
	a = a_real + a_imag * I; // I = sqrt(-1)

	printf("Input Re(b) ->"); while(scanf("%lf", &b_real) < 1);
	printf("Input Im(b) ->"); while(scanf("%lf", &b_imag) < 1);
	b = b_real + b_imag * I;

	// print a and b 
	printf("a = (%+g) + (%+g) * I\n", creal(a), cimag(a));
	printf("b = (%+g) + (%+g) * I\n", creal(b), cimag(b));

	// print a + b, a - b, a * b
	c = a + b;
	printf("a + b = (%+g) + (%+g) * I\n", creal(c), cimag(c));

	c = a - b;
	printf("a - b = (%+g) + (%+g) * I\n", creal(c), cimag(c));

	c = a * b;
	printf("a * b = (%+g) + (%+g) * I\n", creal(c), cimag(c));

	c = a / b;
	printf("a / b = (%+g) + (%+g) * I\n", creal(c), cimag(c));

	return 0;
}
