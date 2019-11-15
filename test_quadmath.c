//******************************************************************************
// test_quadmath.c : Test program of quadmath library with GCC
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
#include <quadmath.h>

int main()
{
	__float128 a, b, c;
	char buf[1024];

	a = sqrtq(2.0q);
	b = sqrtq(3.0q);

	quadmath_snprintf(buf, 1024, "%50.33Qe", a);
	printf("a = %s\n", buf);
	quadmath_snprintf(buf, 1024, "%50.33Qe", b);
	printf("b = %s\n", buf);

	c = a + b;

	quadmath_snprintf(buf, 1024, "%Qe", c);
	printf("c = %s\n", buf);

	return 0;
}