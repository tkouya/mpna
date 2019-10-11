//******************************************************************************
// mpz_input.cpp : Sample code of multiple precision integer
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
#include <iostream>
#include <iomanip>

#include "gmpxx.h"

/*
[tkouya@cs-muse mpna]$ g++ mpz_test.cc -lgmpxx -lgmp
[tkouya@cs-muse mpna]$ ./a.out
123456789012345678901234567890 + 9876543210 = 123456789012345678911111111100
123456789012345678901234567890 - 9876543210 = 123456789012345678891358024680
123456789012345678901234567890 * 9876543210 = 1219326311248285321124828532111263526900
*/

using namespace std;

int main(void)
{
	mpz_class a, b;

	a = "123456789012345678901234567890";
	b = 9876543210;

	cout << a << " + " << b << " = " << a + b << endl;
	cout << a << " - " << b << " = " << a - b << endl;
	cout << a << " * " << b << " = " << a * b << endl;
	cout << a << " / " << b << " = " << a / b << endl;
	cout << a << " % " << b << " = " << a % b << endl;

	return 0;
}
