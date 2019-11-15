//******************************************************************************
// complex_d.cpp : Test program of double complex
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
#include <complex>

using namespace std;

int main(int argc, char* argv[])
{
	// double precision complex
  	complex<double> a, b;

	// input a and b from standard input
	// ex : "(real_part, imag_part)"
	cout << "Input a ->";
	cin >> a;
	cout << "Input b ->";
  	cin >> b;

	// print a and b
	cout << "a = " << a.real() << " + " << a.imag() << " * I" << endl;
	cout << "b = " << b.real() << " + " << b.imag() << " * I" << endl;

	// output a + b, a - b, a * b, and a / b
	cout << a << " + " << b << " = " << a + b << endl;
  	cout << a << " - " << b << " = " << a - b << endl;
  	cout << a << " * " << b << " = " << a * b << endl;
  	cout << a << " / " << b << " = " << a / b << endl;

	return 0;
}
