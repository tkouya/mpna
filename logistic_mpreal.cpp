//******************************************************************************
// logistic_mpreal.cpp : Sequence calculation based on logistic map with mpreal
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

// mpreal class
#include "mpreal.h"

using namespace std;
using namespace mpfr;

int main(int argc, char *argv[])
{
	int i;
	int num_bits, num_decimal;

	// Check arguments
	if(argc <= 1)
	{
		cerr << "Usage: " << argv[0] << " [num_bits]" << endl;
		return 0;
	}

	// Set default prec in bits
	num_bits = atoi(argv[1]);
	if(num_bits <= 24)
		num_bits = 24;

	num_decimal = (int)ceil(log10(2.0) * (double)num_bits);
	mpreal::set_default_prec(num_bits);

	cout << "num_bits = " << num_bits << ", num_decimal = " << 
	num_decimal << endl;

	mpreal x[102];

	// Set an initial value
	x[0] = "0.7501";

	for(i = 0; i <= 100; i++)
	{
		x[i + 1] = 4 * x[i] * (1 - x[i]);
		if((i % 10) == 0)
			cout << setw(5) << scientific << i << setprecision(num_decimal) <<
			 ", " << x[i] << endl;
	}

	return 0;
}
