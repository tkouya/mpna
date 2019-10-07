//******************************************************************************
// logistic_qd.c : Sequence calculation based on logistic map with QD library
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
#include "qd_real.h"

using namespace std;

int main()
{
	int i;
	qd_real x[102]; // Quadruple-double precision

	// set a inital value
	x[0] = "0.7501";

	fpu_fix_start(NULL);

	for(i = 0; i <= 100; i++)
	{
		x[i + 1] = 4 * x[i] * (1 - x[i]);
		if((i % 10) == 0)
			cout << setw(5) << i << setprecision(64) << ", " << x[i] << endl;
	}

	return 0;
}
