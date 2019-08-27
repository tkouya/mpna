//******************************************************************************
// logistic_f.c : Sequence calculation based on logistic map 
//                                                     with single precision
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

int main()
{
	int i;
	float x[102];

	// set a initial value
	x[0] = 0.7501;

	for(i = 0; i <= 100; i++)
	{
		x[i + 1] = 4 * x[i] * (1 - x[i]);
		if((i % 10) == 0)
			printf("%5d, %15.7e\n", i, x[i]);
	}

	return 0;
}
