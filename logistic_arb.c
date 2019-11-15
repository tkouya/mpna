//******************************************************************************
// logistic_arb.c : logistic map with Arb
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
#include "arb.h"

int main()
{
	int i;
	slong prec;
	arb_t x[102];

	// initialize arb variables
	for(i = 0; i < 102; i++)
		arb_init(x[i]);

	printf("prec(bits) = "); while(scanf("%ld", &prec) < 1);

	// set a initial value
	arb_set_str(x[0], "0.7501", prec);

	for(i = 0; i <= 100; i++)
	{
		if((i % 10) == 0)
		{
			printf("%5d, ", i);
			arb_printd(x[i], 17);
			printf("\n");
		}

		//x[i + 1] = 4 * x[i] * (1 - x[i]);
		arb_sub_ui(x[i + 1], x[i], 1UL, prec);
		arb_neg(x[i + 1], x[i + 1]);
		arb_mul(x[i + 1], x[i + 1], x[i], prec);
		arb_mul_ui(x[i + 1], x[i + 1], 4UL, prec);

	}

	// delete arb variables
	for(i = 0; i < 102; i++)
		arb_clear(x[i]);

	return 0;
}
