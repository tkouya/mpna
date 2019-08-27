//******************************************************************************
// logistic_err.c : Sequence calculation based on logistic map 
//                                                        with error estimations
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
#include <fenv.h> // fegetround, fesetround

// Print current rounding mode
int rmode_view(void)
{
	int current_rmode;

	current_rmode = fegetround();

	printf("--- Rounding Mode ");
	switch(current_rmode)
	{
		case FE_TONEAREST:	// RN mode
			printf("RN mode"); break;
		case FE_UPWARD: 	// RP mode
			printf("RP mode"); break;
		case FE_DOWNWARD:	// RM mode
			printf("RM mode"); break;
		case FE_TOWARDZERO: // RZ mode
			printf("RZ mode"); break;
		default:
			printf("Unknown mode");
	}
	printf(" ---\n");

	return current_rmode;
}

/*
$ ./logistic_err
--- Rounding Mode RN mode ---
--- Rounding Mode RP mode ---
--- Rounding Mode RM mode ---
--- Rounding Mode RN mode ---
  i             RM                     RN                       RP
  0,   7.501000000000000e-01   7.501000000000000e-01   7.501000000000000e-01
 10,   8.444959536022354e-01   8.444959536022012e-01   8.444959536022033e-01
 20,   1.429397244947286e-01   1.429397245283995e-01   1.429397245262345e-01
 30,   8.542959855593284e-01   8.542960203146587e-01   8.542960180805079e-01
 40,   7.749537600698271e-01   7.749958851552057e-01   7.749931773385752e-01
 50,   1.096498172466455e-01   7.951287645010524e-02   8.131852429545657e-02
 60,   9.119982359028228e-04   2.731872404408921e-01   5.227061527431969e-01
 70,   2.191116324728341e-01   5.525305620833623e-01   9.110545209083801e-01
 80,   4.050266076761563e-01   2.162556639958134e-01   6.175123937364528e-01
 90,   1.954179677741432e-01   7.874679371884123e-01   7.551369172230279e-01
100,   5.019031920912596e-01   2.697067458876520e-01   5.733967630679918e-01*/

int main()
{
	int i;
	int default_rmode;
	double x[102], x_rp[102], x_rm[102];

	// print the current rounding mode
	default_rmode = rmode_view();

	// set a initial value
	x[0] = 0.7501;
	for(i = 0; i <= 100; i++)
		x[i + 1] = 4 * x[i] * (1 - x[i]);

	// RP mode
	fesetround(FE_UPWARD); rmode_view();

	x_rp[0] = 0.7501;
	for(i = 0; i <= 100; i++)
		x_rp[i + 1] = 4 * x_rp[i] * (1 - x_rp[i]);
	// End of RP mode

	// RM mode
	fesetround(FE_DOWNWARD); rmode_view();

	x_rm[0] = 0.7501;
	for(i = 0; i <= 100; i++)
		x_rm[i + 1] = 4 * x_rm[i] * (1 - x_rm[i]);

	// back to default rounding mode
	fesetround(default_rmode); rmode_view();

	printf("  i             RM                     RN                       RP\n");
	for(i = 0; i <= 100; i++)
	{
		if((i % 10) == 0)
			printf("%3d, %23.15e %23.15e %23.15e\n", i, x_rm[i], x[i], x_rp[i]);
	}

	return 0;
}
