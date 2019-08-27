//******************************************************************************
// logistic_f_err.c : Sequence calculation based on logistic map 
//                                   with single precisin and error estimations
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
$ ./logistic_f_err
--- Rounding Mode RN mode ---
--- Rounding Mode RP mode ---
--- Rounding Mode RM mode ---
--- Rounding Mode RN mode ---
  i        RM             RN              RP
  0,  7.5010002e-01  7.5010002e-01  7.5010002e-01
 10,  8.4451348e-01  8.4451681e-01  8.4449655e-01
 20,  1.2610383e-01  1.2290392e-01  1.4227618e-01
 30,  7.4164075e-01  4.9778005e-01  4.2725600e-02
 40,  3.5835084e-01  5.7825547e-01  1.7286229e-01
 50,  3.1159574e-01  8.7997538e-01  4.6064019e-01
 60,  1.2042288e-02  6.9623333e-01  2.2535935e-01
 70,  2.2944739e-01  5.6387091e-01  4.3607539e-01
 80,  2.8965032e-01  1.4486660e-01  1.0813350e-01
 90,  4.4908887e-01  5.6300414e-01  3.9492986e-01
100,  8.5733098e-01  9.2065656e-01  9.9978036e-01
*/

int main()
{
	int i;
	int default_rmode;
	float x[102], x_rp[102], x_rm[102];

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

	printf("  i        RM             RN              RP\n");
	for(i = 0; i <= 100; i++)
	{
		if((i % 10) == 0)
			printf("%3d, %14.7e %14.7e %14.7e\n", i, x_rm[i], x[i], x_rp[i]);
	}

	return 0;
}
