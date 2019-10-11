//******************************************************************************
// quadratic_eq.c : Algebraic quadratic equation solver (double precision)
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
#include <math.h>
#include <complex.h>

// Sign function
#define SIGN(a) (((a) > 0) ? (+1) : (-1))

int main()
{
	// dcoef[2] * x^2 + dcoef[1] * x + dcoef[0] = 0
	double dcoef[3], dest, dsol[2];
	double complex cdsol[2];

	printf("a = "); while(!scanf("%lf", &dcoef[2]));
	printf("b = "); while(!scanf("%lf", &dcoef[1]));
	printf("c = "); while(!scanf("%lf", &dcoef[0]));

	// d = b^2 - 4ac
	dest = dcoef[1] * dcoef[1] - 4.0 * dcoef[2] * dcoef[0];

	// real solutions
	if(dest >= 0.0)
	{
		// easily loss of siginificant sigits in case of |b| >> |ac|
		//dsol[0] = (-dcoef[1] + sqrt(dest)) / (2 * dcoef[2]);
		//dsol[1] = (-dcoef[1] - sqrt(dest)) / (2 * dcoef[2]);

		// standard way to prevent to loss of siginificant digits
		dsol[0] = (-dcoef[1] - (double)SIGN(dcoef[1]) * sqrt(dest)) / (2 * dcoef[2]);
		dsol[1] = dcoef[0] / (dcoef[2] * dsol[0]);

		printf("x0 = %25.17e\n", dsol[0]);
		printf("x1 = %25.17e\n", dsol[1]);
	}
	// complex solutions
	else
	{
		cdsol[0] = (-(dcoef[1] + 0.0 * I) + csqrt((double complex)(dest + 0.0 * I))) / (2 * (dcoef[2] + 0.0 * I));
		cdsol[1] = (-(dcoef[1] + 0.0 * I) - csqrt((double complex)(dest + 0.0 * I))) / (2 * (dcoef[2] + 0.0 * I));

		printf("x0 = complex(%25.17e, %25.17e)\n", creal(cdsol[0]), cimag(cdsol[0]));
		printf("x1 = complex(%25.17e, %25.17e)\n", creal(cdsol[1]), cimag(cdsol[1]));
	}

	return 0;
}
