//******************************************************************************
// quadratic_eq_mpfi.cpp : Algebraic quadratic equation solver with ARB
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
#include "arb.h" // Arb : Real computation
#include "acb.h" // Acb : Complex Computation

int main()
{
	char str[128];
	
	// dcoef[2] * x^2 + dcoef[1] * x + dcoef[0] = 0
	slong prec;
	arb_t acoef[3], aest, tmp, asol[2];
	acb_t casol[2], ctmp;

	// Initialize all variables
	arb_init(acoef[0]);
	arb_init(acoef[1]);
	arb_init(acoef[2]);
	arb_init(tmp);
	arb_init(asol[0]);
	arb_init(asol[1]);
	arb_init(aest);
	acb_init(casol[0]);
	acb_init(casol[1]);
	acb_init(ctmp);

	printf("prec(bits) = "); while(scanf("%ld", &prec) < 1);
	printf("a = "); while(scanf("%s", str) < 1); arb_set_str(acoef[2], str, prec);
	printf("b = "); while(scanf("%s", str) < 1); arb_set_str(acoef[1], str, prec);
	printf("c = "); while(scanf("%s", str) < 1); arb_set_str(acoef[0], str, prec);

	// d = b^2 - 4ac
	//dest = dcoef[1] * dcoef[1] - 4.0 * dcoef[2] * dcoef[0];
	arb_mul(aest, acoef[1], acoef[1], prec);
	arb_mul(tmp, acoef[2], acoef[0], prec);
	arb_mul_ui(tmp, tmp, 4UL, prec);
	arb_sub(aest, aest, tmp, prec);

	// Real solutions
	// if(dest >= 0)
	if(arb_contains_nonnegative(aest) != 0)
	{
		//dsol[0] = (-dcoef[1] - (double)SIGN(dcoef[1]) * sqrt(dest)) / (2 * dcoef[2]);
		arb_sqrt(asol[0], aest, prec);
		arb_sgn(tmp, acoef[1]);
		arb_mul(tmp, tmp, asol[0], prec);
		arb_neg(asol[0], acoef[1]);
		arb_sub(asol[0], asol[0], tmp, prec);
		arb_mul_ui(tmp, acoef[2], 2UL, prec);
		arb_div(asol[0], asol[0], tmp, prec);

		//dsol[1] = dcoef[0] / (dcoef[2] * dsol[0]);
		arb_mul(tmp, acoef[2], asol[0], prec);
		arb_div(asol[1],acoef[0], tmp, prec);
		
		printf("x0 = "); arb_printd(asol[0], 17); printf("\n");
		printf("x1 = "); arb_printd(asol[1], 17); printf("\n");
	}
	// Complex solutions
	else
	{
		//cdsol[0] = (-(dcoef[1] + 0.0 * I) + csqrt((double complex)(dest + 0.0 * I))) / (2 * (dcoef[2] + 0.0 * I));
		//cdsol[1] = (-(dcoef[1] + 0.0 * I) - csqrt((double complex)(dest + 0.0 * I))) / (2 * (dcoef[2] + 0.0 * I));
		arb_neg(tmp, acoef[1]);
		acb_set_arb(casol[0], tmp);
		acb_set_arb(casol[1], tmp);
		
		acb_set_arb(ctmp, aest);
		acb_sqrt(ctmp, ctmp, prec);
		acb_add(casol[0], casol[0], ctmp, prec);
		acb_sub(casol[1], casol[1], ctmp, prec);
		
		arb_mul_ui(tmp, acoef[2], 2UL, prec);
		acb_div_arb(casol[0], casol[0], tmp, prec);
		acb_div_arb(casol[1], casol[1], tmp, prec);

		printf("x0 = "); acb_printd(casol[0], 17); printf("\n");
		printf("x1 = "); acb_printd(casol[1], 17); printf("\n");
	}

	// Clear all variables
	arb_clear(acoef[0]);
	arb_clear(acoef[1]);
	arb_clear(acoef[2]);
	arb_clear(tmp);
	arb_clear(asol[0]);
	arb_clear(asol[1]);
	arb_clear(aest);
	acb_clear(casol[0]);
	acb_clear(casol[1]);
	acb_clear(ctmp);

	return 0;
}
