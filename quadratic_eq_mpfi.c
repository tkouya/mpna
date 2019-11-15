//******************************************************************************
// quadratic_eq_mpfi.cpp : Quadratic Equation Solver with MPFI
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
#include "mpfr.h" // MPFR
#include "mpfi.h" // MPFI (MPFR + Interval Arithmetic)

int main()
{
	char str[128];
	
	// dcoef[2] * x^2 + dcoef[1] * x + dcoef[0] = 0
	unsigned long prec;
	mpfr_t relerr;
	mpfi_t mcoef[3], mest, tmp, msol[2];
	mpfi_t cmsol_real[2], cmsol_imag[2];

	printf("prec(bits) = "); while(scanf("%ld", &prec) < 1);

	// set prec as default precision
	mpfr_set_default_prec((mp_prec_t)prec);

	mpfr_init(relerr);
	mpfi_init(mcoef[0]);
	mpfi_init(mcoef[1]);
	mpfi_init(mcoef[2]);
	mpfi_init(tmp);
	mpfi_init(msol[0]);
	mpfi_init(msol[1]);
	mpfi_init(mest);

	mpfi_init(cmsol_real[0]);
	mpfi_init(cmsol_imag[0]);
	mpfi_init(cmsol_real[1]);
	mpfi_init(cmsol_imag[1]);

	printf("a = "); while(scanf("%s", str) < 1); mpfi_set_str(mcoef[2], str, 10);
	printf("b = "); while(scanf("%s", str) < 1); mpfi_set_str(mcoef[1], str, 10);
	printf("c = "); while(scanf("%s", str) < 1); mpfi_set_str(mcoef[0], str, 10);

	// d = b^2 - 4ac
	//dest = dcoef[1] * dcoef[1] - 4.0 * dcoef[2] * dcoef[0];
	mpfi_mul(mest, mcoef[1], mcoef[1]);
	mpfi_mul(tmp, mcoef[2], mcoef[0]);
	mpfi_mul_ui(tmp, tmp, 4UL);
	mpfi_sub(mest, mest, tmp);

	// 
	// if(dest >= 0)
	if(mpfi_is_pos(mest) > 0)
	{
		//dsol[0] = (-dcoef[1] - (double)SIGN(dcoef[1]) * sqrt(dest)) / (2 * dcoef[2]);
		mpfi_sqrt(tmp, mest);

		//sign(tmp, mcoef[1]);
		if(mpfi_is_neg(mcoef[1]) > 0) 
			mpfi_neg(tmp, tmp);

		mpfi_neg(msol[0], mcoef[1]);
		mpfi_sub(msol[0], msol[0], tmp);
		mpfi_mul_ui(tmp, mcoef[2], 2UL);
		mpfi_div(msol[0], msol[0], tmp);

		//dsol[1] = dcoef[0] / (dcoef[2] * dsol[0]);
		mpfi_mul(tmp, mcoef[2], msol[0]);
		mpfi_div(msol[1],mcoef[0], tmp);
	
	
		printf("x0 = "); mpfi_out_str(stdout, 10, 17, msol[0]);
		mpfi_diam_rel(relerr, msol[0]); mpfr_printf(" %10.3RNe\n", relerr);
		printf("x1 = "); mpfi_out_str(stdout, 10, 17, msol[1]);
		mpfi_diam_rel(relerr, msol[1]); mpfr_printf(" %10.3RNe\n", relerr);
	}
	// Complex solutions
	else
	{
		//cdsol[0] = (-(dcoef[1] + 0.0 * I) + csqrt((double complex)(dest + 0.0 * I))) / (2 * (dcoef[2] + 0.0 * I));
		//cdsol[1] = (-(dcoef[1] + 0.0 * I) - csqrt((double complex)(dest + 0.0 * I))) / (2 * (dcoef[2] + 0.0 * I));
		// Real part
		mpfi_neg(cmsol_real[0], mcoef[1]);
		mpfi_neg(cmsol_real[1], mcoef[1]);
		mpfi_mul_ui(tmp, mcoef[2], 2UL);
		mpfi_div(cmsol_real[0], cmsol_real[0], tmp);
		mpfi_div(cmsol_real[1], cmsol_real[1], tmp);

		// Imaginary part
		mpfi_neg(cmsol_imag[0], mest);
		mpfi_sqrt(cmsol_imag[0], cmsol_imag[0]);
		mpfi_div(cmsol_imag[0], cmsol_imag[0], tmp);
		mpfi_neg(cmsol_imag[1], cmsol_imag[0]);

		printf("x0 = complex("); mpfi_out_str(stdout, 10, 17, cmsol_real[0]);
		printf(", "); mpfi_out_str(stdout, 10, 17, cmsol_imag[0]); printf(")\n");
		printf("x1 = complex("); mpfi_out_str(stdout, 10, 17, cmsol_real[1]);
		printf(", "); mpfi_out_str(stdout, 10, 17, cmsol_imag[1]); printf(")\n");
	}

	// Clear all variables
	mpfr_clear(relerr);
	mpfi_clear(mcoef[0]);
	mpfi_clear(mcoef[1]);
	mpfi_clear(mcoef[2]);
	mpfi_clear(tmp);
	mpfi_clear(msol[0]);
	mpfi_clear(msol[1]);
	mpfi_clear(mest);
	mpfi_clear(cmsol_real[0]);
	mpfi_clear(cmsol_real[1]);
	mpfi_clear(cmsol_imag[0]);
	mpfi_clear(cmsol_imag[1]);

	return 0;
}
