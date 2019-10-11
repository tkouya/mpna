//******************************************************************************
// mpz_test.c : First sample code of multiple precision integer
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

#include "gmp.h"

int main(void)
{
	mpz_t a, b, c;
	
	mpz_init_set_str(a, "123456789012345678901234567890", 10);
	mpz_init_set_ui(b, 9876543210);
	mpz_init(c);
	
	mpz_add(c, a, b);
	gmp_printf("%Zd + %Zd = %Zd\n", a, b, c);
	
	mpz_sub(c, a, b);
	gmp_printf("%Zd - %Zd = %Zd\n", a, b, c);
	
	mpz_mul(c, a, b);
	gmp_printf("%Zd * %Zd = %Zd\n", a, b, c);
	
	mpz_div(c, a, b);
	gmp_printf("%Zd / %Zd = %Zd\n", a, b, c);
	
	mpz_mod(c, a, b);
	gmp_printf("%Zd %% %Zd = %Zd\n", a, b, c);
	
	mpz_clear(a);
	mpz_clear(b);
	mpz_clear(c);
	
	return 0;
}
