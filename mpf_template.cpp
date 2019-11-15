//******************************************************************************
// mpf_template.cpp : Sample code of mpf_class variables
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
#include "gmpxx.h"

using namespace std;

int main(void)
{
	unsigned long prec;

	cout << "Input default prec in bits: "; cin >> prec;

	// set default precision
	mpf_set_default_prec((mp_bitcnt_t)prec); // in bits

	mpf_class a, b, c; // default precision
	mpf_class ad, bd, cd; // default precision

	// print default precision
	cout << "default prec in bits: " << prec << endl;

	// change twice precision for ad, bd, and cd
	ad.set_prec((mp_bitcnt_t)(prec * 2));
	bd.set_prec((mp_bitcnt_t)(prec * 2));
	cd.set_prec((mp_bitcnt_t)(prec * 2));

	// print each precision of variable
	cout << "prec(ad) = " << ad.get_prec() << endl;

	//----------------------
	// calculations
	//----------------------

	// delete one variable
	delete a; // single

	// delete variables
	delete b, c;
	delete ad, bd, cd;

	return 0;
}
