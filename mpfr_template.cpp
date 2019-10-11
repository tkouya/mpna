//******************************************************************************
// mpfr_template.cpp : Sample code of mpreal class variables
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
#include "mpreal.h" // MPFR C++

using namespace std;
using namespace mpfr; // MPFR

int main(void)
{
	unsigned long prec;

	cout << "Input default prec in bits: "; cin >> prec;

	// set default precision of mpreal class
	mpreal::set_default_prec(prec);

	mpreal a, b, c; // default precision
	mpreal ad, bd, cd; // default precision

	// print default precision
	cout << "default prec in bits: " << mpreal::get_default_prec() << endl;

	// set twice precision for ad, bd, cd
	ad.set_prec((mp_prec_t)(prec * 2));
	bd.set_prec((mp_prec_t)(prec * 2));
	cd.set_prec((mp_prec_t)(prec * 2));

	// print precision of ad
	cout << "prec(ad) = " << ad.get_prec() << endl;

	//----------------------
	// calculations
	//----------------------

	return 0;
}
