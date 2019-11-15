//******************************************************************************
// mpz_input_nextprime.cpp : Get next prime number 
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
#include <cstdlib>
#include "gmpxx.h"

using namespace std;

/*
$ g++ mpz_input_nextprime.cpp -lgmpxx -lgmp
$ ./mpz_input_nextprime_cxx
a = 246857968
a = 246857968
246857968 = 2 * 2 * 2 * 2 * 7 * 73 * 109 * 277
*/

int main(void)
{
	int num_prime;
	mpz_class a, prime, c;

	// input
	cout << "a = ";
	cin >> a;
	cout << "a = " << a << endl;

	if(a < 1UL)
    {
        cerr << "a = " << a << " is less than zero!" << endl;
        exit(EXIT_FAILURE);
    }

    // prime := 2
    prime = 2UL;
    num_prime = 0;
    while(a > prime)
    {
        cout << ++num_prime << ": " << prime << endl;
        mpz_nextprime(prime.get_mpz_t(), prime.get_mpz_t());
    }

    cout << "<= " << a << endl;

	return 0;
}
