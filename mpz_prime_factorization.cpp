//******************************************************************************
// mpz_prime_factorization.cpp : Prime factorization of mpz_class
// integers
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
#include "gmpxx.h"

using namespace std;

/*
$ g++ mpz_prime_factorization.cpp -lgmpxx -lgmp
$ ./a.out
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

	// prime := 2
	prime = 2UL;
	cout << a << " = ";
	num_prime = 0;
	while(a >= prime)
	{
		while(1)
		{
			c = a % prime;
			if(c != 0)
				break;

			a /= prime;
			num_prime++;
			if(num_prime >= 2)
				cout << " * ";

			cout << prime;
		}

		mpz_nextprime(prime.get_mpz_t(), prime.get_mpz_t());
	}
	printf("\n");

	return 0;
}
