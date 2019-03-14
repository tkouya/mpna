/********************************************************************************/
/* dd_test.cc:Double-double and Quadruple precision Linear Computation Library  */
/* Copyright (C) 2015-2016 Tomonori Kouya                                       */
/*                                                                              */
/* This program is free software: you can redistribute it and/or modify it      */
/* under the terms of the GNU Lesser General Public License as published by the */
/* Free Software Foundation, either version 3 of the License or any later       */
/* version.                                                                     */
/*                                                                              */
/* This program is distributed in the hope that it will be useful, but WITHOUT  */
/* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        */
/* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License */
/* for more details.                                                            */
/*                                                                              */
/* You should have received a copy of the GNU Lesser General Public License     */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>.        */
/*                                                                              */
/********************************************************************************/
#include <iostream>
#include <iomanip>

#define QD_INLINE
#include "qd/qd_real.h"
#include "qd/fpu.h"

using namespace std;

#include <cstdio>
#include <cstdlib>
#include <cstring>

extern "C" {
// unsigned char to binary
unsigned char *uc2b(unsigned char *str_ret, unsigned char c)
{
	unsigned tmp_uc;
	unsigned char table[16][5] = {"0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111", "1000", "1001", "1010", "1011", "1100", "1101", "1110", "1111"};

	//printf("c / 16 = %d\n", (unsigned char)(c / 16));
	//printf("c % 16 = %d\n", (unsigned char)(c % 16));

	strcpy((char *)str_ret, (char *)table[(unsigned char)(c / 16)]);
	strcpy((char *)(str_ret + 4), (char *)table[(unsigned char)(c % 16)]);

	str_ret[9] = '\0';

	return str_ret;
}

// print float in binary format
void printb_float(float val)
{
	unsigned char binary[4][16];
	// union to print float data type
	union {
		unsigned char uc_val_f[4]; // 4 bytes = 32 bits
		float val_f;
	} tmp_union_f;

	tmp_union_f.val_f = val;

	// Big Endian
	// sign(1) exp(8) frac(23)+1
	printf("%s %s %s %s\n", 
		uc2b(binary[3], tmp_union_f.uc_val_f[3]),
		uc2b(binary[2], tmp_union_f.uc_val_f[2]),
		uc2b(binary[1], tmp_union_f.uc_val_f[1]),
		uc2b(binary[0], tmp_union_f.uc_val_f[0])
	);

}
// print float in binary format
void printh_float(float val)
{
	// union to print float data type
	union {
		unsigned char uc_val_f[4]; // 4 bytes = 32 bits
		float val_f;
	} tmp_union_f;

	tmp_union_f.val_f = val;

	// Big Endian
	printf("%02x %02x %02x %02x\n", 
		tmp_union_f.uc_val_f[3],
		tmp_union_f.uc_val_f[2],
		tmp_union_f.uc_val_f[1],
		tmp_union_f.uc_val_f[0]
	);
}

// print double in binary format
void printb_double(double val)
{
	int val_exp, implicit_bit;
	double val_significand;
	unsigned char binary[8][16];

	val_significand = frexp(val, &val_exp);
	implicit_bit = ((val_exp >= -1021) && (val_significand != 0)) ? '1' : '0';
	printf("2^(%+-d) * %25.17f\n", val_exp, val_significand);

	// union to print double data type
	union {
		unsigned char uc_val_d[8]; // 8 bytes = 64 bits
		double val_d;
	} tmp_union_d;

	tmp_union_d.val_d = val;

	// Big Endian
	uc2b(binary[7], tmp_union_d.uc_val_d[7]);
	uc2b(binary[6], tmp_union_d.uc_val_d[6]);
	uc2b(binary[5], tmp_union_d.uc_val_d[5]);
	uc2b(binary[4], tmp_union_d.uc_val_d[4]);
	uc2b(binary[3], tmp_union_d.uc_val_d[3]);
	uc2b(binary[2], tmp_union_d.uc_val_d[2]);
	uc2b(binary[1], tmp_union_d.uc_val_d[1]);
	uc2b(binary[0], tmp_union_d.uc_val_d[0]);

	// sign(1) exp(11) frac(52)+1
	printf("%c %c%c%c%c%c%c%c%c%c%c%c %c.%c%c%c%c%c %c%c%c%c%c %c%c%c%c%c %c%c%c%c%c %c%c%c%c%c %c%c%c%c%c %c%c%c%c%c %c%c%c%c%c %c%c%c%c%c %c%c%c%c%c %c%c\n",
		// sign
		binary[7][0], \
		// exponent
		binary[7][1], binary[7][2], binary[7][3], binary[7][4], binary[7][5], binary[7][6], binary[7][7], binary[6][0], binary[6][1], binary[6][2], binary[6][3], \
		 // significand
		implicit_bit, binary[6][4], binary[6][5], binary[6][6], binary[6][7], \
		binary[5][0], binary[5][1], binary[5][2], binary[5][3], binary[5][4], binary[5][5], binary[5][6], binary[5][7], \
		binary[4][0], binary[4][1], binary[4][2], binary[4][3], binary[4][4], binary[4][5], binary[4][6], binary[4][7], \
		binary[3][0], binary[3][1], binary[3][2], binary[3][3], binary[3][4], binary[3][5], binary[3][6], binary[3][7], \
		binary[2][0], binary[2][1], binary[2][2], binary[2][3], binary[2][4], binary[2][5], binary[2][6], binary[2][7], \
		binary[1][0], binary[1][1], binary[1][2], binary[1][3], binary[1][4], binary[1][5], binary[1][6], binary[1][7], \
		binary[0][0], binary[0][1], binary[0][2], binary[0][3], binary[0][4], binary[0][5], binary[0][6], binary[0][7]
	);

/*	printf("%s %s %s %s %s %s %s %s\n", 
		uc2b(binary[7], tmp_union_d.uc_val_d[7]),
		uc2b(binary[6], tmp_union_d.uc_val_d[6]),
		uc2b(binary[5], tmp_union_d.uc_val_d[5]),
		uc2b(binary[4], tmp_union_d.uc_val_d[4]),
		uc2b(binary[3], tmp_union_d.uc_val_d[3]),
		uc2b(binary[2], tmp_union_d.uc_val_d[2]),
		uc2b(binary[1], tmp_union_d.uc_val_d[1]),
		uc2b(binary[0], tmp_union_d.uc_val_d[0])
	);
*/
}
// print double in binary format
void printh_double(double val)
{
	// union to print double data type
	union {
		unsigned char uc_val_d[8]; // 8 bytes = 64 bits
		double val_d;
	} tmp_union_d;

	tmp_union_d.val_d = val;

	// Big Endian
	printf("%02x %02x %02x %02x %02x %02x %02x %02x\n", 
		tmp_union_d.uc_val_d[7],
		tmp_union_d.uc_val_d[6],
		tmp_union_d.uc_val_d[5],
		tmp_union_d.uc_val_d[4],
		tmp_union_d.uc_val_d[3],
		tmp_union_d.uc_val_d[2],
		tmp_union_d.uc_val_d[1],
		tmp_union_d.uc_val_d[0]
	);
}

// print long double in binary format
void printb_long_double(long double val)
{
	unsigned char binary[10][16];

	// union to print double data type
	union {
		unsigned char uc_val_d[10]; // 8 bytes = 64 bits
		double val_d;
	} tmp_union_d;

	tmp_union_d.val_d = val;

	// Big Endian
	printf("%s %s %s %s %s %s %s %s %s %s\n", 
		uc2b(binary[9], tmp_union_d.uc_val_d[9]),
		uc2b(binary[8], tmp_union_d.uc_val_d[8]),
		uc2b(binary[7], tmp_union_d.uc_val_d[7]),
		uc2b(binary[6], tmp_union_d.uc_val_d[6]),
		uc2b(binary[5], tmp_union_d.uc_val_d[5]),
		uc2b(binary[4], tmp_union_d.uc_val_d[4]),
		uc2b(binary[3], tmp_union_d.uc_val_d[3]),
		uc2b(binary[2], tmp_union_d.uc_val_d[2]),
		uc2b(binary[1], tmp_union_d.uc_val_d[1]),
		uc2b(binary[0], tmp_union_d.uc_val_d[0])
	);

}
} // extern "C" 


int main()
{
	int i, j;
	double a, b, w, s, e;
	unsigned int old_cw;

	// DD
	fpu_fix_start(&old_cw);

	a = 1.1111;
	b = 0.23333;

	//s = qd::quick_two_sum(a, b, e);
	s = a + b;
	w = s - a;
	e = b - w;

	printf("a = %25.17e, b = %25.17e\n", a, b);
	printb_double(a); printb_double(b);
	printf("s = %25.17e, w = %25.17e, e = %25.17e\n", s, w, e);
	printb_double(s); printb_double(w);printb_double(e);

	fpu_fix_end(&old_cw);

	return 0;

}
