//******************************************************************************
// rsa.c : RSA Encryption by GNU MP
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
#include <string.h>

#include "gmp.h"

/* Number of character */
//#define BASE 128 // ASCII code
#define BASE 256 // JIS X 201

/* Max length of string */
#ifndef MAXLEN
  #define MAXLEN 256
#endif

/* string to integer */
void str2int(mpz_t ret, char *str)
{
	long int str_len, i;
	unsigned char c;

	/* chop ret code */
	str_len = strlen(str);
	if(str[str_len - 1] == '\n')
		str[str_len - 1] = '\0';

	str_len = strlen(str);

//	printf("%s -> %d \n", str, str_len);

	/* ret = 0 */
	mpz_set_ui(ret, 0UL);

	/* ret = str[str_len - 1] * BASE^(str_len - 1) + ... + str[1] * BASE + str[0] */
	for(i = str_len - 1; i >= 0; i--)
	{
		c = str[i];

		/* Horner's method */
		mpz_mul_ui(ret, ret, (unsigned long)BASE);
		mpz_add_ui(ret, ret, (unsigned long)c);
	}
}

/* integer to string */
void int2str(char *str, mpz_t org_str_int)
{
	long int str_len, i;
	mpz_t max_int, c_int, str_int;

	/* Initialize */
	mpz_init(max_int);
	mpz_init(c_int);
	mpz_init(str_int);

	mpz_set(str_int, org_str_int);

	/* get length of string */
	mpz_set_ui(max_int, 1UL);
	for(i = 0; i < MAXLEN; i++)
	{
//		gmp_printf("%Zd <= %Zd ?\n", str_int, max_int);
		if(mpz_cmp(str_int, max_int) <= 0)
		{
			str_len = i;
			break;
		}

		mpz_mul_ui(max_int, max_int, (unsigned long)BASE);
	}

	/* convert to letter */
	for(i = 0; i < str_len; i++)
	{
		mpz_mod_ui(c_int, str_int, (unsigned long)BASE);
		mpz_sub(str_int, str_int, c_int);
		mpz_tdiv_q_ui(str_int, str_int, (unsigned long)BASE);

		str[i] = mpz_get_ui(c_int);
	}
	str[str_len] = '\0';

//	printf("Length-> %d\n", str_len);

	/* free */
	mpz_clear(max_int);
	mpz_clear(c_int);
	mpz_clear(str_int);
}

/* get public and private keys */
void get_public_private_key(mpz_t n, mpz_t e, mpz_t d, mpz_t org_str_int)
{
	mpz_t p, q, str_int, euler_n, tmp; // prime numbers

	/* Initialize */
	mpz_init(p);
	mpz_init(q);
	mpz_init(str_int);
	mpz_init(euler_n);
	mpz_init(tmp);

	/* tmp = str_int * BASE */
	mpz_mul_ui(str_int, org_str_int, (unsigned long)BASE);

	/* get two prime numbers */
	mpz_nextprime(q, str_int);
	mpz_nextprime(p, q); // p > q > str_int * BASE

	/* n = p * q */
	mpz_mul(n, p, q);

	/* euler(n) */
	mpz_sub_ui(p, p, 1UL);
	mpz_sub_ui(q, q, 1UL);
	mpz_mul(euler_n, p, q); // euler(n) = (p-1)*(q-1)

	/* select e ... very small and not secure!! */
	mpz_set_ui(e, 2UL);
	while(1) {
		mpz_gcd(tmp, euler_n, e);
		if(mpz_cmp_ui(tmp, 1UL) == 0)
			break;
		mpz_nextprime(e, e);
	}

	/* get private key */
	mpz_gcdext(p, d, q, e, euler_n);
	//mpz_gcdext(p, q, d, euler_n, e);
	if(mpz_sgn(d) < 0)
		mpz_add(d, d, euler_n);

	/* Free */
	mpz_clear(p);
	mpz_clear(q);
	mpz_clear(str_int);
	mpz_clear(euler_n);
	mpz_clear(tmp);
}

/* Encryption */
void encrypt(mpz_t enc_str_int, mpz_t str_int, mpz_t public_key_n, mpz_t public_key_e)
{
	/* enc_str_int = str_int^public_key_e mod public_key_n */
	mpz_powm(enc_str_int, str_int, public_key_e, public_key_n);
}

/* Decryption */
void decrypt(mpz_t str_int, mpz_t enc_str_int, mpz_t private_key, mpz_t public_key_n)
{
	/* str_int = enc_str_int^private_key mod public_key_n */
	mpz_powm(str_int, enc_str_int, private_key, public_key_n);
}

int main()
{
	unsigned char string[MAXLEN], out_string[MAXLEN];
	mpz_t string_int, string_int_encrypt;
	mpz_t public_key[2], private_key;

	/* Initialize */
	mpz_init(string_int);
	mpz_init(string_int_encrypt);
	mpz_init(public_key[0]);
	mpz_init(public_key[1]);
	mpz_init(private_key);

	/* Input string */
	printf("Input a string(Max length = %d): ", MAXLEN);
	fgets(string, MAXLEN - 1, stdin);

	/* string -> integer */
	str2int(string_int, string);

#ifdef HEX
	gmp_printf("string(%d)  -> %Zx\n\n", strlen(string), string_int);
#else
	gmp_printf("string(%d)  -> %Zd\n\n", strlen(string), string_int);
#endif

	/* get public keys */
	get_public_private_key(public_key[0], public_key[1], private_key, string_int);

#ifdef HEX
	gmp_printf("Public_key_n : %Zx\n", public_key[0]);
	gmp_printf("Public_key_e : %Zx\n", public_key[1]);
	gmp_printf("Private_key_d: %Zx\n\n", private_key);
#else
	gmp_printf("Public_key_n : %Zd\n", public_key[0]);
	gmp_printf("Public_key_e : %Zd\n", public_key[1]);
	gmp_printf("Private_key_d: %Zd\n\n", private_key);
#endif

	/* Encrypt */
	encrypt(string_int_encrypt, string_int, public_key[0], public_key[1]);

#ifdef HEX
	gmp_printf("Encryption   : %Zx\n", string_int_encrypt);
#else
	gmp_printf("Encryption   : %Zd\n", string_int_encrypt);
#endif

	/* Decode */
	decrypt(string_int, string_int_encrypt, private_key, public_key[0]);

#ifdef HEX
	gmp_printf("Decrypt      : %Zx\n\n", string_int);
#else
	gmp_printf("Decrypt      : %Zd\n\n", string_int);
#endif

	/* integer -> string */
	int2str(out_string, string_int);

#ifdef HEX
	gmp_printf("%Zd\n            -> %s(%d)\n", string_int, out_string, strlen(out_string));
#else
	gmp_printf("%Zd\n            -> %s(%d)\n", string_int, out_string, strlen(out_string));
#endif
	/* free */
	mpz_clear(string_int);
	mpz_clear(string_int_encrypt);
	mpz_clear(public_key[0]);
	mpz_clear(public_key[1]);
	mpz_clear(private_key);

	return 0;
}
