#include <stdio.h>
#include <stdlib.h>

int main()
{
	int a, index;
	long int la;
	unsigned int ua;
	unsigned long int ula;

	// integer
	printf("     int: ");
	for(a = 1, index = 0; ; a *= 2, index++)
	{
		printf("%d ", a);
		if(a == 0) break;
	}
	printf("= 2^%d\n", --index);

	// long integer
	printf("long int: ");
	for(la = 1, index = 0; ; la *= 2, index++)
	{
		printf("%ld ", la);
		if(la == 0) break;
	}
	printf("= 2^%d\n", --index);

	// unsigned integer
	printf("unsigned int: ");
	for(ua = 1, index = 0; ; ua *= 2, index++)
	{
		printf("%u ", ua);
		if(ua == 0) break;
	}
	printf("= 2^%d\n", --index);

	// unsigned long integer
	printf("unsigned long int: ");
	for(ula = 1, index = 0; ; ula *= 2, index++)
	{
		printf("%lu ", ula);
		if(ula == 0) break;
	}
	printf("= 2^%d\n", --index);

	return EXIT_SUCCESS;
}