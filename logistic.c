// logistic Ê‘œ
#include <stdio.h>

int main()
{
	int i;
	double x[102];

	// ‰Šú’l
	x[0] = 0.7501;

	for(i = 0; i <= 100; i++)
	{
		x[i + 1] = 4 * x[i] * (1 - x[i]);
		if((i % 10) == 0)
			printf("%5d, %25.17e\n", i, x[i]);
	}

	return 0;
}
