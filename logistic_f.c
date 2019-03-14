// logistic Ê‘œ
#include <stdio.h>

int main()
{
	int i;
	float x[102];

	// ‰Šú’l
	x[0] = 0.7501;

	for(i = 0; i <= 100; i++)
	{
		x[i + 1] = 4 * x[i] * (1 - x[i]);
		if((i % 10) == 0)
			printf("%5d, %15.7e\n", i, x[i]);
	}

	return 0;
}
