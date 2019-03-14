// logistic Ê‘œ
#include <iostream>
#include <iomanip>
#include "qd_real.h"

using namespace std;

int main()
{
	int i;
	dd_real x[102]; // DD¸“x

	// ‰Šú’l
	x[0] = "0.7501";

	fpu_fix_start(NULL);

	for(i = 0; i <= 100; i++)
	{
		x[i + 1] = 4 * x[i] * (1 - x[i]);
		if((i % 10) == 0)
			cout << setw(5) << i << setprecision(32) << ", " << x[i] << endl;
	}

	return 0;
}
