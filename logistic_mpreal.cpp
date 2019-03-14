// logistic 写像
#include <iostream>
#include <iomanip>

// Multiple precision with MPFR/GMP
#include "mpreal.h"

using namespace std;
using namespace mpfr;

int main(int argc, char *argv[])
{
	int i;
	int num_bits, num_decimal;

	// 引数チェック
	if(argc <= 1)
	{
		cerr << "Usage: " << argv[0] << " [num_bits]" << endl;
		return 0;
	}

	// 計算桁数設定
	num_bits = atoi(argv[1]);
	if(num_bits <= 24)
		num_bits = 24;

	num_decimal = (int)ceil(log10(2.0) * (double)num_bits);
	mpreal::set_default_prec(num_bits);

	cout << "num_bits = " << num_bits << ", num_decimal = " << num_decimal << endl;

	mpreal x[102];

	// 初期値
	x[0] = "0.7501";

	for(i = 0; i <= 100; i++)
	{
		x[i + 1] = 4 * x[i] * (1 - x[i]);
		if((i % 10) == 0)
			cout << setw(5) << scientific << i << setprecision(num_decimal) << ", " << x[i] << endl;
	}

	return 0;
}
