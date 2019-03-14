#include <stdio.h>
#include <complex.h> // C99複素数型

int main(int argc, char* argv[])
{
	// 実部・虚部共にdouble型とする
  	double complex a, b, c;
  	double a_real, a_imag, b_real, b_imag;

	// a, bを標準入力(キーボード)から取り入れる
	printf("Input Re(a) ->"); scanf("%lf", &a_real);
	printf("Input Im(a) ->"); scanf("%lf", &a_imag);
	a = a_real + a_imag * I; // I = sqrt(-1)

	printf("Input Re(b) ->"); scanf("%lf", &b_real);
	printf("Input Im(b) ->"); scanf("%lf", &b_imag);
	b = b_real + b_imag * I;

	// a, bを実数部，虚数部に分けて表示
	printf("a = (%+g) + (%+g) * I\n", creal(a), cimag(a));
	printf("b = (%+g) + (%+g) * I\n", creal(b), cimag(b));

	// 標準出力に四則演算の結果を表示
	c = a + b;
	printf("a  + b = (%+g) + (%+g) * I\n", creal(c), cimag(c));

	c = a - b;
	printf("a  + b = (%+g) + (%+g) * I\n", creal(c), cimag(c));

	c = a * b;
	printf("a  + b = (%+g) + (%+g) * I\n", creal(c), cimag(c));

	c = a / b;
	printf("a  + b = (%+g) + (%+g) * I\n", creal(c), cimag(c));

	// 終了
	return 0;
}
