#include <stdio.h>
#include <quadmath.h> // GCC Quad-Math Library

int main(void)
{
	char input_str[128], output_str[2][128];

	// 実部・虚部共にdouble型とする
  	__complex128 a, b, c;
  	__float128 a_real, a_imag, b_real, b_imag;

	// a, bを標準入力(キーボード)から取り入れる
	printf("Input Re(a) ->"); scanf("%s", input_str); a_real = strtoflt128(input_str, NULL);
	printf("Input Im(a) ->"); scanf("%s", input_str); a_imag = strtoflt128(input_str, NULL);

	// a := (a_real) + (a_imag) * I
	__real__ a = a_real; __imag__ a = a_imag; // GCC独自命令

	printf("Input Re(b) ->"); scanf("%s", input_str); b_real = strtoflt128(input_str, NULL);
	printf("Input Im(b) ->"); scanf("%s", input_str); b_imag = strtoflt128(input_str, NULL);

	// b := (b_real) + (b_imag) * I
	__real__ b = b_real; __imag__ b = b_imag;

	// a, bを実数部，虚数部に分けて表示
	quadmath_snprintf(output_str[0], sizeof(output_str[0]), "%+Qg", crealq(a));
	quadmath_snprintf(output_str[1], sizeof(output_str[1]), "%+Qg", cimagq(a));
	printf("a = (%s) + (%s) * I\n", output_str[0], output_str[1]);

	quadmath_snprintf(output_str[0], sizeof(output_str[0]), "%+Qg", crealq(b));
	quadmath_snprintf(output_str[1], sizeof(output_str[1]), "%+Qg", cimagq(b));
	printf("b = (%s) + (%s) * I\n", output_str[0], output_str[1]);

	// 標準出力に四則演算の結果を表示
	c = a + b;
	quadmath_snprintf(output_str[0], sizeof(output_str[0]), "%+.34Qg", crealq(c));
	quadmath_snprintf(output_str[1], sizeof(output_str[1]), "%+.34Qg", cimagq(c));
	printf("a + b = (%s) + (%s) * I\n", output_str[0], output_str[1]);

	c = a - b;
	quadmath_snprintf(output_str[0], sizeof(output_str[0]), "%+.34Qg", crealq(c));
	quadmath_snprintf(output_str[1], sizeof(output_str[1]), "%+.34Qg", cimagq(c));
	printf("a - b = (%s) + (%s) * I\n", output_str[0], output_str[1]);

	c = a * b;
	quadmath_snprintf(output_str[0], sizeof(output_str[0]), "%+.34Qg", crealq(c));
	quadmath_snprintf(output_str[1], sizeof(output_str[1]), "%+.34Qg", cimagq(c));
	printf("a * b = (%s) + (%s) * I\n", output_str[0], output_str[1]);

	c = a / b;
	quadmath_snprintf(output_str[0], sizeof(output_str[0]), "%+.34Qg", crealq(c));
	quadmath_snprintf(output_str[1], sizeof(output_str[1]), "%+.34Qg", cimagq(c));
	printf("a / b = (%s) + (%s) * I\n", output_str[0], output_str[1]);

	c = csqrtq(a);
	quadmath_snprintf(output_str[0], sizeof(output_str[0]), "%+.34Qg", crealq(c));
	quadmath_snprintf(output_str[1], sizeof(output_str[1]), "%+.34Qg", cimagq(c));
	printf("csqrt(a) = (%s) + (%s) * I\n", output_str[0], output_str[1]);

	c = c * c;
	quadmath_snprintf(output_str[0], sizeof(output_str[0]), "%+.34Qg", crealq(c));
	quadmath_snprintf(output_str[1], sizeof(output_str[1]), "%+.34Qg", cimagq(c));
	printf("(csqrt(a))^2 = (%s) + (%s) * I\n", output_str[0], output_str[1]);

	// 終了
	return 0;
}
