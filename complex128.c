#include <stdio.h>
#include <quadmath.h> // GCC Quad-Math Library

int main(void)
{
	char input_str[128], output_str[2][128];

	// �����E��������double�^�Ƃ���
  	__complex128 a, b, c;
  	__float128 a_real, a_imag, b_real, b_imag;

	// a, b��W������(�L�[�{�[�h)����������
	printf("Input Re(a) ->"); scanf("%s", input_str); a_real = strtoflt128(input_str, NULL);
	printf("Input Im(a) ->"); scanf("%s", input_str); a_imag = strtoflt128(input_str, NULL);

	// a := (a_real) + (a_imag) * I
	__real__ a = a_real; __imag__ a = a_imag; // GCC�Ǝ�����

	printf("Input Re(b) ->"); scanf("%s", input_str); b_real = strtoflt128(input_str, NULL);
	printf("Input Im(b) ->"); scanf("%s", input_str); b_imag = strtoflt128(input_str, NULL);

	// b := (b_real) + (b_imag) * I
	__real__ b = b_real; __imag__ b = b_imag;

	// a, b���������C�������ɕ����ĕ\��
	quadmath_snprintf(output_str[0], sizeof(output_str[0]), "%+Qg", crealq(a));
	quadmath_snprintf(output_str[1], sizeof(output_str[1]), "%+Qg", cimagq(a));
	printf("a = (%s) + (%s) * I\n", output_str[0], output_str[1]);

	quadmath_snprintf(output_str[0], sizeof(output_str[0]), "%+Qg", crealq(b));
	quadmath_snprintf(output_str[1], sizeof(output_str[1]), "%+Qg", cimagq(b));
	printf("b = (%s) + (%s) * I\n", output_str[0], output_str[1]);

	// �W���o�͂Ɏl�����Z�̌��ʂ�\��
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

	// �I��
	return 0;
}
