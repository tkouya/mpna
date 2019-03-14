#include <stdio.h>
#include <quadmath.h> // GCC Quad-Math Library

int main(void)
{
	char input_str[128], output_str[128];

  	__float128 a, b, c;

	// a, b��W������(�L�[�{�[�h)����������
	printf("Input a ->"); scanf("%s", input_str); a = strtoflt128(input_str, NULL);
	printf("Input b ->"); scanf("%s", input_str); b = strtoflt128(input_str, NULL);

	// a, b��\��
	quadmath_snprintf(output_str, sizeof(output_str), "%+Qg", a);
	printf("a = %s\n", output_str);
	quadmath_snprintf(output_str, sizeof(output_str), "%+Qa", a);
	printf("a = %s\n", output_str);

	quadmath_snprintf(output_str, sizeof(output_str), "%+Qg", b);
	printf("b = %s\n", output_str);

	// �W���o�͂Ɏl�����Z�̌��ʂ�\��
	c = a + b;
	quadmath_snprintf(output_str, sizeof(output_str), "%+.34Qg", c);
	printf("a + b = %s\n", output_str);

	c = a - b;
	quadmath_snprintf(output_str, sizeof(output_str), "%+.34Qg", c);
	printf("a - b = %s\n", output_str);

	c = a * b;
	quadmath_snprintf(output_str, sizeof(output_str), "%+.34Qg", c);
	printf("a * b = %s\n", output_str);

	c = a / b;
	quadmath_snprintf(output_str, sizeof(output_str), "%+.34Qg", c);
	printf("a / b = %s\n", output_str);

	c = sqrtq(a);
	quadmath_snprintf(output_str, sizeof(output_str), "%+.34Qg", c);
	printf("sqrt(a) = %s\n", output_str);

	c = c * c;
	quadmath_snprintf(output_str, sizeof(output_str), "%+.34Qg", c);
	printf("(sqrt(a))^2 = %s\n", output_str);

	// �I��
	return 0;
}
