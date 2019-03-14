#include <stdio.h>
#include <complex.h> // C99���f���^

int main(int argc, char* argv[])
{
	// �����E��������double�^�Ƃ���
  	double complex a, b, c;
  	double a_real, a_imag, b_real, b_imag;

	// a, b��W������(�L�[�{�[�h)����������
	printf("Input Re(a) ->"); scanf("%lf", &a_real);
	printf("Input Im(a) ->"); scanf("%lf", &a_imag);
	a = a_real + a_imag * I; // I = sqrt(-1)

	printf("Input Re(b) ->"); scanf("%lf", &b_real);
	printf("Input Im(b) ->"); scanf("%lf", &b_imag);
	b = b_real + b_imag * I;

	// a, b���������C�������ɕ����ĕ\��
	printf("a = (%+g) + (%+g) * I\n", creal(a), cimag(a));
	printf("b = (%+g) + (%+g) * I\n", creal(b), cimag(b));

	// �W���o�͂Ɏl�����Z�̌��ʂ�\��
	c = a + b;
	printf("a  + b = (%+g) + (%+g) * I\n", creal(c), cimag(c));

	c = a - b;
	printf("a  + b = (%+g) + (%+g) * I\n", creal(c), cimag(c));

	c = a * b;
	printf("a  + b = (%+g) + (%+g) * I\n", creal(c), cimag(c));

	c = a / b;
	printf("a  + b = (%+g) + (%+g) * I\n", creal(c), cimag(c));

	// �I��
	return 0;
}
