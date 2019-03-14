#include <stdio.h>

#include "mpfr.h"

// Emulation program of IEEE754 double precision arithmetic
// ref. MPFR 4.0.2 Mamual

int main(void)
{
	mpfr_t xa, xb; int i;
	float fa, fb;
	double a, b;

// �P���x�G�~�����[�g

	// ������ 24 bits
	mpfr_set_default_prec(24);

	// �A���_�[�t���[�̋��E�w����: 2^(-125-23)
	mpfr_set_emin(-148);

	// �I�[�o�[�t���[�̋��E�w����: 2^(+126)
	mpfr_set_emax(126);

	// �f�t�H���g���x���ŏ�����
	mpfr_init(xa); mpfr_init(xb);

	// �P���x�l���Z�b�g
	fb = 34.3f; mpfr_set_d(xb, (double)fb, MPFR_RNDN);

	// a = hex(1.1235)�~2^(-125)
 	fa = 0x1.1235P-125; mpfr_set_d(xa, (double)fa, MPFR_RNDN);

	// �P���x�v�Z
	fa /= fb;
	printf("float: a / b = %15.7e\n", fa); 

	// MPFR�G�~�����[�g�{���x�v�Z
	i = mpfr_div(xa, xa, xb, MPFR_RNDN);
	mpfr_printf("MPFR emulated float: xa / xb = %15.7RNe\n", xa);
	mpfr_printf("                             = %15.6RNa\n", xa);
	i = mpfr_subnormalize (xa, i, MPFR_RNDN); /* new ternary value */
	mpfr_printf("              subnormalized -> %15.7RNe\n", xa);
	mpfr_printf("                             = %15.6RNa\n", xa);

	// ����
	mpfr_clear(xa); mpfr_clear(xb);


// �{���x�G�~�����[�g

	// ������ 53 bits
	mpfr_set_default_prec(53);

	// �A���_�[�t���[�̋��E�w����: 2^(-1021-52)
	mpfr_set_emin(-1073);

	// �I�[�o�[�t���[�̋��E�w����: 2^(+1024)
	mpfr_set_emax(1024);

	// �f�t�H���g���x���ŏ�����
	mpfr_init(xa); mpfr_init(xb);

	// �{���x�l���Z�b�g
	b = 34.3; mpfr_set_d(xb, b, MPFR_RNDN);

	// a = hex(1.1235)�~2^(-1021)
 	a = 0x1.1235P-1021; mpfr_set_d(xa, a, MPFR_RNDN);

	// �{���x�v�Z
	a /= b;
	printf("double: a / b = %25.17e\n",a); 

	// MPFR�G�~�����[�g�{���x�v�Z
	i = mpfr_div(xa, xa, xb, MPFR_RNDN);
	mpfr_printf("MPFR emulated double: xa / xb  = %25.17RNe\n", xa);
	mpfr_printf("                               = %25.12RNa\n", xa);
	i = mpfr_subnormalize (xa, i, MPFR_RNDN); /* new ternary value */
	mpfr_printf("                subnormalized -> %25.17RNe\n", xa);
	mpfr_printf("                               = %25.12RNa\n", xa);

	// ����
	mpfr_clear(xa); mpfr_clear(xb);

	return 0;
}
