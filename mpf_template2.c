// mpf_template.c: mpf_t�^�̎g����

#include <stdio.h>
#include <math.h>
#include "gmp.h"

int main(void)
{
	unsigned long prec, dec_prec;
	char format[128] = "";
	mpf_t a, b, c; // �f�t�H���g�̐��x
	mpf_t ad, bd, cd; // a, b, c��蒷����

	printf("Input default prec in bits: "); scanf("%ld", &prec);

	// �f�t�H���g�̉������r�b�g�����Z�b�g
	mpf_set_default_prec((mp_bitcnt_t)prec); // in bits

	// �f�t�H���g�̐��x��\��
	printf("default prec in bits: %ld\n", prec);

	mpf_init(a); // �f�t�H���g���x
	mpf_inits(b, c, NULL); // �f�t�H���g���x�ł܂Ƃ߂ď�����

	// �f�t�H���g��������2�{�ɐݒ�
	mpf_init2(ad, (mp_bitcnt_t)(prec * 2));
	mpf_init2(bd, (mp_bitcnt_t)(prec * 2));
	mpf_init2(cd, (mp_bitcnt_t)(prec * 2));

	// �ϐ����Ƃ̉������r�b�g����\��
	printf("prec(ad) = %ld\n", mpf_get_prec(ad));

	// �l���Z�b�g���ĕ\��
	mpf_set_ui(a, 5UL);
	mpf_set_d(b, sqrt(2.0));
	mpf_set_str(c, "3.1415926535897932384626433832795e+1568759", 10);
	gmp_printf("a = %Fe, b = %Fe\n", a, b);
	gmp_printf("c = %50.43Fe\n", c);

	// �W���^�ɕϊ����ĕ\��
	printf("a = %ld\n", mpf_get_ui(a));
	printf("b = %f\n", mpf_get_d(b));
	printf("c = %25.17e\n", mpf_get_d(c));

	// �t���\��
	dec_prec = (unsigned long)ceil(log10(2.0) * prec);
	sprintf(format, "%%%ld.%ldFe", dec_prec + 8, dec_prec);

	mpf_set_ui(ad, 2UL); mpf_sqrt(ad, ad);
	printf("ad = "); gmp_printf(format, ad); printf("\n");

	// �ϐ�����
	mpf_clear(a); // �P��

	// �܂Ƃ߂ď���
	mpf_clears(b, c, NULL);       
	mpf_clears(ad, bd, cd, NULL); 

	return 0;
}
