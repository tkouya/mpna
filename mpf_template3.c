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

	//a = sqrt(2);
	//b = sqrt(3);
	mpf_set_ui(a, 2UL); mpf_sqrt(a, a);
	mpf_set_ui(b, 3UL); mpf_sqrt(b, b);
	
	mpf_add(c, a, b); // c = a + b;
	gmp_printf("%50.43Fe + %50.43Fe = %50.43Fe\n", a, b, c);
	
	mpf_sub(c, a, b); // c = a - b;
	gmp_printf("%50.43Fe - %50.43Fe = %50.43Fe\n", a, b, c);
	
	mpf_mul(c, a, b); // c = a * b;
	gmp_printf("%50.43Fe * %50.43Fe = %50.43Fe\n", a, b, c);
	
	mpf_div(c, a, b); // c = a / b;
	gmp_printf("%50.43Fe / %50.43Fe = %50.43Fe\n", a, b, c);

	mpf_sqrt(c, a); // c = sqrt(a);
	gmp_printf("sqrt(%50.43Fe) = %50.43Fe\n", a, c);

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
