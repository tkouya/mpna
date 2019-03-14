// mpf_template.c: mpf_t�^�̎g����

#include <stdio.h>
#include "gmp.h"

int main(void)
{
	unsigned long prec;
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

	// �E�E�E�E�E�E�E�E�E�E
	// �E�E�E���Z�����E�E�E
	// �E�E�E�E�E�E�E�E�E�E

	// �ϐ�����
	mpf_clear(a); // �P��

	// �܂Ƃ߂ď���
	mpf_clears(b, c, NULL);       
	mpf_clears(ad, bd, cd, NULL); 

	return 0;
}
