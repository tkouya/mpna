#include <iostream>
#include <complex>

// ���O��Ԃ�std���g�p
using namespace std;

int main(int argc, char* argv[])
{
	// �����E��������double�^�Ƃ���
  	complex<double> a, b;

	// a, b��W������(�L�[�{�[�h)����������
	// ���͎���"(������,������)"�Ǝw�肷�邱�ƁI
	cout << "Input a ->";
	cin >> a;
	cout << "Input b ->";
  	cin >> b;

	// a, b���������C�������ɕ����ĕ\��
	cout << "a = " << a.real() << " + " << a.imag() << " * I" << endl;
	cout << "b = " << b.real() << " + " << b.imag() << " * I" << endl;

	// �W���o�͂Ɏl�����Z�̌��ʂ�\��
	cout << a << " + " << b << " = " << a + b << endl;
  	cout << a << " - " << b << " = " << a - b << endl;
  	cout << a << " * " << b << " = " << a * b << endl;
  	cout << a << " / " << b << " = " << a / b << endl;

	// �I��
	return 0;
}