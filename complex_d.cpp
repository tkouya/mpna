#include <iostream>
#include <complex>

// 名前空間はstdを使用
using namespace std;

int main(int argc, char* argv[])
{
	// 実部・虚部共にdouble型とする
  	complex<double> a, b;

	// a, bを標準入力(キーボード)から取り入れる
	// 入力時は"(実数部,虚数部)"と指定すること！
	cout << "Input a ->";
	cin >> a;
	cout << "Input b ->";
  	cin >> b;

	// a, bを実数部，虚数部に分けて表示
	cout << "a = " << a.real() << " + " << a.imag() << " * I" << endl;
	cout << "b = " << b.real() << " + " << b.imag() << " * I" << endl;

	// 標準出力に四則演算の結果を表示
	cout << a << " + " << b << " = " << a + b << endl;
  	cout << a << " - " << b << " = " << a - b << endl;
  	cout << a << " * " << b << " = " << a * b << endl;
  	cout << a << " / " << b << " = " << a / b << endl;

	// 終了
	return 0;
}