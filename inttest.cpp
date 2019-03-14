#include <iostream>
#include <iomanip>

using namespace std;

/*
[tkouya@cs-muse mpna]$ ./a.out
a = 2147500032
b = 2147483648
2147500032 + 2147483648 = 16384
2147500032 - 2147483648 = 16384
*/

int main()
{
	unsigned int a, b;
	int ia, ib;

	cout << "a = ";	cin >> a;
	cout << "b = "; cin >> b;

	cout << a << " + " << b << " = " << (unsigned int)(a + b) << endl;
	cout << a << " - " << b << " = " << (unsigned int)(a - b) << endl;
	cout << a << " * " << b << " = " << (unsigned int)(a * b) << endl;

	ia = a;
	ib = b;

	cout << ia << " + " << ib << " = " << (int)(ia + ib) << endl;
	cout << ia << " - " << ib << " = " << (int)(ia - ib) << endl;
	cout << ia << " * " << ib << " = " << (int)(ia + ib) << endl;
	cout << hex << showbase << ia << " + " << ib << " = " << (int)(ia + ib) << endl;
	cout << hex << showbase << ia << " - " << ib << " = " << (int)(ia - ib) << endl;
	cout << hex << showbase << ia << " * " << ib << " = " << (int)(ia * ib) << endl;

	return 0;
}
