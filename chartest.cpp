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
	unsigned char a, b;
	char ia, ib;

	cout << "a = ";	cin >> a;
	cout << "b = "; cin >> b;

	cout << a << " + " << b << " = " << (unsigned char)(a + b) << endl;
	cout << a << " - " << b << " = " << (unsigned char)(a - b) << endl;

	ia = a;
	ib = b;

	cout <<  ia << " + " << ib << " = " << (char)(ia + ib) << endl;
	cout << ia << " - " << ib << " = " << (char)(ia - ib) << endl;
	cout << hex << showbase << ia << " + " << ib << " = " << (char)(ia + ib) << endl;
	cout << hex << showbase << ia << " - " << ib << " = " << (char)(ia - ib) << endl;

	return 0;
}
