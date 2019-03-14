#include <iostream>
#include <iomanip>

using namespace std;

int main()
{
	int i;
	unsigned int uia, uiam1, uiap1;
	int ia, iap1, iam1;
	unsigned long int ula, ulap1, ulam1;
	long int la, lap1, lam1;

	ia = 1;
	la = 1;

	for(i = 0; i < 64; i++)
	{
		ia *= 2; la *= 2; uia = ia; ula = la;
		iap1 = ia + 1; lap1 = la + 1; uiap1 = iap1; ulap1 = lap1;
		iam1 = ia - 1; lam1 = la - 1; uiam1 = iam1; ulam1 = lam1;

		cout << "2^" <<  i - 1 << " - 1 = " << iam1 << ", " << uiam1 << ", " << lam1 << ", " << ulam1 << endl;
		cout << "2^" <<  i     << "     = " << ia   << ", " << uia   << ", " << la   << ", " << ula   << endl;
		cout << "2^" <<  i + 1 << " + 1 = " << iap1 << ", " << uiap1 << ", " << lap1 << ", " << ulap1 << endl;
	}

	return 0;
}
