#include <iostream>
#include <iomanip>

#include "gmpxx.h"

/*
[tkouya@cs-muse mpna]$ g++ mpf_relerr.cc -lgmpxx -lgmp
[tkouya@cs-muse mpna]$ ./a.out
 1.4142135623730950488016887242096980785690000e+00 +  1.7320508075688772935274463415058723669420000e+00 =  3.1462643699419723423291350657155704455110000e+00
 1.4142135623730950488016887242096980785690000e+00 -  1.7320508075688772935274463415058723669420000e+00 = -3.1783724519578224472575761729617428837340000e-01
 1.4142135623730950488016887242096980785690000e+00 *  1.7320508075688772935274463415058723669420000e+00 =  2.4494897427831780981972840747058913919640000e+00
 1.4142135623730950488016887242096980785690000e+00 /  1.7320508075688772935274463415058723669420000e+00 =  8.1649658092772603273242802490196379732180000e-01
*/

using namespace std;

// 相対誤差
void relerr(mpf_class &c, mpf_class a, mpf_class ad)
{
	c = abs((a - ad) / ad);
}

int main(void)
{
	unsigned long prec; // 精度

	cout << " Input default prec in bits: "; cin >> prec;

	mpf_set_default_prec(prec); // in bits

	mpf_class a, b, c, d;
	mpf_class ad(0, prec * 2), bd(0, prec * 2), cd(0, prec * 2);

	//a = sqrt(2);
	//b = sqrt(3);
	//mpf_sqrt_ui(a.get_mpf_t(), 2UL);
	//mpf_sqrt_ui(b.get_mpf_t(), 3UL);
	mpf_sqrt_ui(a.get_mpf_t(), 5UL);
	mpf_sqrt_ui(b.get_mpf_t(), 7UL);

	mpf_sqrt_ui(ad.get_mpf_t(), 5UL);
	mpf_sqrt_ui(bd.get_mpf_t(), 7UL);

	//gmp_printf("ad = %97.90Fe\n", ad.get_mpf_t());
	relerr(c, a, ad);
	gmp_printf("relerr(a) = %10.3Fe\n", c.get_mpf_t());
	relerr(c, b, bd);
	gmp_printf("relerr(b) = %10.3Fe\n", c.get_mpf_t());

// unworking setprecision(xx) for mpf_class variables !
//	cout << setprecision(50) << a << " + " << b << " = " << a + b << endl;

	c = a + b; // short
	cd = ad + bd; // long
	relerr(c, c, cd);
	gmp_printf("relerr(a + b) = %10.3Fe\n", c.get_mpf_t());

	c = a - b; // short 
	cd = ad - bd; // long
	relerr(c, c, cd);
	gmp_printf("relerr(a - b) = %10.3Fe\n", c.get_mpf_t());

	c = a * b; // short
	cd = ad * bd; // long
	relerr(c, c, cd);
	gmp_printf("relerr(a * b) = %10.3Fe\n", c.get_mpf_t());

	c = a / b; // short
	cd = ad / bd; // long
	relerr(c, c, cd);
	gmp_printf("relerr(a / b) = %10.3Fe\n", c.get_mpf_t());

	return 0;
}
