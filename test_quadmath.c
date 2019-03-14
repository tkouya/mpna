#include <stdio.h>
#include <quadmath.h>

int main()
{
	__float128 a, b, c;
	char buf[1024];

	a = sqrtq(2.0q);
	b = sqrtq(3.0q);

	quadmath_snprintf(buf, 1024, "%50.33Qe", a);
	printf("a = %s\n", buf);
	quadmath_snprintf(buf, 1024, "%50.33Qe", b);
	printf("b = %s\n", buf);

	c = a + b;

	quadmath_snprintf(buf, 1024, "%Qe", c);
	printf("c = %s\n", buf);

	return 0;
}

