#include <stdio.h>
#include <math.h>

int main(void)
{
	float a, b, c;
	double da, db, dc;

	a = NAN;
	b = +INFINITY;
	c = -INFINITY;

	da = NAN;
	db = +INFINITY;
	dc = -INFINITY;

	printf("a, b, c = %f, %f, %f\n", a, b, c);
	printf("da, db, dc = %lf, %lf, %lf\n", da, db, dc);

	printf("(+INF) * (-INF) = %f, %lf\n", b * c, db * dc);
 	printf("(-INF) * (-INF) = %f, %lf\n", c * c, dc * dc);
	printf("(+INF) / (-INF) = %f, %lf\n", b / c, db / dc);

	a = +0.0;
	b = -0.0;

	printf("zeros = %+f, %+f\n", a, b);

	return 0;

}