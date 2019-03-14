#include <stdio.h>
#include <math.h>
#include <float.h>
#include <fenv.h>

// View current rounding mode
void view_rnmode(void)
{
	// rounding mode
	switch(fegetround())
	{
		case FE_TONEAREST:
			printf("RN mode\n"); break;

		case FE_TOWARDZERO:
			printf("RZ mode\n"); break;

		case FE_UPWARD:
			printf("RP mode\n"); break;

		case FE_DOWNWARD:
			printf("RM mode\n"); break;
	}
}

int main(void)
{

	// Float
	printf("FLT_MIN     = %+15.7e\n", FLT_MIN);
	printf("FLT_MAX     = %+15.7e\n", FLT_MAX);
	printf("FLT_EPSILON = %+15.7e\n", FLT_EPSILON);

	// Double
	printf("DBL_MIN     = %+25.17e\n", DBL_MIN);
	printf("DBL_MAX     = %+25.17e\n", DBL_MAX);
	printf("DBL_EPSILON = %+25.17e\n", DBL_EPSILON);

#pragma STDC FENV_ACCESS ON
	// rounding mode
	fesetround(FE_TONEAREST);
	view_rnmode();
	printf("RN(sqrt(2) * sqrt(3)) = %25.17e\n", acosf(-1));
	printf("RN(sqrt(2) * sqrt(3)) = %25.17e\n", acos(-1));

	fesetround(FE_TOWARDZERO);
	view_rnmode();
	printf("RZ(sqrt(2) * sqrt(3)) = %25.17e\n", acosf(-1));
	printf("RZ(sqrt(2) * sqrt(3)) = %25.17e\n", acos(-1));

	fesetround(FE_UPWARD);
	view_rnmode();
	printf("RP(sqrt(2) * sqrt(3)) = %25.17e\n", acosf(-1));
	printf("RP(sqrt(2) * sqrt(3)) = %25.17e\n", acos(-1));

	fesetround(FE_DOWNWARD);
	view_rnmode();
	printf("RM(sqrt(2) * sqrt(3)) = %25.17e\n", acosf(-1));
	printf("RM(sqrt(2) * sqrt(3)) = %25.17e\n", acos(-1));

	return 0;
}
