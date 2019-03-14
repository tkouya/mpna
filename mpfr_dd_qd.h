#ifndef __MPFR_DD_QD_H__
#define __MPFR_DD_QD_H__

// ---------------------------------------------------------------------------------
// mpfr_dd_qd.h
// These following codes are originally written by members of MPFR development team.
// Modified by Tomonori Kouya
// All codes are used under MPFR's license.
// ---------------------------------------------------------------------------------

extern "C" {

// based on get_ld.c code in MPFR 4.0.0 
#define _IN_H_MPFR_EXP_ZERO (mpfr_get_emin() + 1)
#define _IN_H_MPFR_EXP_NAN  (mpfr_get_emin() + 2)
#define _IN_H_MPFR_EXP_INF  (mpfr_get_emin() + 3)
#define _IN_H_MPFR_EXP_UBF  (mpfr_get_emin() + 4)

#define _IN_H_MPFR_SIGN(x) ((x)->_mpfr_sign)
#define _IN_H_MPFR_SIGN_POS (1)
#define _IN_H_MPFR_SIGN_NEG (-1)

#define _IN_H_MPFR_IS_NEG(x) (MPFR_SIGN(x) < 0)
#define _IN_H_MPFR_IS_POS(x) (MPFR_SIGN(x) > 0)

#define _IN_H_MPFR_SET_POS(x) (MPFR_SIGN(x) = MPFR_SIGN_POS)
#define _IN_H_MPFR_SET_NEG(x) (MPFR_SIGN(x) = MPFR_SIGN_NEG)

// FALSE(=0), TRUE(=1)
//#define MPFR_IS_SINGULAR(x) (MPFR_EXP(x) <= MPFR_EXP_INF)
int mpfr_is_singular(mpfr_srcptr x)
{
	return (mpfr_get_exp(x) <= _IN_H_MPFR_EXP_INF);
}

/* generic code */
// ret_dd[2] = ret_dd[high == 0], ret_dd[low == 1]
void mpfr_get_dd(double ret_dd[2], mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
	//if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (x)))
	if (mpfr_is_singular(x))
	{
		//return (long double) mpfr_get_d (x, rnd_mode);
		ret_dd[0] = mpfr_get_d(x, rnd_mode);
		ret_dd[1] = (double)0.0;
	}
	else /* now x is a normal non-zero number */
	{
		// long double r; /* result */
		double s; /* part of result */

		//MPFR_SAVE_EXPO_DECL (expo);
		//MPFR_SAVE_EXPO_MARK (expo);

		//#if defined(HAVE_LDOUBLE_MAYBE_DOUBLE_DOUBLE)
		//if (MPFR_LDBL_MANT_DIG == 106)
		//{
		/* Assume double-double format (as found with the PowerPC ABI).
		   The generic code below isn't used because numbers with
		   precision > 106 would not be supported. */
		s = mpfr_get_d (x, MPFR_RNDN); /* high part of x */
		/* Let's first consider special cases separately. The test for
		   infinity is really needed to avoid a NaN result. The test
		   for NaN is mainly for optimization. The test for 0 is useful
		   to get the correct sign (assuming mpfr_get_d supports signed
		   zeros on the implementation). */
		//if (s == 0 || DOUBLE_ISNAN (s) || DOUBLE_ISINF (s))
		if (s == 0 || std::isnan(s) || std::isinf(s))
		{
			// r = (long double) s;
			ret_dd[0] = s;
			ret_dd[1] = 0;
		}
		else
		{
			mpfr_t y, z;
			
			mpfr_init2 (y, mpfr_get_prec (x));
			mpfr_init2 (z, 53); /* keep the precision small */
			mpfr_set_d (z, s, MPFR_RNDN);  /* exact */

			// y := x - z = x - s
			mpfr_sub (y, x, z, MPFR_RNDN); /* exact */

			/* Add the second part of y (in the correct rounding mode). */
			//r = (long double) s + (long double) mpfr_get_d (y, rnd_mode);
			ret_dd[0] = s;
			ret_dd[1] = mpfr_get_d(y, rnd_mode);

			mpfr_clear (z);
			mpfr_clear (y);
		}
		//}
		//MPFR_SAVE_EXPO_FREE (expo);
		//return r;
	}
}

/* double-double code */
int mpfr_set_dd(mpfr_ptr r, double d[2], mpfr_rnd_t rnd_mode)
{
	mpfr_t t, u;
	int inexact, shift_exp;
	double h, l;
	//  MPFR_SAVE_EXPO_DECL (expo);
	
	/* Check for NAN */
	//ONGDOUBLE_NAN_ACTION (d, goto nan);
	if(std::isnan(d[0])) goto nan;

	/* Check for INF */
	//  if (d > MPFR_LDBL_MAX)
	if (d[0] > DBL_MAX)
	{
		mpfr_set_inf (r, 1);
		return 0;
	}
	else if (d[0] < -DBL_MAX)   //else if (d < -MPFR_LDBL_MAX)
	{
		mpfr_set_inf (r, -1);
	    return 0;
	}
	/* Check for ZERO */
	else if (d[0] == 0.0)
		return mpfr_set_d (r, (double) d[0], rnd_mode);

	// if (d >= (long double) MPFR_LDBL_MAX || d <= (long double) -MPFR_LDBL_MAX)
	//   h = (d >= (long double) MPFR_LDBL_MAX) ? MPFR_LDBL_MAX : -MPFR_LDBL_MAX;
	if (d[0] >= DBL_MAX || d[0] <= -DBL_MAX)
		h = (d[0] >= DBL_MAX) ? DBL_MAX : -DBL_MAX;
	else
		h = (double) d[0]; /* should not overflow */

	l = d[1]; //l = (double) (d - (long double) h);

	//MPFR_SAVE_EXPO_MARK (expo);

	mpfr_init2 (t, 53);
	mpfr_init2 (u, 53);

	inexact = mpfr_set_d (t, h, MPFR_RNDN);
	//MPFR_ASSERTN(inexact == 0);
	inexact = mpfr_set_d (u, l, MPFR_RNDN);
	//MPFR_ASSERTN(inexact == 0);
	inexact = mpfr_add (r, t, u, rnd_mode);

	mpfr_clear (t);
	mpfr_clear (u);

	//  MPFR_SAVE_EXPO_FREE (expo);
	inexact = mpfr_check_range (r, inexact, rnd_mode);
	return inexact;
	
nan:
	mpfr_set_nan(r); //MPFR_SET_NAN(r);
	//return MPFR_RET_NAN;
	return 0;

}

/* generic code */
// ret_qd[4] = ret_qd[high == 0], ret_qd[1], ret_qd[2], ret_qd[3]
void mpfr_get_qd(double ret_qd[4], mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
	int i;
	double s;

	//if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (x)))
	if (mpfr_is_singular(x))
	{
		//return (long double) mpfr_get_d (x, rnd_mode);
		ret_qd[0] = mpfr_get_d(x, rnd_mode);
		ret_qd[1] = (double)0.0;
		ret_qd[2] = (double)0.0;
		ret_qd[3] = (double)0.0;
	}
	else /* now x is a normal non-zero number */
	{
		//#if defined(HAVE_LDOUBLE_MAYBE_DOUBLE_DOUBLE)
		//if (MPFR_LDBL_MANT_DIG == 106)
		//{
		/* Assume double-double format (as found with the PowerPC ABI).
		   The generic code below isn't used because numbers with
		   precision > 106 would not be supported. */

		s = mpfr_get_d (x, MPFR_RNDN); /* high part of x */

		/* Let's first consider special cases separately. The test for
		   infinity is really needed to avoid a NaN result. The test
		   for NaN is mainly for optimization. The test for 0 is useful
		   to get the correct sign (assuming mpfr_get_d supports signed
		   zeros on the implementation). */
		//if (s == 0 || DOUBLE_ISNAN (s) || DOUBLE_ISINF (s))

		if (s == 0 || std::isnan(s) || std::isinf(s))
		{
			// r = (long double) s;
			ret_qd[0] = s;
			ret_qd[1] = 0;
			ret_qd[1] = 0;
			ret_qd[1] = 0;
		}
		else
		{
			mpfr_t y, z[3];
			
			mpfr_init2 (y, mpfr_get_prec (x));
			mpfr_init2 (z[0], 53); /* keep the precision small */
			mpfr_init2 (z[1], 53); /* keep the precision small */
			mpfr_init2 (z[2], 53); /* keep the precision small */
			
			// y[0] := x - z = x - s
			mpfr_set_d (z[0], s, MPFR_RNDN);  /* exact */
			mpfr_sub (y, x, z[0], MPFR_RNDN); /* exact */
			ret_qd[0] = s;
			
			// y[1] := y[0] - z = x - s
			s = mpfr_get_d (y, MPFR_RNDN); /* high part of y */
			mpfr_set_d (z[1], s, MPFR_RNDN);  /* exact */
			mpfr_sub (y, x, z[0], MPFR_RNDN); /* exact */
			mpfr_sub (y, y, z[1], MPFR_RNDN); /* exact */
			ret_qd[1] = s;
			
			// y[2] := y[1] - z = x - s
			s = mpfr_get_d (y, MPFR_RNDN); /* high part of y */
			mpfr_set_d (z[2], s, MPFR_RNDN);  /* exact */
			mpfr_sub (y, x, z[0], MPFR_RNDN); /* exact */
			mpfr_sub (y, y, z[1], MPFR_RNDN); /* exact */
			mpfr_sub (y, y, z[2], MPFR_RNDN); /* exact */
			ret_qd[2] = s;

			// y[3] := x - z = x - s
			ret_qd[3] = mpfr_get_d(y, rnd_mode);

			mpfr_clear (z[0]);
			mpfr_clear (z[1]);
			mpfr_clear (z[2]);
			mpfr_clear (y);
		}
		//}
		//MPFR_SAVE_EXPO_FREE (expo);
		//return r;
	}
}

/* double-double code */
int mpfr_set_qd(mpfr_ptr r, double d[4], mpfr_rnd_t rnd_mode)
{
	mpfr_t t, u, v, w;
	int inexact, shift_exp;
	double h, l[3];
	//  MPFR_SAVE_EXPO_DECL (expo);
	
	/* Check for NAN */
	//ONGDOUBLE_NAN_ACTION (d, goto nan);
	if(std::isnan(d[0])) goto nan;

	/* Check for INF */
	//  if (d > MPFR_LDBL_MAX)
	if (d[0] > DBL_MAX)
	{
		mpfr_set_inf (r, 1);
		return 0;
	}
	else if (d[0] < -DBL_MAX)	//else if (d < -MPFR_LDBL_MAX)
	{
		mpfr_set_inf (r, -1);
	    return 0;
	}
	/* Check for ZERO */
	else if (d[0] == 0.0)
		return mpfr_set_d (r, (double) d[0], rnd_mode);

	// if (d >= (long double) MPFR_LDBL_MAX || d <= (long double) -MPFR_LDBL_MAX)
	//   h = (d >= (long double) MPFR_LDBL_MAX) ? MPFR_LDBL_MAX : -MPFR_LDBL_MAX;
	if (d[0] >= DBL_MAX || d[0] <= -DBL_MAX)
		h = (d[0] >= DBL_MAX) ? DBL_MAX : -DBL_MAX;
	else
		h = (double) d[0]; /* should not overflow */

	//l = (double) (d - (long double) h);
	l[0] = d[1];
	l[1] = d[2];
	l[2] = d[3];
	
	//MPFR_SAVE_EXPO_MARK (expo);
	
	mpfr_init2 (t, 53);
	mpfr_init2 (u, 53);
	mpfr_init2 (v, 53);
	mpfr_init2 (w, 53);
	
	inexact = mpfr_set_d (t, h, MPFR_RNDN);
	inexact = mpfr_set_d (u, l[0], MPFR_RNDN);
	inexact = mpfr_set_d (v, l[1], MPFR_RNDN);
	inexact = mpfr_set_d (w, l[2], MPFR_RNDN);
	
	inexact = mpfr_add (r, t, u, rnd_mode);
	inexact = mpfr_add (r, r, v, rnd_mode);
	inexact = mpfr_add (r, r, w, rnd_mode);
	
	mpfr_clear (t);
	mpfr_clear (u);
	mpfr_clear (v);
	mpfr_clear (w);
	
	//  MPFR_SAVE_EXPO_FREE (expo);
	inexact = mpfr_check_range (r, inexact, rnd_mode);
	return inexact;

nan:
	mpfr_set_nan(r); //MPFR_SET_NAN(r);
	//return MPFR_RET_NAN;

	return 0;
}

} //end of extern "C" {

#endif // __MPFR_DD_QD_H__
