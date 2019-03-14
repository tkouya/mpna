// logistic 写像
#include <iostream>
#include <iomanip>

// QD
#include "qd_real.h"

// Multiple precision with MPFR/GMP
#include "mpreal.h"

using namespace std;
using namespace mpfr;


// based on get_ld.c code in MPFR 4.0.0 
#define MPFR_EXP_ZERO (mpfr_get_emin() + 1)
#define MPFR_EXP_NAN  (mpfr_get_emin() + 2)
#define MPFR_EXP_INF  (mpfr_get_emin() + 3)
#define MPFR_EXP_UBF  (mpfr_get_emin() + 4)

// FALSE(=0), TRUE(=1)
//#define MPFR_IS_SINGULAR(x) (MPFR_EXP(x) <= MPFR_EXP_INF)
int mpfr_is_singular(mpfr_srcptr x)
{
	return (mpfr_get_exp(x) <= MPFR_EXP_INF);
}

/* generic code */
// ret_dd[2] = ret_dd[high == 0], ret_dd[low == 1]
void mpfr_get_dd (double ret_dd[2], mpfr_srcptr x, mpfr_rnd_t rnd_mode)
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
          if (s == 0 || isnan(s) || isinf(s))
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
int mpfr_set_dd (mpfr_ptr r, double d[2], mpfr_rnd_t rnd_mode)
{
  mpfr_t t, u;
  int inexact, shift_exp;
  double h, l;
//  MPFR_SAVE_EXPO_DECL (expo);

  /* Check for NAN */
  //ONGDOUBLE_NAN_ACTION (d, goto nan);
	if(isnan(d[0])) goto nan;

  /* Check for INF */
//  if (d > MPFR_LDBL_MAX)
  if (d[0] > DBL_MAX)
    {
      mpfr_set_inf (r, 1);
      return 0;
    }
  //else if (d < -MPFR_LDBL_MAX)
  else if (d[0] < -DBL_MAX)
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
	l = d[1];

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
void mpfr_get_qd (double ret_qd[4], mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
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
          if (s == 0 || isnan(s) || isinf(s))
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
int mpfr_set_qd (mpfr_ptr r, double d[4], mpfr_rnd_t rnd_mode)
{
  mpfr_t t, u, v, w;
  int inexact, shift_exp;
  double h, l[3];
//  MPFR_SAVE_EXPO_DECL (expo);

  /* Check for NAN */
  //ONGDOUBLE_NAN_ACTION (d, goto nan);
	if(isnan(d[0])) goto nan;

  /* Check for INF */
//  if (d > MPFR_LDBL_MAX)
  if (d[0] > DBL_MAX)
    {
      mpfr_set_inf (r, 1);
      return 0;
    }
  //else if (d < -MPFR_LDBL_MAX)
  else if (d[0] < -DBL_MAX)
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
  inexact = mpfr_add (r, v, w, rnd_mode);

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


int main(int argc, char *argv[])
{
	int i;
	int num_bits, num_decimal;

	// 引数チェック
	if(argc <= 1)
	{
		cerr << "Usage: " << argv[0] << " [num_bits]" << endl;
		return 0;
	}

	// 計算桁数設定
	num_bits = atoi(argv[1]);
	if(num_bits <= 24)
		num_bits = 24;

	num_decimal = (int)ceil(log10(2.0) * (double)num_bits);
	mpreal::set_default_prec(num_bits);

	cout << "num_bits = " << num_bits << ", num_decimal = " << num_decimal << endl;

	dd_real dd_x[102], mpfr_ddreal; // dd
	qd_real qd_x[102], mpfr_qdreal; // qd
	mpreal x[102], dd_mpfr, qd_mpfr, reldiff_dd, reldiff_qd, reldiff_mpfr; // mpfr
	mpreal long_x[102]; // longer
	double mpfr_dd[2], mpfr_qd[4];

	// 2倍の桁数で初期化
	for(i = 0; i < 102; i++)
		long_x[i].set_prec(num_bits * 2);

	// fix FPU mode for QD
	fpu_fix_start(NULL);

	// 初期値
	x[0] = "0.7501";
	dd_x[0] = "0.7501";
	qd_x[0] = "0.7501";
	long_x[0] = "0.7501";

	for(i = 0; i <= 100; i++)
	{
		x[i + 1] = 4 * x[i] * (1 - x[i]);
		dd_x[i + 1] = 4 * dd_x[i] * (1 - dd_x[i]);
		qd_x[i + 1] = 4 * qd_x[i] * (1 - qd_x[i]);
		long_x[i + 1] = 4 * long_x[i] * (1 - long_x[i]);

		if((i % 10) == 0)
		{
			//cout << setw(5) << scientific << i << setprecision(num_decimal) << ", " << x[i] << setprecision(32) << ", " << dd_x[i] << endl;
			//cout << setw(5) << scientific << i << setprecision(num_decimal) << ", " << x[i] << setprecision(32) << ", " << dd_x[i] << ", " << qd_x[i] << endl;

			// dd := mpfr
			mpfr_get_dd(mpfr_dd, x[i].mpfr_srcptr(), MPFR_RNDN);
			mpfr_ddreal = dd_real(mpfr_dd[0], mpfr_dd[1]);
			//cout << setw(5) << scientific << i << setprecision(32) << ", " << mpfr_ddreal << setprecision(32) << ", " << dd_x[i] << endl;
			// mpfr := dd
			mpfr_set_dd(dd_mpfr.mpfr_ptr(), dd_x[i].x, MPFR_RNDN);

			// qd := mpfr
			mpfr_get_qd(mpfr_qd, x[i].mpfr_srcptr(), MPFR_RNDN);
			mpfr_qdreal = qd_real(mpfr_qd[0], mpfr_qd[1], mpfr_qd[2], mpfr_qd[3]);

			cout << setw(5) << scientific << i << setprecision(64) << ", " << mpfr_qdreal << endl;
			cout << setw(5) << scientific << i << setprecision(64) << ", " << qd_x[i] << endl;

			// mpfr := qd
			mpfr_set_qd(qd_mpfr.mpfr_ptr(), qd_x[i].x, MPFR_RNDN);

			// reldiff_xd := |x[i] - xd_mpfr| / x[i]
			mpfr_reldiff(reldiff_dd.mpfr_ptr(), long_x[i].mpfr_srcptr(), dd_mpfr.mpfr_srcptr(), MPFR_RNDN);
			mpfr_reldiff(reldiff_qd.mpfr_ptr(), long_x[i].mpfr_srcptr(), qd_mpfr.mpfr_srcptr(), MPFR_RNDN);
			mpfr_reldiff(reldiff_mpfr.mpfr_ptr(), long_x[i].mpfr_srcptr(), x[i].mpfr_srcptr(), MPFR_RNDN);

			//cout << setw(5) << scientific << i << setprecision(5) << ", " << reldiff_mpfr << ", " << reldiff_dd << ", " << reldiff_qd << endl;
		}
	}

	return 0;
}
