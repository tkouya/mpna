// 2-fold precision

#ifndef DDSIZE
	#define DDSIZE 2
#endif // DDSIZE

#ifndef QDSIZE
	#define QDSIZE 4
#endif // QDSIZE

#define _QD_SPLITTER 134217729.0               // = 2^27 + 1
#define _QD_SPLIT_THRESH 6.69692879491417e+299 // = 2^996

/*********** Basic Functions ************/
/* Computes fl(a+b) and err(a+b).  Assumes |a| >= |b|. */
double quick_two_sum(double a, double b, double *err)
{
	double s = a + b;
	*err = b - (s - a);
	return s;
}

/* Computes fl(a-b) and err(a-b).  Assumes |a| >= |b| */
double quick_two_diff(double a, double b, double *err)
{
	double s = a - b;
	*err = (a - s) - b;
	return s;
}

/* Computes fl(a+b) and err(a+b).  */
double two_sum(double a, double b, double *err)
{
	double s = a + b;
	double bb = s - a;
	*err = (a - (s - bb)) + (b - bb);
	return s;
}

/* Computes fl(a-b) and err(a-b).  */
double two_diff(double a, double b, double *err)
{
	double s = a - b;
	double bb = s - a;
	*err = (a - (s - bb)) - (b + bb);
	return s;
}

#ifndef QD_FMS
/* Computes high word and lo word of a */
void split(double a, double *hi, double *lo)
{
	double temp;

	if (a > _QD_SPLIT_THRESH || a < -_QD_SPLIT_THRESH)
	{
		a *= 3.7252902984619140625e-09;  // 2^-28
		temp = _QD_SPLITTER * a;
		*hi = temp - (temp - a);
		*lo = a - (*hi);
		*hi *= 268435456.0;          // 2^28
		*lo *= 268435456.0;          // 2^28
	}
	else
	{
		temp = _QD_SPLITTER * a;
		*hi = temp - (temp - a);
		*lo = a - (*hi);
	}
}
#endif // QD_FMS

// check FP_FAST_FMA
#ifdef FP_FAST_FMA
	#define QD_FMA(a, b, c) fma((a), (b), (c)) // a * b + c
	#define QD_FMS(a, b, c) QD_FMA((a), (b), -(c)) // a * b - c
#endif // FP_FAST_FMA


/* Computes fl(a*b) and err(a*b). */
double two_prod(double a, double b, double *err)
{
#ifdef QD_FMS
	double p = a * b;

	*err = QD_FMS(a, b, p);

	return p;
#else
	double a_hi, a_lo, b_hi, b_lo;
	double p = a * b;

	split(a, &a_hi, &a_lo);
	split(b, &b_hi, &b_lo);
	*err = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;

	return p;
#endif
}

#ifdef X86
#ifdef  _WIN32
#include <float.h>
#else

#ifdef HAVE_FPU_CONTROL_H
#include <fpu_control.h>
#endif

#ifndef _FPU_GETCW
#define _FPU_GETCW(x) asm volatile ("fnstcw %0":"=m" (x));
#endif

#ifndef _FPU_SETCW
#define _FPU_SETCW(x) asm volatile ("fldcw %0": :"m" (x));
#endif

#ifndef _FPU_EXTENDED
#define _FPU_EXTENDED 0x0300
#endif

#ifndef _FPU_DOUBLE
#define _FPU_DOUBLE 0x0200
#endif

#endif
#endif /* X86 */

void fpu_fix_start(unsigned int *old_cw) {
#ifdef X86
#ifdef _WIN32
#ifdef __BORLANDC__
  /* Win 32 Borland C */
  unsigned short cw = _control87(0, 0);
  _control87(0x0200, 0x0300);
  if (old_cw) {
    *old_cw = cw;
  }
#else
  /* Win 32 MSVC */
  unsigned int cw = _control87(0, 0);
  _control87(0x00010000, 0x00030000);
  if (old_cw) {
    *old_cw = cw;
  }
#endif
#else
  /* Linux */
  volatile unsigned short cw, new_cw;
  _FPU_GETCW(cw);

  new_cw = (cw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
  _FPU_SETCW(new_cw);
  
  if (old_cw) {
    *old_cw = cw;
  }
#endif
#endif
}

void fpu_fix_end(unsigned int *old_cw) {
#ifdef X86
#ifdef _WIN32

#ifdef __BORLANDC__
  /* Win 32 Borland C */
  if (old_cw) {
    unsigned short cw = (unsigned short) *old_cw;
    _control87(cw, 0xFFFF);
  }
#else
  /* Win 32 MSVC */
  if (old_cw) {
    _control87(*old_cw, 0xFFFFFFFF);
  }
#endif

#else
  /* Linux */
  if (old_cw) {
    int cw;
    cw = *old_cw;
    _FPU_SETCW(cw);
  }
#endif
#endif
}
/* Computes fl(a*a) and err(a*a).  Faster than the above method. */
double two_sqr(double a, double *err)
{
#ifdef QD_FMS
	double p = a * a;

	*err = QD_FMS(a, a, p);

	return p;
#else
	double hi, lo;
	double q = a * a;

	split(a, &hi, &lo);
	*err = ((hi * hi - q) + 2.0 * hi * lo) + lo * lo;

	return q;
#endif
}
