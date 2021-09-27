//Doug - this is code from R - to condense into one file, I have manually substituted includes
//also changed / commented at Lines 1082 and 1093

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998-2016  The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/* Private header file for use during compilation of Mathlib */
#ifndef MATHLIB_PRIVATE_H
#define MATHLIB_PRIVATE_H

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

/* Required by C99 but might be slow */
#ifdef HAVE_LONG_DOUBLE
#  define LDOUBLE long double
#else
#  define LDOUBLE double
#endif

/* To ensure atanpi, cospi,  sinpi, tanpi are defined */
# ifndef __STDC_WANT_IEC_60559_FUNCS_EXT__
#  define __STDC_WANT_IEC_60559_FUNCS_EXT__ 1
# endif

#include <math.h>
#include <float.h> /* DBL_MIN etc */

//#include <Rconfig.h>
////////////////////DOUG SUBSTITUTED IN START////////////////////////
#ifndef R_RCONFIG_H
#define R_RCONFIG_H

#ifndef R_CONFIG_H

#define HAVE_F77_UNDERSCORE 1
#define IEEE_754 1
/* #undef HAVE_VISIBILITY_ATTRIBUTE */
#define ENABLE_NLS 1

/* All of the following assume a correct C99 implementation. */
#define HAVE_NEARBYINT 1
#define HAVE_RINT 1
#define HAVE_WORKING_ISFINITE 1
#define HAVE_HYPOT 1
#define HAVE_EXPM1 1
#define HAVE_LOG1P 1
#define HAVE_WORKING_LOG1P 1

#endif /* not R_CONFIG_H */

#endif /* not R_RCONFIG_H */
////////////////////DOUG SUBSTITUTED IN END////////////////////////

//#include <Rmath.h>
////////////////////DOUG SUBSTITUTED IN START////////////////////////
#ifndef RMATH_H
#define RMATH_H

/* needed for cospi etc */
#ifndef __STDC_WANT_IEC_60559_FUNCS_EXT__
# define __STDC_WANT_IEC_60559_FUNCS_EXT__ 1
#endif
#if defined(__cplusplus) && !defined(DO_NOT_USE_CXX_HEADERS)
# include <cmath>
// See comment in R.h
# ifdef __SUNPRO_CC
using namespace std;
# endif
#else
# include <math.h>
#endif

#ifdef NO_C_HEADERS
# warning "use of NO_C_HEADERS is defunct and will be ignored"
#endif

/*-- Mathlib as part of R --  define this for standalone : */
/* #undef MATHLIB_STANDALONE */

#define R_VERSION_STRING "@PACKAGE_VERSION@"

#ifndef HAVE_EXPM1
@RMATH_HAVE_EXPM1@
#endif

#ifndef HAVE_HYPOT
@RMATH_HAVE_HYPOT@
#endif

#ifndef HAVE_LOG1P
@RMATH_HAVE_LOG1P@
#endif

#ifndef HAVE_WORKING_LOG1P
@RMATH_HAVE_WORKING_LOG1P@
#endif

#if defined(HAVE_LOG1P) && !defined(HAVE_WORKING_LOG1P)
/* remap to avoid problems with getting the right entry point */
double  Rlog1p(double);
#define log1p Rlog1p
#endif


	/* Undo SGI Madness */

#ifdef __sgi
# ifdef ftrunc
#  undef ftrunc
# endif
# ifdef qexp
#  undef qexp
# endif
# ifdef qgamma
#  undef qgamma
# endif
#endif


/* ----- The following constants and entry points are part of the R API ---- */

/* 30 Decimal-place constants */
/* Computed with bc -l (scale=32; proper round) */

/* SVID & X/Open Constants */
/* Names from Solaris math.h */

#ifndef M_E
#define M_E		2.718281828459045235360287471353	/* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E		1.442695040888963407359924681002	/* log2(e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E	0.434294481903251827651128918917	/* log10(e) */
#endif

#ifndef M_LN2
#define M_LN2		0.693147180559945309417232121458	/* ln(2) */
#endif

#ifndef M_LN10
#define M_LN10		2.302585092994045684017991454684	/* ln(10) */
#endif

#ifndef M_PI
#define M_PI		3.141592653589793238462643383280	/* pi */
#endif

#ifndef M_2PI
#define M_2PI		6.283185307179586476925286766559	/* 2*pi */
#endif

#ifndef M_PI_2
#define M_PI_2		1.570796326794896619231321691640	/* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4		0.785398163397448309615660845820	/* pi/4 */
#endif

#ifndef M_1_PI
#define M_1_PI		0.318309886183790671537767526745	/* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI		0.636619772367581343075535053490	/* 2/pi */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI	1.128379167095512573896158903122	/* 2/sqrt(pi) */
#endif

#ifndef M_SQRT2
#define M_SQRT2		1.414213562373095048801688724210	/* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2	0.707106781186547524400844362105	/* 1/sqrt(2) */
#endif

/* R-Specific Constants */

#ifndef M_SQRT_3
#define M_SQRT_3	1.732050807568877293527446341506	/* sqrt(3) */
#endif

#ifndef M_SQRT_32
#define M_SQRT_32	5.656854249492380195206754896838	/* sqrt(32) */
#endif

#ifndef M_LOG10_2
#define M_LOG10_2	0.301029995663981195213738894724	/* log10(2) */
#endif

#ifndef M_SQRT_PI
#define M_SQRT_PI	1.772453850905516027298167483341	/* sqrt(pi) */
#endif

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI	0.398942280401432677939946059934	/* 1/sqrt(2pi) */
#endif

#ifndef M_SQRT_2dPI
#define M_SQRT_2dPI	0.797884560802865355879892119869	/* sqrt(2/pi) */
#endif


#ifndef M_LN_2PI
#define M_LN_2PI	1.837877066409345483560659472811	/* log(2*pi) */
#endif

#ifndef M_LN_SQRT_PI
#define M_LN_SQRT_PI	0.572364942924700087071713675677	/* log(sqrt(pi))
								   == log(pi)/2 */
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI	0.918938533204672741780329736406	/* log(sqrt(2*pi))
								 == log(2*pi)/2 */
#endif

#ifndef M_LN_SQRT_PId2
#define M_LN_SQRT_PId2	0.225791352644727432363097614947	/* log(sqrt(pi/2))
								   == log(pi/2)/2 */
#endif


#ifdef MATHLIB_STANDALONE
# ifndef R_EXT_BOOLEAN_H_
/* "copy-paste" R_ext/Boolean.h if not already included: */
 #define R_EXT_BOOLEAN_H_
 #undef FALSE
 #undef TRUE
 typedef enum { FALSE = 0, TRUE } Rboolean;
# endif
#else
//# include <R_ext/Boolean.h>	//this code is copied within Rmath.h
#ifndef R_EXT_BOOLEAN_H_
#define R_EXT_BOOLEAN_H_

#undef FALSE
#undef TRUE

#ifdef  __cplusplus
extern "C" {
#endif
typedef enum { FALSE = 0, TRUE /*, MAYBE */ } Rboolean;

#ifdef  __cplusplus
}
#endif

#endif /* R_EXT_BOOLEAN_H_ */
#endif


#if !defined(MATHLIB_STANDALONE) && !defined(R_NO_REMAP_RMATH)
#define bessel_i	Rf_bessel_i
#define bessel_j	Rf_bessel_j
#define bessel_k	Rf_bessel_k
#define bessel_y	Rf_bessel_y
#define bessel_i_ex	Rf_bessel_i_ex
#define bessel_j_ex	Rf_bessel_j_ex
#define bessel_k_ex	Rf_bessel_k_ex
#define bessel_y_ex	Rf_bessel_y_ex
#define beta		Rf_beta
#define choose		Rf_choose
#define dbeta		Rf_dbeta
#define dbinom		Rf_dbinom
#define dbinom_raw	Rf_dbinom_raw
#define dcauchy		Rf_dcauchy
#define dchisq		Rf_dchisq
#define dexp		Rf_dexp
#define df		Rf_df
#define dgamma		Rf_dgamma
#define dgeom		Rf_dgeom
#define dhyper		Rf_dhyper
#define digamma		Rf_digamma
#define dlnorm		Rf_dlnorm
#define dlogis		Rf_dlogis
#define dnbeta		Rf_dnbeta
#define dnbinom		Rf_dnbinom
#define dnbinom_mu	Rf_dnbinom_mu
#define dnchisq		Rf_dnchisq
#define dnf		Rf_dnf
#define dnorm4		Rf_dnorm4
#define dnt		Rf_dnt
#define dpois_raw	Rf_dpois_raw
#define dpois		Rf_dpois
#define dpsifn		Rf_dpsifn
#define dsignrank	Rf_dsignrank
#define dt		Rf_dt
#define dtukey		Rf_dtukey
#define dunif		Rf_dunif
#define dweibull	Rf_dweibull
#define dwilcox		Rf_dwilcox
#define fmax2		Rf_fmax2
#define fmin2		Rf_fmin2
#define fprec		Rf_fprec
#define fround		Rf_fround
#define ftrunc		Rf_ftrunc
#define fsign		Rf_fsign
#define gammafn		Rf_gammafn
#define imax2		Rf_imax2
#define imin2		Rf_imin2
#define lbeta		Rf_lbeta
#define lchoose		Rf_lchoose
#define lgammafn	Rf_lgammafn
#define lgammafn_sign	Rf_lgammafn_sign
#define lgamma1p	Rf_lgamma1p
#define log1pexp       	Rf_log1pexp
#define log1pmx		Rf_log1pmx
#define logspace_add	Rf_logspace_add
#define logspace_sub	Rf_logspace_sub
#define logspace_sum	Rf_logspace_sum
#define pbeta		Rf_pbeta
#define pbeta_raw	Rf_pbeta_raw
#define pbinom		Rf_pbinom
#define pcauchy		Rf_pcauchy
#define pchisq		Rf_pchisq
#define pentagamma	Rf_pentagamma
#define pexp		Rf_pexp
#define pf		Rf_pf
#define pgamma		Rf_pgamma
#define pgeom		Rf_pgeom
#define phyper		Rf_phyper
#define plnorm		Rf_plnorm
#define plogis		Rf_plogis
#define pnbeta		Rf_pnbeta
#define pnbinom		Rf_pnbinom
#define pnbinom_mu     	Rf_pnbinom_mu
#define pnchisq		Rf_pnchisq
#define pnf		Rf_pnf
#define pnorm5		Rf_pnorm5
#define pnorm_both	Rf_pnorm_both
#define pnt		Rf_pnt
#define ppois		Rf_ppois
#define psignrank	Rf_psignrank
#define psigamma	Rf_psigamma
#define pt		Rf_pt
#define ptukey		Rf_ptukey
#define punif		Rf_punif
#define pythag		Rf_pythag
#define pweibull	Rf_pweibull
#define pwilcox		Rf_pwilcox
#define qbeta		Rf_qbeta
#define qbinom		Rf_qbinom
#define qcauchy		Rf_qcauchy
#define qchisq		Rf_qchisq
#define qchisq_appr	Rf_qchisq_appr
#define qexp		Rf_qexp
#define qf		Rf_qf
#define qgamma		Rf_qgamma
#define qgeom		Rf_qgeom
#define qhyper		Rf_qhyper
#define qlnorm		Rf_qlnorm
#define qlogis		Rf_qlogis
#define qnbeta		Rf_qnbeta
#define qnbinom		Rf_qnbinom
#define qnbinom_mu     	Rf_qnbinom_mu
#define qnchisq		Rf_qnchisq
#define qnf		Rf_qnf
#define qnorm5		Rf_qnorm5
#define qnt		Rf_qnt
#define qpois		Rf_qpois
#define qsignrank	Rf_qsignrank
#define qt		Rf_qt
#define qtukey		Rf_qtukey
#define qunif		Rf_qunif
#define qweibull	Rf_qweibull
#define qwilcox		Rf_qwilcox
#define rbeta		Rf_rbeta
#define rbinom		Rf_rbinom
#define rcauchy		Rf_rcauchy
#define rchisq		Rf_rchisq
#define rexp		Rf_rexp
#define rf		Rf_rf
#define rgamma		Rf_rgamma
#define rgeom		Rf_rgeom
#define rhyper		Rf_rhyper
#define rlnorm		Rf_rlnorm
#define rlogis		Rf_rlogis
#define rmultinom	Rf_rmultinom
#define rnbeta		Rf_rnbeta
#define rnbinom		Rf_rnbinom
#define rnbinom_mu     	Rf_rnbinom_mu
#define rnchisq		Rf_rnchisq
#define rnf		Rf_rnf
#define rnorm		Rf_rnorm
#define rnt		Rf_rnt
#define rpois		Rf_rpois
#define rsignrank	Rf_rsignrank
#define rt		Rf_rt
#define rtukey		Rf_rtukey
#define runif		Rf_runif
#define rweibull	Rf_rweibull
#define rwilcox		Rf_rwilcox
#define sign		Rf_sign
#define tetragamma	Rf_tetragamma
#define trigamma	Rf_trigamma
#endif

#define dnorm dnorm4
#define pnorm pnorm5
#define qnorm qnorm5

#ifdef  __cplusplus
extern "C" {
#endif
	/* R's versions with !R_FINITE checks */

double R_pow(double x, double y);
double R_pow_di(double, int);

	/* Random Number Generators */

double	norm_rand(void);
double	unif_rand(void);
double	exp_rand(void);
#ifdef MATHLIB_STANDALONE
void	set_seed(unsigned int, unsigned int);
void	get_seed(unsigned int *, unsigned int *);
#endif

	/* Normal Distribution */

double	dnorm(double, double, double, int);
double	pnorm(double, double, double, int, int);
double	qnorm(double, double, double, int, int);
double	rnorm(double, double);
void	pnorm_both(double, double *, double *, int, int);/* both tails */

	/* Uniform Distribution */

double	dunif(double, double, double, int);
double	punif(double, double, double, int, int);
double	qunif(double, double, double, int, int);
double	runif(double, double);

	/* Gamma Distribution */

double	dgamma(double, double, double, int);
double	pgamma(double, double, double, int, int);
double	qgamma(double, double, double, int, int);
double	rgamma(double, double);

double  log1pmx(double);
double  log1pexp(double); // <-- ../nmath/plogis.c
double  lgamma1p(double);
double  logspace_add(double, double);
double  logspace_sub(double, double);
double  logspace_sum(const double *, int);

	/* Beta Distribution */

double	dbeta(double, double, double, int);
double	pbeta(double, double, double, int, int);
double	qbeta(double, double, double, int, int);
double	rbeta(double, double);

	/* Lognormal Distribution */

double	dlnorm(double, double, double, int);
double	plnorm(double, double, double, int, int);
double	qlnorm(double, double, double, int, int);
double	rlnorm(double, double);

	/* Chi-squared Distribution */

double	dchisq(double, double, int);
double	pchisq(double, double, int, int);
double	qchisq(double, double, int, int);
double	rchisq(double);

	/* Non-central Chi-squared Distribution */

double	dnchisq(double, double, double, int);
double	pnchisq(double, double, double, int, int);
double	qnchisq(double, double, double, int, int);
double	rnchisq(double, double);

	/* F Distibution */

double	df(double, double, double, int);
double	pf(double, double, double, int, int);
double	qf(double, double, double, int, int);
double	rf(double, double);

	/* Student t Distibution */

double	dt(double, double, int);
double	pt(double, double, int, int);
double	qt(double, double, int, int);
double	rt(double);

	/* Binomial Distribution */

double  dbinom_raw(double x, double n, double p, double q, int give_log);
double	dbinom(double, double, double, int);
double	pbinom(double, double, double, int, int);
double	qbinom(double, double, double, int, int);
double	rbinom(double, double);

	/* Multnomial Distribution */

void	rmultinom(int, double*, int, int*);

	/* Cauchy Distribution */

double	dcauchy(double, double, double, int);
double	pcauchy(double, double, double, int, int);
double	qcauchy(double, double, double, int, int);
double	rcauchy(double, double);

	/* Exponential Distribution */

double	dexp(double, double, int);
double	pexp(double, double, int, int);
double	qexp(double, double, int, int);
double	rexp(double);

	/* Geometric Distribution */

double	dgeom(double, double, int);
double	pgeom(double, double, int, int);
double	qgeom(double, double, int, int);
double	rgeom(double);

	/* Hypergeometric Distibution */

double	dhyper(double, double, double, double, int);
double	phyper(double, double, double, double, int, int);
double	qhyper(double, double, double, double, int, int);
double	rhyper(double, double, double);

	/* Negative Binomial Distribution */

double	dnbinom(double, double, double, int);
double	pnbinom(double, double, double, int, int);
double	qnbinom(double, double, double, int, int);
double	rnbinom(double, double);

double	dnbinom_mu(double, double, double, int);
double	pnbinom_mu(double, double, double, int, int);
double	qnbinom_mu(double, double, double, int, int);
double	rnbinom_mu(double, double);

	/* Poisson Distribution */

double	dpois_raw (double, double, int);
double	dpois(double, double, int);
double	ppois(double, double, int, int);
double	qpois(double, double, int, int);
double	rpois(double);

	/* Weibull Distribution */

double	dweibull(double, double, double, int);
double	pweibull(double, double, double, int, int);
double	qweibull(double, double, double, int, int);
double	rweibull(double, double);

	/* Logistic Distribution */

double	dlogis(double, double, double, int);
double	plogis(double, double, double, int, int);
double	qlogis(double, double, double, int, int);
double	rlogis(double, double);

	/* Non-central Beta Distribution */

double	dnbeta(double, double, double, double, int);
double	pnbeta(double, double, double, double, int, int);
double	qnbeta(double, double, double, double, int, int);
double	rnbeta(double, double, double);

	/* Non-central F Distribution */

double  dnf(double, double, double, double, int);
double	pnf(double, double, double, double, int, int);
double	qnf(double, double, double, double, int, int);

	/* Non-central Student t Distribution */

double	dnt(double, double, double, int);
double	pnt(double, double, double, int, int);
double	qnt(double, double, double, int, int);

	/* Studentized Range Distribution */

double	ptukey(double, double, double, double, int, int);
double	qtukey(double, double, double, double, int, int);

	/* Wilcoxon Rank Sum Distribution */

double dwilcox(double, double, double, int);
double pwilcox(double, double, double, int, int);
double qwilcox(double, double, double, int, int);
double rwilcox(double, double);

	/* Wilcoxon Signed Rank Distribution */

double dsignrank(double, double, int);
double psignrank(double, double, int, int);
double qsignrank(double, double, int, int);
double rsignrank(double);

	/* Gamma and Related Functions */
double	gammafn(double);
double	lgammafn(double);
double	lgammafn_sign(double, int*);
void    dpsifn(double, int, int, int, double*, int*, int*);
double	psigamma(double, double);
double	digamma(double);
double	trigamma(double);
double	tetragamma(double);
double	pentagamma(double);

double	beta(double, double);
double	lbeta(double, double);

double	choose(double, double);
double	lchoose(double, double);

	/* Bessel Functions */

double	bessel_i(double, double, double);
double	bessel_j(double, double);
double	bessel_k(double, double, double);
double	bessel_y(double, double);
double	bessel_i_ex(double, double, double, double *);
double	bessel_j_ex(double, double, double *);
double	bessel_k_ex(double, double, double, double *);
double	bessel_y_ex(double, double, double *);


	/* General Support Functions */

#ifndef HAVE_HYPOT
double 	hypot(double, double);
#endif
double 	pythag(double, double);
#ifndef HAVE_EXPM1
double  expm1(double); /* = exp(x)-1 {care for small x} */
#endif
#ifndef HAVE_LOG1P
double  log1p(double); /* = log(1+x) {care for small x} */
#endif
int	imax2(int, int);
int	imin2(int, int);
double	fmax2(double, double);
double	fmin2(double, double);
double	sign(double);
double	fprec(double, double);
double	fround(double, double);
double	fsign(double, double);
double	ftrunc(double);

double  log1pmx(double); /* Accurate log(1+x) - x, {care for small x} */
double  lgamma1p(double);/* accurate log(gamma(x+1)), small x (0 < x < 0.5) */

/* More accurate cos(pi*x), sin(pi*x), tan(pi*x)

   These declarations might clash with system headers if someone had
   already included math.h with __STDC_WANT_IEC_60559_FUNCS_EXT__
   defined (and we try, above).
   We can add a check for that via the value of
   __STDC_IEC_60559_FUNCS__ (>= 201506L).
*/
#if !(defined(__STDC_IEC_60559_FUNCS__) && __STDC_IEC_60559_FUNCS__ >= 201506L)
double cospi(double);
double sinpi(double);
double tanpi(double);
#endif

/* Compute the log of a sum or difference from logs of terms, i.e.,
 *
 *     log (exp (logx) + exp (logy))
 * or  log (exp (logx) - exp (logy))
 *
 * without causing overflows or throwing away too much accuracy:
 */
double  logspace_add(double logx, double logy);
double  logspace_sub(double logx, double logy);


/* ----------------- Private part of the header file ------------------- */

	/* old-R Compatibility */

#ifdef OLD_RMATH_COMPAT
# define snorm	norm_rand
# define sunif	unif_rand
# define sexp	exp_rand
#endif

#if defined(MATHLIB_STANDALONE) && !defined(MATHLIB_PRIVATE_H)
/* second is defined by nmath.h */

/* If isnan is a macro, as C99 specifies, the C++
   math header will undefine it. This happens on macOS */
# ifdef __cplusplus
  int R_isnancpp(double); /* in mlutils.c */
#  define ISNAN(x)     R_isnancpp(x)
# else
#  define ISNAN(x)     (isnan(x)!=0)
# endif

# define R_FINITE(x)    R_finite(x)
int R_finite(double);

# ifdef _WIN32  /* not Win32 as no config information */
#  ifdef RMATH_DLL
#   define R_EXTERN extern __declspec(dllimport)
#  else
#   define R_EXTERN extern
#  endif
R_EXTERN double NA_REAL;
R_EXTERN double R_PosInf;
R_EXTERN double R_NegInf;
R_EXTERN int N01_kind;
#  undef R_EXTERN
#else
extern int N01_kind;
# endif

#endif /* MATHLIB_STANDALONE */

#ifdef  __cplusplus
}
#endif

#endif /* RMATH_H */
////////////////////DOUG SUBSTITUTED IN END////////////////////////

/* Used internally only */
double  Rf_d1mach(int);
double	Rf_gamma_cody(double);

//#include <R_ext/RS.h>
////////////////////DOUG SUBSTITUTED IN START////////////////////////
#ifndef R_RS_H
#define R_RS_H

#if defined(__cplusplus) && !defined(DO_NOT_USE_CXX_HEADERS)
# include <cstring>
# include <cstddef>
# define R_SIZE_T std::size_t
#else
# include <string.h>		/* for memcpy, memset */
# include <stddef.h> /* for size_t */
# define R_SIZE_T size_t
#endif

//#include <Rconfig.h>		/* for F77_APPEND_UNDERSCORE */	//this code is copied within RS.h
#ifndef R_RCONFIG_H
#define R_RCONFIG_H

#ifndef R_CONFIG_H

#define HAVE_F77_UNDERSCORE 1
#define IEEE_754 1
/* #undef HAVE_VISIBILITY_ATTRIBUTE */
#define ENABLE_NLS 1

/* All of the following assume a correct C99 implementation. */
#define HAVE_NEARBYINT 1
#define HAVE_RINT 1
#define HAVE_WORKING_ISFINITE 1
#define HAVE_HYPOT 1
#define HAVE_EXPM1 1
#define HAVE_LOG1P 1
#define HAVE_WORKING_LOG1P 1

#endif /* not R_CONFIG_H */

#endif /* not R_RCONFIG_H */

#ifdef  __cplusplus
extern "C" {
#endif

/* S Like Error Handling */

//#include <R_ext/Error.h>	/* for error and warning */	//this code is copied within RS.h
#ifndef R_ERROR_H_
#define R_ERROR_H_

#ifdef  __cplusplus
extern "C" {
#endif

#if defined(__GNUC__) && __GNUC__ >= 3
#define NORET __attribute__((noreturn))
#else
#define NORET
#endif

void NORET Rf_error(const char *, ...);
void NORET UNIMPLEMENTED(const char *);
void NORET WrongArgCount(const char *);

void	Rf_warning(const char *, ...);
void 	R_ShowMessage(const char *s);
    

#ifdef  __cplusplus
}
#endif

#ifndef R_NO_REMAP
#define error Rf_error
#define warning Rf_warning
#endif


#endif /* R_ERROR_H_ */

#ifndef STRICT_R_HEADERS

#define R_PROBLEM_BUFSIZE	4096
/* Parentheses added for FC4 with gcc4 and -D_FORTIFY_SOURCE=2 */
#define PROBLEM			{char R_problem_buf[R_PROBLEM_BUFSIZE];(sprintf)(R_problem_buf,
#define MESSAGE                 {char R_problem_buf[R_PROBLEM_BUFSIZE];(sprintf)(R_problem_buf,
#define ERROR			),error(R_problem_buf);}
#define RECOVER(x)		),error(R_problem_buf);}
#define WARNING(x)		),warning(R_problem_buf);}
#define LOCAL_EVALUATOR		/**/
#define NULL_ENTRY		/**/
#define WARN			WARNING(NULL)

#endif

/* S Like Memory Management */

extern void *R_chk_calloc(R_SIZE_T, R_SIZE_T);
extern void *R_chk_realloc(void *, R_SIZE_T);
extern void R_chk_free(void *);

#ifndef STRICT_R_HEADERS
/* S-PLUS 3.x but not 5.x NULLs the pointer in the following */
#define Calloc(n, t)   (t *) R_chk_calloc( (R_SIZE_T) (n), sizeof(t) )
#define Realloc(p,n,t) (t *) R_chk_realloc( (void *)(p), (R_SIZE_T)((n) * sizeof(t)) )
#define Free(p)        (R_chk_free( (void *)(p) ), (p) = NULL)
#endif
#define R_Calloc(n, t)   (t *) R_chk_calloc( (R_SIZE_T) (n), sizeof(t) )
#define R_Realloc(p,n,t) (t *) R_chk_realloc( (void *)(p), (R_SIZE_T)((n) * sizeof(t)) )
#define R_Free(p)      (R_chk_free( (void *)(p) ), (p) = NULL)

#define Memcpy(p,q,n)  memcpy( p, q, (R_SIZE_T)(n) * sizeof(*p) )

/* added for 3.0.0 */
#define Memzero(p,n)  memset(p, 0, (R_SIZE_T)(n) * sizeof(*p))

#define CallocCharBuf(n) (char *) R_chk_calloc((R_SIZE_T) ((n)+1), sizeof(char))

/* S Like Fortran Interface */
/* These may not be adequate everywhere. Convex had _ prepending common
   blocks, and some compilers may need to specify Fortran linkage */

#ifdef HAVE_F77_UNDERSCORE
# define F77_CALL(x)	x ## _
#else
# define F77_CALL(x)	x
#endif
#define F77_NAME(x)    F77_CALL(x)
#define F77_SUB(x)     F77_CALL(x)
#define F77_COM(x)     F77_CALL(x)
#define F77_COMDECL(x) F77_CALL(x)

#ifndef NO_CALL_R
void	call_R(char*, long, void**, char**, long*, char**, long, char**);
#endif

#ifdef  __cplusplus
}
#endif

#endif /* R_RS_H */
////////////////////DOUG SUBSTITUTED IN END////////////////////////

/* possibly needed for debugging */
//#include <R_ext/Print.h>
////////////////////DOUG SUBSTITUTED IN START////////////////////////
#ifndef R_EXT_PRINT_H_
#define R_EXT_PRINT_H_

#ifdef  __cplusplus
/* If the vprintf interface is defined at all in C++ it may only be
   defined in namespace std.  It is part of the C++11 standard. */
# ifdef R_USE_C99_IN_CXX
#  include <cstdarg>
#  define R_VA_LIST std::va_list
# endif
extern "C" {
#else
# include <stdarg.h>
# define R_VA_LIST va_list
#endif

void Rprintf(const char *, ...);
void REprintf(const char *, ...);
#if !defined(__cplusplus) || defined R_USE_C99_IN_CXX
void Rvprintf(const char *, R_VA_LIST);
void REvprintf(const char *, R_VA_LIST);
#endif

#ifdef  __cplusplus
}
#endif

#endif /* R_EXT_PRINT_H_ */
////////////////////DOUG SUBSTITUTED IN END////////////////////////

/* moved from dpq.h */
#ifdef HAVE_NEARYINT
# define R_forceint(x)   nearbyint()
#else
# define R_forceint(x)   round(x)
#endif
//R >= 3.1.0: # define R_nonint(x) 	  (fabs((x) - R_forceint(x)) > 1e-7)
# define R_nonint(x) 	  (fabs((x) - R_forceint(x)) > 1e-7*fmax2(1., fabs(x)))

#ifndef MATHLIB_STANDALONE

//#include <R_ext/Error.h>
////////////////////DOUG SUBSTITUTED IN START////////////////////////
#ifndef R_ERROR_H_
#define R_ERROR_H_

#ifdef  __cplusplus
extern "C" {
#endif

#if defined(__GNUC__) && __GNUC__ >= 3
#define NORET __attribute__((noreturn))
#else
#define NORET
#endif

void NORET Rf_error(const char *, ...);
void NORET UNIMPLEMENTED(const char *);
void NORET WrongArgCount(const char *);

void	Rf_warning(const char *, ...);
void 	R_ShowMessage(const char *s);
    

#ifdef  __cplusplus
}
#endif

#ifndef R_NO_REMAP
#define error Rf_error
#define warning Rf_warning
#endif


#endif /* R_ERROR_H_ */
////////////////////DOUG SUBSTITUTED IN END////////////////////////

# define MATHLIB_ERROR(fmt,x)		error(fmt,x);
# define MATHLIB_WARNING(fmt,x)		warning(fmt,x)
# define MATHLIB_WARNING2(fmt,x,x2)	warning(fmt,x,x2)
# define MATHLIB_WARNING3(fmt,x,x2,x3)	warning(fmt,x,x2,x3)
# define MATHLIB_WARNING4(fmt,x,x2,x3,x4) warning(fmt,x,x2,x3,x4)
# define MATHLIB_WARNING5(fmt,x,x2,x3,x4,x5) warning(fmt,x,x2,x3,x4,x5)

//#include <R_ext/Arith.h>
////////////////////DOUG SUBSTITUTED IN START////////////////////////
#ifndef R_ARITH_H_
#define R_ARITH_H_

/* 
   This used to define _BSD_SOURCE to make declarations of isfinite
   and isnan visible in glibc.  But that was deprecated in glibc 2.20,
   and --std=c99 suffices nowadays.
*/

//#include <R_ext/libextern.h>	//this code is copied within Arith.h
#undef LibExtern
#undef LibImport
#undef LibExport

/* Don't try to include CYGWIN here: decorating some symbols breaks
   the auto-export that it relies on, even if R_DLL_BUILD were set. */
#ifdef _WIN32 /* _WIN32 as does not depend on config.h */
#define LibImport __declspec(dllimport)
#define LibExport __declspec(dllexport)
#else
#define LibImport
#define LibExport
#endif

#ifdef __MAIN__
#define LibExtern LibExport
#define extern
#elif defined(R_DLL_BUILD)
#define LibExtern extern
#else
#define LibExtern extern LibImport
#endif

#ifdef  __cplusplus
extern "C" {
#else
/* needed for isnan and isfinite, neither of which are used under C++ */
# include <math.h>
#endif

/* implementation of these : ../../main/arithmetic.c */
LibExtern double R_NaN;		/* IEEE NaN */
LibExtern double R_PosInf;	/* IEEE Inf */
LibExtern double R_NegInf;	/* IEEE -Inf */
LibExtern double R_NaReal;	/* NA_REAL: IEEE */
LibExtern int	 R_NaInt;	/* NA_INTEGER:= INT_MIN currently */
#ifdef __MAIN__
#undef extern
#undef LibExtern
#endif

#define NA_LOGICAL	R_NaInt
#define NA_INTEGER	R_NaInt
/* #define NA_FACTOR	R_NaInt  unused */
#define NA_REAL		R_NaReal
/* NA_STRING is a SEXP, so defined in Rinternals.h */

int R_IsNA(double);		/* True for R's NA only */
int R_IsNaN(double);		/* True for special NaN, *not* for NA */
int R_finite(double);		/* True if none of NA, NaN, +/-Inf */
#define ISNA(x)	       R_IsNA(x)

/* ISNAN(): True for *both* NA and NaN.
   NOTE: some systems do not return 1 for TRUE.
   Also note that C++ math headers specifically undefine
   isnan if it is a macro (it is on macOS and in C99),
   hence the workaround.  This code also appears in Rmath.h
*/
#ifdef __cplusplus
  int R_isnancpp(double); /* in arithmetic.c */
#  define ISNAN(x)     R_isnancpp(x)
#else
#  define ISNAN(x)     (isnan(x)!=0)
#endif

/* The following is only defined inside R */
#ifdef HAVE_WORKING_ISFINITE
/* isfinite is defined in <math.h> according to C99 */
# define R_FINITE(x)    isfinite(x)
#else
# define R_FINITE(x)    R_finite(x)
#endif

#ifdef  __cplusplus
}
#endif

#endif /* R_ARITH_H_ */
////////////////////DOUG SUBSTITUTED IN END////////////////////////

//Doug has redefined these
//#define ML_POSINF	R_PosInf
//#define ML_NEGINF	R_NegInf
//#define ML_NAN		R_NaN
#define ML_POSINF	INFINITY
#define ML_NEGINF	-INFINITY
#define ML_NAN		NAN

void R_CheckUserInterrupt(void);
/* Ei-ji Nakama reported that AIX 5.2 has calloc as a macro and objected
   to redefining it.  Tests added for 2.2.1 */
//Doug has edited these out
//#ifdef calloc
//# undef calloc
//#endif
//#define calloc R_chk_calloc
//#ifdef free
//# undef free
//#endif
//#define free R_chk_free

#ifdef ENABLE_NLS
#include <libintl.h>
////////////////////DOUG SUBSTITUTED IN START////////////////////////
#undef LibExtern
#undef LibImport
#undef LibExport

/* Don't try to include CYGWIN here: decorating some symbols breaks
   the auto-export that it relies on, even if R_DLL_BUILD were set. */
#ifdef _WIN32 /* _WIN32 as does not depend on config.h */
#define LibImport __declspec(dllimport)
#define LibExport __declspec(dllexport)
#else
#define LibImport
#define LibExport
#endif

#ifdef __MAIN__
#define LibExtern LibExport
#define extern
#elif defined(R_DLL_BUILD)
#define LibExtern extern
#else
#define LibExtern extern LibImport
#endif
////////////////////DOUG SUBSTITUTED IN END////////////////////////

#define _(String) gettext (String)
#else
#define _(String) (String)
#endif

#else
/* Mathlib standalone */

#include <stdio.h>
#include <stdlib.h> /* for exit */
#define MATHLIB_ERROR(fmt,x)	{ printf(fmt,x); exit(1); }
#define MATHLIB_WARNING(fmt,x)		printf(fmt,x)
#define MATHLIB_WARNING2(fmt,x,x2)	printf(fmt,x,x2)
#define MATHLIB_WARNING3(fmt,x,x2,x3)	printf(fmt,x,x2,x3)
#define MATHLIB_WARNING4(fmt,x,x2,x3,x4) printf(fmt,x,x2,x3,x4)
#define MATHLIB_WARNING5(fmt,x,x2,x3,x4,x5) printf(fmt,x,x2,x3,x4,x5)

#define ISNAN(x) (isnan(x)!=0)
// Arith.h defines it
#ifndef R_FINITE
#ifdef HAVE_WORKING_ISFINITE
/* isfinite is defined in <math.h> according to C99 */
# define R_FINITE(x)    isfinite(x)
#else
# define R_FINITE(x)    R_finite(x)
#endif
#endif
int R_finite(double);

#define ML_POSINF	(1.0 / 0.0)
#define ML_NEGINF	((-1.0) / 0.0)
#define ML_NAN		(0.0 / 0.0)

#define _(String) String
#endif /* standalone */

#define ML_VALID(x)	(!ISNAN(x))

#define ME_NONE		0
/*	no error */
#define ME_DOMAIN	1
/*	argument out of domain */
#define ME_RANGE	2
/*	value out of range */
#define ME_NOCONV	4
/*	process did not converge */
#define ME_PRECISION	8
/*	does not have "full" precision */
#define ME_UNDERFLOW	16
/*	and underflow occured (important for IEEE)*/

#define ML_ERR_return_NAN { ML_ERROR(ME_DOMAIN, ""); return ML_NAN; }

/* For a long time prior to R 2.3.0 ML_ERROR did nothing.
   We don't report ME_DOMAIN errors as the callers collect ML_NANs into
   a single warning.
 */
#define ML_ERROR(x, s) { \
   if(x > ME_DOMAIN) { \
       char *msg = ""; \
       switch(x) { \
       case ME_DOMAIN: \
	   msg = _("argument out of domain in '%s'\n");	\
	   break; \
       case ME_RANGE: \
	   msg = _("value out of range in '%s'\n");	\
	   break; \
       case ME_NOCONV: \
	   msg = _("convergence failed in '%s'\n");	\
	   break; \
       case ME_PRECISION: \
	   msg = _("full precision may not have been achieved in '%s'\n"); \
	   break; \
       case ME_UNDERFLOW: \
	   msg = _("underflow occurred in '%s'\n");	\
	   break; \
       } \
       MATHLIB_WARNING(msg, s); \
   } \
}

/* Wilcoxon Rank Sum Distribution */

#define WILCOX_MAX 50

#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
# define attribute_hidden
#endif

/* Formerly private part of Mathlib.h */

/* always remap internal functions */
#define bd0       	Rf_bd0
#define chebyshev_eval	Rf_chebyshev_eval
#define chebyshev_init	Rf_chebyshev_init
#define gammalims	Rf_gammalims
#define lfastchoose	Rf_lfastchoose
#define lgammacor	Rf_lgammacor
#define stirlerr       	Rf_stirlerr
#define pnchisq_raw   	Rf_pnchisq_raw
#define pgamma_raw   	Rf_pgamma_raw
#define pnbeta_raw   	Rf_pnbeta_raw
#define pnbeta2       	Rf_pnbeta2
#define bratio       	Rf_bratio

	/* Chebyshev Series */

int	attribute_hidden chebyshev_init(double*, int, double);
double	attribute_hidden chebyshev_eval(double, const double *, const int);

	/* Gamma and Related Functions */

void	attribute_hidden gammalims(double*, double*);
double	attribute_hidden lgammacor(double); /* log(gamma) correction */
double  attribute_hidden stirlerr(double);  /* Stirling expansion "error" */

double	attribute_hidden lfastchoose(double, double);

double  attribute_hidden bd0(double, double);

double  attribute_hidden pnchisq_raw(double, double, double, double, double,
				     int, Rboolean, Rboolean);
double  attribute_hidden pgamma_raw(double, double, int, int);
double	attribute_hidden pbeta_raw(double, double, double, int, int);
double  attribute_hidden qchisq_appr(double, double, double, int, int, double tol);
LDOUBLE attribute_hidden pnbeta_raw(double, double, double, double, double);
double	attribute_hidden pnbeta2(double, double, double, double, double, int, int);

int	Rf_i1mach(int);

/* From toms708.c */
void attribute_hidden bratio(double a, double b, double x, double y,
	    		     double *w, double *w1, int *ierr, int log_p);


#endif /* MATHLIB_PRIVATE_H */
