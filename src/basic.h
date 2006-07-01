
#ifndef GSL_VS_R_H
#define GSL_VS_R_H 1

#define NDEBUG 1
#include <assert.h>


#ifdef RF_GSL /* RF_GSL */


#ifndef MATHLIB_STANDALONE
#define MATHLIB_STANDALONE 1
#endif


//#include <R.h>
//#include <Rinternals.h>
//#include <Rmath.h>
//include <errno.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
//#include  "/usr/local/lib/R/include/Rmath.h" -- is in standalone.h
// problems: pnorm (extremes, MPPFcts), chol(direct), dsv(direct)

// #include <bits/nan.h> ???

#define RANDOMNUMBERGENERATOR gsl_rng_mt19937
//                            gsl_rng_taus,...
#define END_WITH_GSL ,gsl_rng *
#define END_WITH_RANDOM ,RANDOM
#define END_WITH_GSLRNG ,gsl_rng *RANDOM
#define EXTERN
#define PRINTF printf 
#define RF_NAN 1E308
#define RF_NEGINF -1E308
#define RF_INF 1E308
#define RF_ISNA ISNAN
#define PI 3.14159265358979323846264338328
#define GAUSS_RANDOM(SIGMA) gsl_rng_gaussian(RANDOM, SIGMA)
#define UNIFORM_RANDOM gsl_rng_uniform(RANDOM)
#define POISSON_RANDOM(x) gsl_ran_poisson(RANDOM, x)
#define SQRT2 1.414213562373095048801688724210
#define SQRTPI 1.772453850905516027298167483341
#define INVPI 0.318309886183790671537767526745
#define TWOPI 6.283185307179586476925286766559
#define PIHALF 1.570796326794896619231321691640
extern gsl_rng *RANDOM;


#else        /* RF_GSL */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#define END_WITH_GSL
#define END_WITH_GSLRNG
#define END_WITH_RANDOM 
#define EXTERN extern "C"
#define PRINTF Rprintf
#define RF_NAN NA_REAL 
#define RF_NEGINF R_NegInf
#define RF_INF R_PosInf
#define RF_ISNA ISNAN
#define GAUSS_RANDOM(SIGMA) rnorm(0.0, SIGMA)
#define UNIFORM_RANDOM unif_rand()
#define POISSON_RANDOM(x) rpois(x)
#define SQRT2 M_SQRT2
#define SQRTPI M_SQRT_PI
#define INVPI M_1_PI
//#define TWOPI M_2PI
#define TWOPI 6.283185307179586476925286766559
#define PIHALF M_PI_2 
#endif      /* RF_GSL */

#endif /* GSL_VS_R_H */


