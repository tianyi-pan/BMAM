//
// src.c
// To compile for R use: R CMD SHLIB DeconvolveGH.c
//
// Contains all C fns used in Simulations for Bayesian MMM.
//
// The function is from package MMLB (https://github.com/mercaldo/MMLB)

// Include the R headers
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#include <R_ext/Parse.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#define ETA_TOLERANCE 1e-7
#define ETA_MAX_ITER 100
#define ETA_EPS 1e-7

/*--------------------------------------------------------------*/
static double expit ( double x)
{
  double out, expX;
  expX = exp(x);
  out = expX / (1.0+expX);
  return out;
}
/*--------------------------------------------------------------*/

// utility function to eval R functions
// NOTE: for information about how variable arguments work, look here:
//       http://www.gnu.org/software/libc/manual/html_node/Variadic-Functions.html
// parameters:
//    func  - the function name
//    env   - the environment you want the function to be called in
//            (most of the time this will probably be R_GlobalEnv)
//    nArgs - the number of arguments to the function
//    ...   - comma-seperated list of SEXP objects
SEXP evalR(const char *func, SEXP env, int nArgs, ...) {
  SEXP Rfunc, expr, retval, tail;
  va_list args;
  int i;

  // initialize variable-argument pointer
  va_start(args, nArgs);

  // this part finds the R object that the function lives in
  Rfunc = findFun(install(func), env);

  // set up an R function call
  PROTECT(expr = allocVector(LANGSXP, nArgs+1));
  SETCAR(expr, Rfunc);
  tail = CDR(expr);
  for(i = 1; i <= nArgs; i++) {
    SETCAR(tail, va_arg(args, SEXP));
    tail = CDR(tail);
  }

  // call the function
  retval = eval(expr, env);
  UNPROTECT(1);   // unprotects expr
  va_end(args);

  return retval;
}

// utility function to print an R object
void printR(SEXP obj) { evalR("print", R_GlobalEnv, 1, obj); }

// utility function to print C arrays as R vectors/matrices
void printMat(const double *array, int nrows, int ncols) {
  SEXP SEXP_tmp;
  int i;

  // assume that array is a regular vector if nrows == 1
  if (nrows == 1)
    PROTECT(SEXP_tmp = allocVector(REALSXP, ncols));
  else
    PROTECT(SEXP_tmp = allocMatrix(REALSXP, nrows, ncols));
  // populate the R object
  for(i = 0; i < nrows*ncols; i++) REAL(SEXP_tmp)[i] = array[i];

  // print it
  printR(SEXP_tmp);
  UNPROTECT(1);
}

//
// Using vectors eta, gamma and sigma, calculate delta using the convolution equation
// Here, gamma and sigma are vectors of length n and are Xgam %*% gam and Xsig %*% sig
//
SEXP DeconvolveGH_CALL1(SEXP SEXP_eta,
                       SEXP SEXP_gamma,
                       SEXP SEXP_sigma,
                       SEXP SEXP_z,
                       SEXP SEXP_w){
  //                           double offset ){
  double etai, gamma, sigma, z, w;
  int q_points = 0;
  int n, s, t;
  double deltai, dmuiddeltai, new_mui, delmu, mu;
  int j, count, eta_converge;
  double h_0, h_1, p_z;
  double *p_z_lag, p_z_lag_q; //, p_z_lag_start

  /*** Added for the change to .Call ***/
  double *S_eta, *S_gamma, *S_sigma, *S_z, *S_w;
  double *S_Delta_C; //*S_p_zlag, *S_flag,
  SEXP SEXP_Delta_C;

  PROTECT(SEXP_eta = coerceVector(SEXP_eta, REALSXP));
  S_eta = REAL(SEXP_eta);
  n     = LENGTH(SEXP_eta);

  PROTECT(SEXP_gamma = coerceVector(SEXP_gamma, REALSXP));
  S_gamma = REAL(SEXP_gamma);

  PROTECT(SEXP_sigma = coerceVector(SEXP_sigma, REALSXP));
  S_sigma = REAL(SEXP_sigma);

  PROTECT(SEXP_z = coerceVector(SEXP_z, REALSXP));
  S_z = REAL(SEXP_z);
  q_points = LENGTH(SEXP_z);

  PROTECT(SEXP_w = coerceVector(SEXP_w, REALSXP));
  S_w = REAL(SEXP_w);

  PROTECT(SEXP_Delta_C = allocVector(REALSXP,n));
  S_Delta_C = REAL(SEXP_Delta_C);

  p_z_lag = (double *)malloc((unsigned) (q_points)*sizeof(double));
  for(t=0; t<q_points; t++) *(p_z_lag+t) = 0.0;                     // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 0.5

  //printf("blah --------- %f ", offset);
  /*************************************/
  for (s=0; s<n; s++){

    gamma = *(S_gamma+s);
    sigma = *(S_sigma+s);
    etai  = *(S_eta+s);

    deltai = etai;
    mu     = expit(etai);

    eta_converge  = 0;
    count         = 0;

    do{

      dmuiddeltai = 0.0;
      new_mui = 0.0;

      for(j=0; j<q_points; j++){

        z          = *(S_z+j);
        w          = *(S_w+j);
        p_z_lag_q  = *(p_z_lag+j);

        h_0 = expit( deltai +         sigma*z );
        h_1 = expit( deltai + gamma + sigma*z );

        p_z = h_0 * (1.0-p_z_lag_q) + h_1 * p_z_lag_q;

        new_mui = new_mui + (w * p_z);

        dmuiddeltai = dmuiddeltai + (w * ( h_0 * (1-h_0) * (1.0-p_z_lag_q) +
          h_1 * (1-h_1) * (    p_z_lag_q) ));
      }

      delmu = (new_mui - mu) / (dmuiddeltai+ETA_EPS);

      if (delmu< -0.5) delmu = -0.5;
      if (delmu>  0.5) delmu =  0.5;

      deltai= deltai - delmu;

      if( fabs(delmu) < ETA_TOLERANCE ) eta_converge=1;
      count++;

    } while(count<ETA_MAX_ITER && !eta_converge );

    *(S_Delta_C+s) = deltai;

    for(j=0; j<q_points; j++){
      z            = *(S_z+j);
      p_z_lag_q    = *(p_z_lag+j);
      *(p_z_lag+j) = (expit( deltai +         sigma*z) *(1.0-p_z_lag_q))+
        (expit( deltai + gamma + sigma*z) *(    p_z_lag_q)) ;
    }

  }
  /*printf("%f", p_z_lag);*/
  free(p_z_lag);
  UNPROTECT(6);
  return SEXP_Delta_C;
}
