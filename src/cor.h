#ifndef _ConAn_COR_H
#define _ConAn_COR_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include <math.h>
#include <Rcpp.h>
#include <stdlib.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// Functions to export to R
RcppExport SEXP S_pcor(SEXP R_x);
RcppExport SEXP S_erase_vals(SEXP R_x, SEXP R_ix);
RcppExport SEXP S_cv(SEXP R_x);
RcppExport SEXP S_bgcv(SEXP R_x, SEXP R_ix);
RcppExport SEXP S_mc(SEXP R_x);
RcppExport SEXP S_bgmc(SEXP R_x, SEXP R_ix);
RcppExport SEXP S_mdc(SEXP R_xr, SEXP R_xt, SEXP R_bgr, SEXP R_bgt);

// Functions to be used within C++
mat pcor(mat x);
mat erase_vals(mat x, NumericVector ix);

NumericVector extract_cv(mat x);
NumericVector cv(mat x);
NumericVector bgcv(mat x, NumericVector ix);

double mean_cv(mat x);
double mc(mat x);
double bgmc(mat x, NumericVector ix);

#endif
