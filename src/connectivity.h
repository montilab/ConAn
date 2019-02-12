#ifndef _ConAn_CONNECTIVITY_H
#define _ConAn_CONNECTIVITY_H

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
RcppExport SEXP S_erase_mods(SEXP R_x, SEXP R_ix);
RcppExport SEXP S_atanh_lower_tri_pcor(SEXP R_x);
RcppExport SEXP S_atanh_lower_tri_erase_mods_pcor(SEXP R_x, SEXP R_ix);
RcppExport SEXP S_mean_atanh_lower_tri_erase_mods_pcor(SEXP R_x, SEXP R_ix);
RcppExport SEXP S_bg_corrected_atanh_lower_tri_pcor(SEXP R_x, SEXP R_bg);
RcppExport SEXP S_mean_bg_corrected_atanh_lower_tri_pcor(SEXP R_x, SEXP R_bg);
RcppExport SEXP S_modular_differential_connectivity(SEXP R_xr, SEXP R_xt, SEXP R_bgr, SEXP R_bgt, SEXP R_type);

// Functions to be used within C++
mat           pcor(const mat& x);
mat           erase_mods(mat x, const NumericVector& ix);
NumericVector atanh_lower_tri(const mat& x);
NumericVector atanh_lower_tri_pcor(const mat& x);
NumericVector atanh_lower_tri_erase_mods_pcor(const mat& x, const NumericVector& ix);
NumericVector bg_corrected_atanh_lower_tri_pcor(const mat& x, const double& bg);
double        mean_atanh_lower_tri_erase_mods_pcor(const mat& x, const NumericVector& ix);
double        mean_bg_corrected_atanh_lower_tri_pcor(const mat& x, const double& bg);
double        modular_differential_connectivity(const mat& xr, const mat& xt, const double& bgr, const double& bgt, const int &type);

#endif
