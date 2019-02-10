#include <R.h>
#include "connectivity.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// Pearson correlation
mat pcor(mat x) {
    return cor(x);
}

// Erase matrix values
mat erase_vals(mat x, NumericVector ix) {
    int len = ix.size();
    for (int i = 0; i < len; i++) { 
         x[ix[i]] = R_NaN;
    }
    return x;
}

// Extract connectivity vector
NumericVector extract_cv(mat x) {
    int n = x.n_cols;
    NumericVector v((n*n-n)/2);
    int ix = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
            if (!ISNAN(x[i*n+j])) {
                v[ix] = atanh(x[i*n+j]); ix++;
            }
        }
    }
    return v;
}

// Connectivity vector
NumericVector cv(mat x) {
    return extract_cv(cor(x));
}

// Background connectivity vector
NumericVector bgcv(mat x, NumericVector ix) {
    return extract_cv(erase_vals(cor(x), ix));
}

// Mean connectivity
double mc(mat x) {
    return mean(cv(x));
}

// Background mean connectivity
double bgmc(mat x, NumericVector ix) {
    return mean(bgcv(x, ix));
}

// Modular connectivity minus background
double mc_mbg(mat x, double bg) {
    return mean( pow ( tanh( cv(x) - bg ), 2) );
}

// Modular differential connectivity
double mdc(mat xr, mat xt, double bgr, double bgt) {
    double cxr = mc_mbg(xr, bgr);
    double cxt = mc_mbg(xt, bgt);
    double dc = cxt/cxr;
    return dc;
}

//
// Rccp Wrappers
//

// Pearson correlation
SEXP S_pcor(SEXP R_x) {
    NumericMatrix x(R_x);
    arma::mat m = Rcpp::as<arma::mat>(x);
    return wrap(pcor(m));
}

// Erase matrix values
SEXP S_erase_vals(SEXP R_x, SEXP R_ix) {
    NumericMatrix x(R_x);
    arma::mat m = Rcpp::as<arma::mat>(x);
    NumericVector ix(R_ix);
    return wrap(erase_vals(m, ix));
}

// Connectiviy vector
SEXP S_cv(SEXP R_x) {
    NumericMatrix x(R_x);
    arma::mat m = Rcpp::as<arma::mat>(x);
    return wrap(cv(m));
}

// Background connectivity vector
SEXP S_bgcv(SEXP R_x, SEXP R_ix) {
    NumericMatrix x(R_x);
    arma::mat m = Rcpp::as<arma::mat>(x);
    NumericVector ix(R_ix);
    return wrap(bgcv(m, ix));
}

// Mean connectivity
SEXP S_mc(SEXP R_x) {
    NumericMatrix x(R_x);
    arma::mat m = Rcpp::as<arma::mat>(x);
    return wrap(mc(m));
}

// Background mean connectivity
SEXP S_bgmc(SEXP R_x, SEXP R_ix) {
    NumericMatrix x(R_x);
    arma::mat m = Rcpp::as<arma::mat>(x);
    NumericVector ix(R_ix);
    return wrap(bgmc(m, ix));
}

// Modular differential connectivity
SEXP S_mdc(SEXP R_xr, SEXP R_xt, SEXP R_bgr, SEXP R_bgt) {
    NumericMatrix xr(R_xr);
    arma::mat mr = Rcpp::as<arma::mat>(xr);

    NumericMatrix xt(R_xt);
    arma::mat mt = Rcpp::as<arma::mat>(xt);

    double bgr = as<double>(R_bgr);
    double bgt = as<double>(R_bgt);

    return wrap(mdc(mr, mt, bgr, bgt));
}