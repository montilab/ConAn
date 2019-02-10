#include <R.h>
#include "cor.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// Erase Matrix Values
mat erase_vals(mat x, NumericVector ix) {
    int len = ix.size();
    for (int i = 0; i < len; i++) { 
         x[ix[i]] = R_NaN;
    }
    return x;
}

// Extract Connectivity Vector
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

// Mean of Connectivity Vector
double mean_cv(mat x) {
    int n = x.n_cols;
    int nval = 0;
    double usum = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
            if (!ISNAN(x[i*n+j])) {
                usum += abs(x[i*n+j]);
                nval++;
            }
        }
    }
    return usum/nval;
}

// Pearson Correlation
mat pcor(mat x) {
    return cor(x);
}

// Connectivity Vector
NumericVector cv(mat x) {
    return extract_cv(cor(x));
}

// Background Connectivity Vector
NumericVector bgcv(mat x, NumericVector ix) {
    return extract_cv(erase_vals(cor(x), ix));
}

// Mean Connectivity
double mc(mat x) {
    return mean_cv(cor(x));
}

// Background mean connectivity
double bgmc(mat x, NumericVector ix) {
    return mean_cv(erase_vals(cor(x), ix));
}

// Module connectivity minus background
double mc_mbg(mat x, double bg) {
    return mean( mean ( square( tanh( cor(x) - bg ) ) ) );
}

// Module differential connectivity
double mdc(mat xr, mat xt, double bgr, double bgt) {
    double cxr = mc_mbg(xr, bgr);
    double cxt = mc_mbg(xt, bgt);

    double dc = cxr/cxt;

    return dc;
}

// Module Differential Connectivity
SEXP S_mdc(SEXP R_xr, SEXP R_xt, SEXP R_bgr, SEXP R_bgt) {
    NumericMatrix xr(R_xr);
    arma::mat mr = Rcpp::as<arma::mat>(xr);

    NumericMatrix xt(R_xt);
    arma::mat mt = Rcpp::as<arma::mat>(xt);

    double bgr = as<double>(R_bgr);
    double bgt = as<double>(R_bgt);

    return wrap(mdc(mr, mt, bgr, bgt));
}

//
// Rccp Wrappers
//

// Pearson Correlation Wrapper
SEXP S_pcor(SEXP R_x) {
    NumericMatrix x(R_x);
    arma::mat m = Rcpp::as<arma::mat>(x);
    return wrap(pcor(m));
}

// Background Mean Connectivity Wrapper
SEXP S_erase_vals(SEXP R_x, SEXP R_ix) {
    NumericMatrix x(R_x);
    arma::mat m = Rcpp::as<arma::mat>(x);
    NumericVector ix(R_ix);
    return wrap(erase_vals(m, ix));
}

// Mean Connectivity Wrapper
SEXP S_cv(SEXP R_x) {
    NumericMatrix x(R_x);
    arma::mat m = Rcpp::as<arma::mat>(x);
    return wrap(cv(m));
}

// Background Mean Connectivity Wrapper
SEXP S_bgcv(SEXP R_x, SEXP R_ix) {
    NumericMatrix x(R_x);
    arma::mat m = Rcpp::as<arma::mat>(x);
    NumericVector ix(R_ix);
    return wrap(bgcv(m, ix));
}

// Mean Connectivity Wrapper
SEXP S_mc(SEXP R_x) {
    NumericMatrix x(R_x);
    arma::mat m = Rcpp::as<arma::mat>(x);
    return wrap(mc(m));
}

// Background Mean Connectivity Wrapper
SEXP S_bgmc(SEXP R_x, SEXP R_ix) {
    NumericMatrix x(R_x);
    arma::mat m = Rcpp::as<arma::mat>(x);
    NumericVector ix(R_ix);
    return wrap(bgmc(m, ix));
}
