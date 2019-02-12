#include <R.h>
#include "connectivity.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// Pearson correlation
mat pcor(const mat& x) {
    return cor(x);
}

// Erase pairwaise module genes from matrix
mat erase_mods(mat x, const NumericVector& ix) {
    int len = ix.size();
    for (int i = 0; i < len; i++) { 
        x[ix[i]] = R_NaN;
    }
    return x;
}

// Z-transformation of lower triangle of matrix
NumericVector atanh_lower_tri(const mat& x) {
    int n = x.n_cols, vl = (n*n-n)/2;
    NumericVector v(vl);
    int ix = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
            if (!ISNAN( x[i*n+j] )) {
                v[ix] = atanh( x[i*n+j] ); ix++;
            }
        }
    }
    v.erase (v.begin()+ix, v.end()); // Delete empty end of vector
    return v;
}

// Connectivity vector
NumericVector atanh_lower_tri_pcor(const mat& x) {
    return atanh_lower_tri( pcor(x) );
}

// Background connectivity vector
NumericVector atanh_lower_tri_erase_mods_pcor(const mat& x, const NumericVector& ix) {
    return atanh_lower_tri( erase_mods( pcor(x), ix) );
}

// Mean of background connevtivity vector
double mean_atanh_lower_tri_erase_mods_pcor(const mat& x, const NumericVector& ix) {
    return mean(atanh_lower_tri_erase_mods_pcor(x, ix));
}

// Background-corrected modular connectivity
NumericVector bg_corrected_atanh_lower_tri_pcor(const mat& x, const double& bg) {
    return atanh_lower_tri_pcor(x) - bg;
}

// Mean of Background-corrected modular connectivity
double mean_bg_corrected_atanh_lower_tri_pcor(const mat& x, const double& bg) {
    return mean( bg_corrected_atanh_lower_tri_pcor(x, bg) );
}

// Modular differential connectivity
double modular_differential_connectivity(const mat& xr, const mat& xt, const double& bgr, const double& bgt, const int &type) {
    double cxr = mean( pow ( tanh( bg_corrected_atanh_lower_tri_pcor(xr, bgr) ), 2) );
    double cxt = mean( pow ( tanh( bg_corrected_atanh_lower_tri_pcor(xt, bgt) ), 2) );
    if (type == 1) {
        return cxt/cxr;
    }
    if (type == 2) {
        return cxt-cxr;
    }
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
SEXP S_erase_mods(SEXP R_x, SEXP R_ix) {
    NumericMatrix x(R_x);
    arma::mat m = Rcpp::as<arma::mat>(x);
    NumericVector ix(R_ix);
    return wrap(erase_mods(m, ix));
}

// Connectiviy vector
SEXP S_atanh_lower_tri_pcor(SEXP R_x) {
    NumericMatrix x(R_x);
    arma::mat m = Rcpp::as<arma::mat>(x);
    return wrap(atanh_lower_tri_pcor(m));
}

// Background connectivity vector
SEXP S_atanh_lower_tri_erase_mods_pcor(SEXP R_x, SEXP R_ix) {
    NumericMatrix x(R_x);
    arma::mat m = Rcpp::as<arma::mat>(x);
    NumericVector ix(R_ix);
    return wrap(atanh_lower_tri_erase_mods_pcor(m, ix));
}

// Background mean connectivity
SEXP S_mean_atanh_lower_tri_erase_mods_pcor(SEXP R_x, SEXP R_ix) {
    NumericMatrix x(R_x);
    arma::mat m = Rcpp::as<arma::mat>(x);
    NumericVector ix(R_ix);
    return wrap(mean_atanh_lower_tri_erase_mods_pcor(m, ix));
}

SEXP S_bg_corrected_atanh_lower_tri_pcor(SEXP R_x, SEXP R_bg) {
    NumericMatrix x(R_x);
    arma::mat m = Rcpp::as<arma::mat>(x);
    double bg = as<double>(R_bg);   
    return wrap(bg_corrected_atanh_lower_tri_pcor(m, bg));
}

SEXP S_mean_bg_corrected_atanh_lower_tri_pcor(SEXP R_x, SEXP R_bg) {
    NumericMatrix x(R_x);
    arma::mat m = Rcpp::as<arma::mat>(x);
    double bg = as<double>(R_bg);   
    return wrap(mean_bg_corrected_atanh_lower_tri_pcor(m, bg));
}

// Modular differential connectivity
SEXP S_modular_differential_connectivity(SEXP R_xr, SEXP R_xt, SEXP R_bgr, SEXP R_bgt, SEXP R_type) {
    NumericMatrix xr(R_xr);
    arma::mat mr = Rcpp::as<arma::mat>(xr);

    NumericMatrix xt(R_xt);
    arma::mat mt = Rcpp::as<arma::mat>(xt);

    double bgr = as<double>(R_bgr);
    double bgt = as<double>(R_bgt);

    int type = as<int>(R_type);

    return wrap(modular_differential_connectivity(mr, mt, bgr, bgt, type));
}
