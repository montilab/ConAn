#include <R.h>
#include <omp.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdlib.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

//' Pearson correlation
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat pcor(const arma::mat& x) {
    return cor(x);
}

//' Erase pairwaise module genes from matrix
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat erase_mods(arma::mat x, const NumericVector& ix) {
    int len = ix.size();
    for (int i = 0; i < len; i++) { 
        x[ix[i]] = R_NaN;
    }
    return x;
}

//' Z-transformation of lower triangle of matrix
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector atanh_lower_tri(const arma::mat& x) {
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

//' Connectivity vector
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector atanh_lower_tri_pcor(const arma::mat& x) {
    return atanh_lower_tri( pcor(x) );
}

//' Background connectivity vector
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector atanh_lower_tri_erase_mods_pcor(const arma::mat& x, const NumericVector& ix) {
    return atanh_lower_tri( erase_mods( pcor(x), ix) );
}

//' Mean of background connevtivity vector
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double mean_atanh_lower_tri_erase_mods_pcor(const arma::mat& x, const NumericVector& ix) {
    return mean(atanh_lower_tri_erase_mods_pcor(x, ix));
}

//' Background-corrected modular connectivity
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector bg_corrected_atanh_lower_tri_pcor(const arma::mat& x, const double& bg) {
    return atanh_lower_tri_pcor(x) - bg;
}

//' Mean of Background-corrected modular connectivity
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double mean_bg_corrected_atanh_lower_tri_pcor(const arma::mat& x, const double& bg) {
    return mean( bg_corrected_atanh_lower_tri_pcor(x, bg) );
}

//' Modular differential connectivity
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double modular_differential_connectivity(const arma::mat& xr, const arma::mat& xt, const double& bgr, const double& bgt, const int &type) {
    double cxr = mean( pow ( tanh( bg_corrected_atanh_lower_tri_pcor(xr, bgr) ), 2) );
    double cxt = mean( pow ( tanh( bg_corrected_atanh_lower_tri_pcor(xt, bgt) ), 2) );
    if (type == 1) {
        return cxt/cxr;
    }
    if (type == 2) {
        return cxt-cxr;
    }
}
