#include <RcppEigen.h>

#include "qshapr_utils.h"

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::VectorXd r_inv_binom_coef(int d) {
    return inv_binom_coef(d);
}

// [[Rcpp::export]]
Eigen::MatrixXcd r_complex_v_invc_degree(int d) {
    return complex_v_invc_degree(d);
}
