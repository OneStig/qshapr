#include <RcppEigen.h>

#include "qshapr_utils.h"

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::VectorXd r_inv_binom_coef(int d) {
    return inv_binom_coef(d);
}
