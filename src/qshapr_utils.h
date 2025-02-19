#ifndef QSHAPR_UTILS_H
#define QSHAPR_UTILS_H

#include <Eigen/Dense>
#include <complex>
#include <cmath>

Eigen::VectorXd inv_binom_coef(int d);

Eigen::MatrixXcd complex_v_invc_degree(int d);

#endif