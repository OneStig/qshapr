#ifndef QSHAP
#define QSHAP

#include <RcppEigen.h>
#include <complex>
#include <vector>

#include "utils.h"

void T2_sample(
    int i, 
    const Eigen::MatrixXd& w_matrix, 
    const Eigen::MatrixXi& w_ind, 
    const Eigen::VectorXd& init_prediction, 
    const Eigen::MatrixXcd& store_v_invc, 
    const Eigen::MatrixXcd& store_z, 
    Eigen::MatrixXd& shap_value, 
    const Eigen::VectorXi& feature_uniq
);

Eigen::MatrixXd T2(
    const Eigen::MatrixXd& x, 
    const TreeSummary& summary_tree,
    const Eigen::MatrixXcd& store_v_invc, 
    const Eigen::MatrixXcd& store_z,
    bool parallel = true
);

Eigen::MatrixXd loss_treeshap(
        const Eigen::MatrixXd& x,
        const Eigen::VectorXd& y,
        const TreeSummary& summary_tree,
        const Eigen::MatrixXcd& store_v_invc,
        const Eigen::MatrixXcd& store_z,
        const Eigen::MatrixXd& T0_x,
        double learning_rate = 1.0
);
    
#endif