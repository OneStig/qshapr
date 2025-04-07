#ifndef QSHAPR_UTILS_H
#define QSHAPR_UTILS_H

#include <Eigen/Dense>
#include <complex>
#include <cmath>

struct TreeSummary {
    Eigen::VectorXi children_left;
    Eigen::VectorXi children_right;
    Eigen::VectorXi feature;
    Eigen::VectorXi feature_uniq;
    Eigen::VectorXd threshold;
    int max_depth;
    Eigen::VectorXd sample_weight;
    Eigen::VectorXd init_prediction;
    int node_count;
};

struct SimpleTree {
    Eigen::VectorXi children_left;
    Eigen::VectorXi children_right;
    Eigen::VectorXi feature;
    Eigen::VectorXd threshold;
    int max_depth;
    Eigen::VectorXd n_node_samples;
    Eigen::VectorXd value;
    int node_count;
};

Eigen::VectorXd inv_binom_coef(int d);

Eigen::MatrixXcd complex_v_invc_degree(int d);

Eigen::MatrixXcd store_complex_v_invc(int d);

Eigen::MatrixXcd store_complex_root(int d);

double complex_dot_v2(const Eigen::VectorXcd& p, const Eigen::VectorXcd& v_invc, int d);

void traversal_weight(
    const Eigen::VectorXd& x,
    int v,
    Eigen::VectorXd& w,
    const Eigen::VectorXi& children_left,
    const Eigen::VectorXi& children_right,
    const Eigen::VectorXi& feature,
    const Eigen::VectorXd& threshold,
    const Eigen::VectorXd& sample_weight,
    const Eigen::VectorXi& leaf_ind,
    Eigen::MatrixXd& w_res,
    Eigen::MatrixXi& w_ind,
    int depth,
    Eigen::VectorXi& met_feature
);

std::pair<Eigen::MatrixXd, Eigen::MatrixXi> weight(
    const Eigen::VectorXd& x,
    const TreeSummary& summary_tree
);

#endif
