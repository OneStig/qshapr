#include <RcppEigen.h>
#include "utils.h"
#include "tree_explainer.h"

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::MatrixXcd r_store_complex_v_invc(int d) {
    return store_complex_v_invc(d);
}

// [[Rcpp::export]]
Eigen::MatrixXcd r_store_complex_root(int d) {
    return store_complex_root(d);
}

// [[Rcpp::export]]
double r_complex_dot_v2(const Eigen::VectorXcd& p, const Eigen::VectorXcd& v_invc, int d) {
    return complex_dot_v2(p, v_invc, d);
}

// [[Rcpp::export]]
SEXP create_tree_explainer_cpp(const Rcpp::List& tree_model) {
    Rcpp::XPtr<TreeExplainer> ptr(new TreeExplainer(tree_model), true);
    return ptr;
}

// [[Rcpp::export]]
Rcpp::List get_tree_summary(SEXP explainer_ptr) {
    Rcpp::XPtr<TreeExplainer> ptr(explainer_ptr);
    TreeSummary summary = ptr->get_summary();
    
    Rcpp::List result;
    result["children_left"] = Rcpp::wrap(summary.children_left);
    result["children_right"] = Rcpp::wrap(summary.children_right);
    result["feature"] = Rcpp::wrap(summary.feature);
    result["feature_uniq"] = Rcpp::wrap(summary.feature_uniq);
    result["threshold"] = Rcpp::wrap(summary.threshold);
    result["max_depth"] = summary.max_depth;
    result["sample_weight"] = Rcpp::wrap(summary.sample_weight);
    result["init_prediction"] = Rcpp::wrap(summary.init_prediction);
    result["node_count"] = summary.node_count;
    
    return result;
}

// [[Rcpp::export]]
Eigen::VectorXd rsq(
    SEXP explainer_ptr, 
    const Eigen::MatrixXd& x,
    const Eigen::VectorXd& y,
    const Eigen::MatrixXd& T0_x
) {
    Rcpp::XPtr<TreeExplainer> ptr(explainer_ptr);
    return ptr->get_rsq(x, y, T0_x);
}
