#include <RcppEigen.h>
#include "utils.h"
#include "tree_explainer.h"

// [[Rcpp::depends(RcppEigen)]]


// // [[Rcpp::export]]
// SEXP create_tree_explainer_cpp(const Rcpp::List& tree_model) {
//     Rcpp::XPtr<TreeExplainer> ptr(new TreeExplainer(tree_model), true);
//     return ptr;
// }


// // [[Rcpp::export]]
// SEXP create_xgboost_explainer(int max_depth, long double base_score, const Rcpp::List& xgb_trees) {
//     TreeExplainer* explainer = new TreeExplainer();

//     int tree_count = xgb_trees.size();
//     std::vector<SimpleTree> xgb_trees_vec;
//     xgb_trees_vec.reserve(tree_count);

//     for (int i = 0; i < tree_count; i++) {
//         Rcpp::List cur_tree = xgb_trees[i];

//         Eigen::VectorXi cl = Rcpp::as<Eigen::VectorXi>( cur_tree["children_left"] );
//         Eigen::VectorXi cr = Rcpp::as<Eigen::VectorXi>( cur_tree["children_right"] );
//         Eigen::VectorXi feat = Rcpp::as<Eigen::VectorXi>( cur_tree["feature"] );
//         Eigen::VectorXd thr = Rcpp::as<Eigen::VectorXd>( cur_tree["threshold"] );
//         int md = Rcpp::as<int>( cur_tree["max_depth"] );
//         Eigen::VectorXd ns = Rcpp::as<Eigen::VectorXd>( cur_tree["n_node_samples"] );
//         Eigen::VectorXd val = Rcpp::as<Eigen::VectorXd>( cur_tree["value"] );
//         int nn = Rcpp::as<int>( cur_tree["node_count"] );

//         SimpleTree cur_tree_obj = {cl, cr, feat, thr, md, ns, val, nn};
//         xgb_trees_vec.push_back(cur_tree_obj);
//     }

//     explainer->initialize_xgboost(max_depth, base_score, xgb_trees_vec);

//     Rcpp::XPtr<TreeExplainer> ptr(explainer, true);
//     return ptr;
// }

// // [[Rcpp::export]]
// Rcpp::List get_tree_summary(SEXP explainer_ptr) {
//     Rcpp::XPtr<TreeExplainer> ptr(explainer_ptr);
//     TreeSummary summary = ptr->get_summary();
    
//     Rcpp::List result;
//     result["children_left"] = Rcpp::wrap(summary.children_left);
//     result["children_right"] = Rcpp::wrap(summary.children_right);
//     result["feature"] = Rcpp::wrap(summary.feature);
//     result["feature_uniq"] = Rcpp::wrap(summary.feature_uniq);
//     result["threshold"] = Rcpp::wrap(summary.threshold);
//     result["max_depth"] = summary.max_depth;
//     result["sample_weight"] = Rcpp::wrap(summary.sample_weight);
//     result["init_prediction"] = Rcpp::wrap(summary.init_prediction);
//     result["node_count"] = summary.node_count;
    
//     return result;
// }

// // [[Rcpp::export]]
// Eigen::VectorXd rsq(
//     SEXP explainer_ptr, 
//     const Eigen::MatrixXd& x,
//     const Eigen::VectorXd& y,
//     const Eigen::MatrixXd& T0_x
// ) {
//     Rcpp::XPtr<TreeExplainer> ptr(explainer_ptr);
//     return ptr->get_rsq(x, y, T0_x);
// }
