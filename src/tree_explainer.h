// #ifndef QSHAPR_TREE_PARSER_H
// #define QSHAPR_TREE_PARSER_H

// #include <RcppEigen.h>
// #include <vector>
// #include <string>
// #include <map>
// #include <optional>

// #include "utils.h"
// #include "qshap.h"

// class TreeExplainer {
// public:
//     enum class TreeType {
//         XGBOOST,
//         LIGHTGBM,
//         GBM
//     };

// private:
//     int max_depth;

//     TreeSummary summary;
//     Eigen::MatrixXcd store_v_invc;
//     Eigen::MatrixXcd store_z;
//     TreeType model_type;

//     // Specfiic to XGBoost
//     long double base_score;
//     std::vector<SimpleTree> xgb_trees;
    
// public:
//     TreeExplainer();
//     TreeExplainer(const Rcpp::List& tree_model);
//     TreeSummary get_summary() const { return summary; }
//     // ~TreeExplainer();
    
//     // XGBoost specific methods
//     void initialize_xgboost(const int max_depth, const long double base_score, const std::vector<SimpleTree> xgb_trees);
    
//     Eigen::MatrixXd get_loss(
//         const Eigen::MatrixXd& x,
//         const Eigen::VectorXd& y,
//         const Eigen::MatrixXd& T0_x,
//         double y_mean_ori = std::numeric_limits<double>::quiet_NaN()
//     );
    
//     Eigen::VectorXd get_rsq(
//         const Eigen::MatrixXd& x,
//         const Eigen::VectorXd& y,
//         const Eigen::MatrixXd& T0_x,
//         bool loss_out = false,
//         int ncore = 1,
//         int nsample = -1,
//         double nfrac = -1.0,
//         int random_state = 42
//     );
// };

// #endif 
