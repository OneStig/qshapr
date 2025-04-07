#ifndef QSHAPR_TREE_PARSER_H
#define QSHAPR_TREE_PARSER_H

#include <RcppEigen.h>
#include <vector>
#include <string>
#include <map>
#include <optional>

#include "utils.h"
#include "qshap.h"

class TreeExplainer {
private:
    TreeSummary summary;
    Eigen::MatrixXcd store_v_invc;
    Eigen::MatrixXcd store_z;
    
public:
    TreeExplainer();
    TreeExplainer(const Rcpp::List& tree_model);
    TreeSummary get_summary() const { return summary; }
    // ~TreeExplainer();
    
    Eigen::MatrixXd get_loss(
        const Eigen::MatrixXd& x,
        const Eigen::VectorXd& y,
        const Eigen::MatrixXd& T0_x,
        double y_mean_ori = std::numeric_limits<double>::quiet_NaN()
    );
    
    Eigen::VectorXd get_rsq(
        const Eigen::MatrixXd& x,
        const Eigen::VectorXd& y,
        const Eigen::MatrixXd& T0_x,
        bool loss_out = false,
        int ncore = 1,
        int nsample = -1,
        double nfrac = -1.0,
        int random_state = 42
    );
};

#endif 
