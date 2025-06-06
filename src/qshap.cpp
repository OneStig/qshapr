#include "qshap.h"


Eigen::MatrixXd T2(
        const Eigen::MatrixXd& x, 
        const Rcpp::List& tree_summary,
        const Eigen::MatrixXcd& store_v_invc, 
        const Eigen::MatrixXcd& store_z,
        bool parallel
) {
    TreeSummary summary_tree = list_to_tree_summary(tree_summary);
    
    // Extract init_prediction from summary_tree and filter non-zero values
    std::vector<double> init_prediction_vec;
    for (int i = 0; i < summary_tree.init_prediction.size(); i++) {
        if (summary_tree.init_prediction(i) != 0) {
            init_prediction_vec.push_back(summary_tree.init_prediction(i));
        }
    }
    Eigen::Map<Eigen::VectorXd> init_prediction(init_prediction_vec.data(), init_prediction_vec.size());
    
    // Initialize shap_value
    Eigen::MatrixXd shap_value = Eigen::MatrixXd::Zero(x.rows(), x.cols());
    
    // Process each sample
    for (int i = 0; i < x.rows(); i++) {
        Eigen::VectorXd xi = x.row(i);
        
        // Call weight function from utils.h
        std::pair<Eigen::MatrixXd, Eigen::MatrixXi> w = weight(xi, summary_tree);
        Eigen::MatrixXd w_matrix = w.first;
        Eigen::MatrixXi w_ind = w.second;
        
        // Calculate T2 for this sample
        T2_sample(i, w_matrix, w_ind, init_prediction, store_v_invc, store_z, shap_value, summary_tree.feature_uniq);
    }
    
    return shap_value;
}

void T2_sample(
    int i, 
    const Eigen::MatrixXd& w_matrix, 
    const Eigen::MatrixXi& w_ind, 
    const Eigen::VectorXd& init_prediction, 
    const Eigen::MatrixXcd& store_v_invc, 
    const Eigen::MatrixXcd& store_z, 
    Eigen::MatrixXd& shap_value, 
    const Eigen::VectorXi& feature_uniq
) {
    // Calculate T2 for each sample
    int L = w_matrix.rows();
    
    for (int l1 = 0; l1 < L; l1++) {
        for (int l2 = l1; l2 < L; l2++) {
            double init_prediction_product = init_prediction(l1) * init_prediction(l2);
            
            std::vector<int> union_f12_vec;
            for (int j_idx = 0; j_idx < feature_uniq.size(); j_idx++) { // Renamed inner loop variable
                int feat = feature_uniq(j_idx);
                if ((w_ind(l1, feat) + w_ind(l2, feat)) >= 1) {
                    union_f12_vec.push_back(feat);
                }
            }
            
            if (union_f12_vec.empty()) continue; // Avoid issues if union is empty

            Eigen::Map<Eigen::VectorXi> union_f12(union_f12_vec.data(), union_f12_vec.size());
            
            int n12 = union_f12.size();
            if (n12 == 0) continue; // Should be caught by previous check, but good for safety

            int n12_c = n12 / 2 + 1;
            
            Eigen::VectorXcd v_invc = store_v_invc.row(n12).head(n12_c);
            Eigen::VectorXcd z_roots = store_z.row(n12).head(n12_c); // Renamed z to z_roots to avoid conflict
            
            // Pre-calculate p_z(k) for the current (l1, l2) pair, as in Python
            // This p_z is product over all features in union_f12
            Eigen::VectorXcd p_z_overall = Eigen::VectorXcd::Zero(n12_c);
            for (int k_idx = 0; k_idx < n12_c; k_idx++) {
                std::complex<double> prod_overall(1.0, 0.0);
                for (int feat_idx = 0; feat_idx < union_f12.size(); feat_idx++) {
                    int current_feat = union_f12(feat_idx);
                    prod_overall *= z_roots(k_idx) + w_matrix(l1, current_feat) * w_matrix(l2, current_feat);
                }
                p_z_overall(k_idx) = prod_overall;
            }

            for (int idx = 0; idx < union_f12.size(); idx++) {
                int j_feature = union_f12(idx); // Current feature to calculate contribution for
                
                Eigen::VectorXcd tmp_p_z_for_j = Eigen::VectorXcd::Zero(n12_c);
                for (int k_idx = 0; k_idx < n12_c; k_idx++) {
                    std::complex<double> denominator = z_roots(k_idx) + w_matrix(l1, j_feature) * w_matrix(l2, j_feature);
                    if (std::abs(denominator) < 1e-9) { // Avoid division by zero or near-zero
                        // This case needs careful handling; Python's NumPy might handle it implicitly.
                        // For now, let's assume it doesn't happen or set to a large number if p_z_overall is non-zero
                        // or zero if p_z_overall is zero. This depends on expected mathematical properties.
                        // A robust solution might involve limits or specific handling if z_roots(k_idx) is -w_matrix*w_matrix
                        tmp_p_z_for_j(k_idx) = 0; // Placeholder, needs careful thought if this occurs
                    } else {
                        tmp_p_z_for_j(k_idx) = p_z_overall(k_idx) / denominator;
                    }
                }
                
                double w_factor = w_matrix(l1, j_feature) * w_matrix(l2, j_feature) - 1.0;
                double contribution_val = complex_dot_v2(tmp_p_z_for_j, v_invc, n12);
                double final_contribution = w_factor * contribution_val * init_prediction_product;
                                
                if (l1 == l2) {
                    shap_value(i, j_feature) += final_contribution;
                } else {
                    shap_value(i, j_feature) += 2.0 * final_contribution;
                }
            }
        }
    }
}



Eigen::MatrixXd loss_treeshap(
        const Eigen::MatrixXd& x,
        const Eigen::VectorXd& y,
        const Rcpp::List& tree_summary,
        const Eigen::MatrixXcd& store_v_invc,
        const Eigen::MatrixXcd& store_z,
        const Eigen::MatrixXd& T0_x,  // Pre-computed SHAP values
        double learning_rate
) {
    // TreeSummary summary_tree_obj = list_to_tree_summary(tree_summary);


    int n_samples = x.rows();
    int n_features = T0_x.cols(); // Use T0_x_tree for dimensions
    
    Eigen::MatrixXd loss = Eigen::MatrixXd::Zero(n_samples, n_features);
    
    // Calculate T2 values for the current tree
    // The 'parallel' argument for T2 was not in the original signature shown, assuming default or not needed here.
    // If T2 is defined as: Eigen::MatrixXd T2(const Eigen::MatrixXd& x, const Rcpp::List& tree_summary, const Eigen::MatrixXcd& store_v_invc, const Eigen::MatrixXcd& store_z, bool parallel = false)
    // then call as:
    Eigen::MatrixXd T2_values = T2(x, tree_summary, store_v_invc, store_z, false); // Assuming T2 takes Rcpp::List

    Eigen::MatrixXd square_treeshap_term = T2_values * learning_rate * learning_rate;
    Eigen::MatrixXd scaled_T0_x_term = T0_x * learning_rate;
    
    for (int i = 0; i < n_samples; i++) {
        for (int j = 0; j < n_features; j++) {
            // y(i) is the residual for instance i, for the current tree
            // T0_x_tree(i, j) is the SHAP value for instance i, feature j, for the current tree
            loss(i, j) = square_treeshap_term(i, j) - 2 * y(i) * scaled_T0_x_term(i, j);
        }
    }
    
    return loss;
}