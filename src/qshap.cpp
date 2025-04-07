#include "qshap.h"

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
            
            // Find union_f12 (features where either l1 or l2 has value >= 1)
            std::vector<int> union_f12_vec;
            for (int j = 0; j < feature_uniq.size(); j++) {
                int feat = feature_uniq(j);
                if ((w_ind(l1, feat) + w_ind(l2, feat)) >= 1) {
                    union_f12_vec.push_back(feat);
                }
            }
            
            // Convert to Eigen::VectorXi for easier handling
            Eigen::Map<Eigen::VectorXi> union_f12(union_f12_vec.data(), union_f12_vec.size());
            
            int n12 = union_f12.size();
            // Begin to use the property of complex conjugate
            int n12_c = n12 / 2 + 1;
            
            // Extract relevant parts of store_v_invc and store_z
            Eigen::VectorXcd v_invc = store_v_invc.row(n12).head(n12_c);
            Eigen::VectorXcd z = store_z.row(n12).head(n12_c);
            
            Eigen::VectorXcd p_z = Eigen::VectorXcd::Zero(n12_c);
            Eigen::VectorXcd tmp_p_z = Eigen::VectorXcd::Zero(n12_c);
            
            // Calculate p_z values
            for (int k = 0; k < n12_c; k++) {
                std::complex<double> prod(1.0, 0.0);
                for (int idx = 0; idx < union_f12.size(); idx++) {
                    int j = union_f12(idx);
                    prod *= z(k) + w_matrix(l1, j) * w_matrix(l2, j);
                }
                p_z(k) = prod;
            }
            
            // Update only when feature j belongs to the union of f1 and f2
            for (int idx = 0; idx < union_f12.size(); idx++) {
                int j = union_f12(idx);
                
                // Remove the operation for j by dividing
                for (int k = 0; k < n12_c; k++) {
                    tmp_p_z(k) = p_z(k) / (z(k) + w_matrix(l1, j) * w_matrix(l2, j));
                }
                
                if (l1 != l2) {
                    shap_value(i, j) += 2.0 * (w_matrix(l1, j) * w_matrix(l2, j) - 1.0) * init_prediction_product * 
                        complex_dot_v2(tmp_p_z, v_invc, n12);
                } else {
                    shap_value(i, j) += (w_matrix(l1, j) * w_matrix(l2, j) - 1.0) * init_prediction_product * 
                        complex_dot_v2(tmp_p_z, v_invc, n12);
                }
            }
        }
    }
}

Eigen::MatrixXd T2(
    const Eigen::MatrixXd& x, 
    const TreeSummary& summary_tree,
    const Eigen::MatrixXcd& store_v_invc, 
    const Eigen::MatrixXcd& store_z,
    bool parallel
) {
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


Eigen::MatrixXd loss_treeshap(
        const Eigen::MatrixXd& x,
        const Eigen::VectorXd& y,
        const TreeSummary& summary_tree,
        const Eigen::MatrixXcd& store_v_invc,
        const Eigen::MatrixXcd& store_z,
        const Eigen::MatrixXd& T0_x,  // Pre-computed SHAP values
        double learning_rate
) {
    // Calculate T2 values and apply learning rate
    Eigen::MatrixXd square_treeshap_x = T2(x, summary_tree, store_v_invc, store_z) * 
        std::pow(learning_rate, 2);
    
    // Apply learning rate to pre-computed SHAP values
    Eigen::MatrixXd T0_x_scaled = T0_x * learning_rate;
    
    // Calculate the final result: square_treeshap_x - 2 * (y * T0_x.T).T
    // In Eigen, we need to handle the matrix multiplication carefully
    Eigen::MatrixXd res = square_treeshap_x - 2.0 * (y.asDiagonal() * T0_x_scaled);
    
    return res;
}