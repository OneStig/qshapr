// #include "tree_explainer.h"
// #include <Rcpp.h>

// TreeExplainer::TreeExplainer() {
//     // remove this stub later
// }

// TreeExplainer::TreeExplainer(const Rcpp::List& tree_model) {
    
//     // Get the frame data.frame from the rpart object
//     Rcpp::DataFrame frame = Rcpp::as<Rcpp::DataFrame>(tree_model["frame"]);
    
//     // Extract var column which contains splitting variables (or "<leaf>" for leaf nodes)
//     Rcpp::CharacterVector var = frame["var"];
//     int node_count = var.size();
    
//     // Get yval (predicted values) and n (number of samples per node)
//     Rcpp::NumericVector yval = frame["yval"];
//     Rcpp::IntegerVector n = frame["n"];
    
//     // Get control parameters
//     Rcpp::List control = tree_model["control"];
//     int max_depth = Rcpp::as<int>(control["maxdepth"]);
    
//     // We need to extract splits information to get feature indices and thresholds
//     Rcpp::NumericMatrix splits = Rcpp::as<Rcpp::NumericMatrix>(tree_model["splits"]);
    
//     // Initialize vectors for TreeSummary
//     Eigen::VectorXi children_left = Eigen::VectorXi::Constant(node_count, -1);
//     Eigen::VectorXi children_right = Eigen::VectorXi::Constant(node_count, -1);
//     Eigen::VectorXi feature = Eigen::VectorXi::Constant(node_count, -1);
//     Eigen::VectorXd threshold = Eigen::VectorXd::Zero(node_count);
//     Eigen::VectorXd sample_weight = Eigen::VectorXd::Ones(node_count);
//     Eigen::VectorXd init_prediction = Eigen::VectorXd::Zero(node_count);
    
//     // Build the tree structure
//     // In rpart, we need to build parent-child relationships
//     std::vector<int> split_index(node_count, -1);
    
//     // First pass: identify split nodes and their features/thresholds
//     int split_idx = 0;
//     for (int i = 0; i < node_count; ++i) {
//         std::string var_name = Rcpp::as<std::string>(var[i]);
//         if (var_name != "<leaf>") {
//             // This is a split node
//             feature[i] = split_idx; // temporary storing the split index
//             if (split_idx < splits.nrow()) {
//                 threshold[i] = splits(split_idx, 3); // index 3 contains the split threshold
//                 split_index[i] = split_idx++;
//             }
//         }
//     }
    
//     // Second pass: establish parent-child connections
//     // rpart nodes follow a depth-first order where each split node is followed by its left child subtree, then right child subtree
//     std::vector<int> stack;
//     stack.push_back(0); // Start with the root node
    
//     for (int i = 1; i < node_count; ++i) {
//         while (!stack.empty()) {
//             int parent = stack.back();
//             if (children_left[parent] == -1) {
//                 // Left child not assigned yet
//                 children_left[parent] = i;
//                 break;
//             } else if (children_right[parent] == -1) {
//                 // Right child not assigned yet
//                 children_right[parent] = i;
//                 break;
//             } else {
//                 // Both children assigned, pop this node
//                 stack.pop_back();
//             }
//         }
        
//         // If not a leaf, add to stack
//         if (Rcpp::as<std::string>(var[i]) != "<leaf>") {
//             stack.push_back(i);
//         }
//     }
    
//     // Third pass: map feature indices to actual feature indices
//     // For simplicity, we'll map internal split indices to sequential indices (0, 1, 2, ...)
//     // and build feature_uniq vector
//     std::map<std::string, int> feature_map;
//     int next_feature_idx = 0;
    
//     for (int i = 0; i < node_count; ++i) {
//         std::string var_name = Rcpp::as<std::string>(var[i]);
//         if (var_name != "<leaf>") {
//             if (feature_map.find(var_name) == feature_map.end()) {
//                 feature_map[var_name] = next_feature_idx++;
//             }
//             feature[i] = feature_map[var_name];
//         }
//     }
    
//     // Create feature_uniq vector
//     Eigen::VectorXi feature_uniq(feature_map.size());
//     for (const auto& pair : feature_map) {
//         feature_uniq[pair.second] = pair.second;
//     }
    
//     // Calculate init_prediction and sample_weight as in the Python implementation
//     double n_root = n[0];
    
//     // Define helper function for traversal, using C++11 lambda
//     std::function<void(int)> traversal_summarize_tree = [&](int v) {
//         int v_l = children_left[v];
//         int v_r = children_right[v];
        
//         double n_v = n[v];
        
//         if (v_l < 0) {  // leaf
//             init_prediction[v] = yval[v] * n_v / n_root;
//         } else {
//             double n_l = n[v_l];
//             double n_r = n[v_r];
            
//             sample_weight[v_l] = n_v / n_l;
//             sample_weight[v_r] = n_v / n_r;
            
//             traversal_summarize_tree(v_l);
//             traversal_summarize_tree(v_r);
//         }
//     };
    
//     // Start traversal from the root
//     traversal_summarize_tree(0);
    
//     // Assign values to TreeSummary struct
//     summary.children_left = children_left;
//     summary.children_right = children_right;
//     summary.feature = feature;
//     summary.feature_uniq = feature_uniq;
//     summary.threshold = threshold;
//     summary.max_depth = max_depth;
//     summary.sample_weight = sample_weight;
//     summary.init_prediction = init_prediction;
//     summary.node_count = node_count;
    
//     // Store precomputed store_v_invc, store_z
//     store_v_invc = store_complex_v_invc(max_depth * 2);
//     store_z = store_complex_root(max_depth * 2);
// }

// void TreeExplainer::initialize_xgboost(const int max_depth, const long double base_score, const std::vector<SimpleTree> xgb_trees) {
//     this->model_type = TreeType::XGBOOST;

//     this->max_depth = max_depth;
//     this->base_score = base_score;
//     this->xgb_trees = xgb_trees;
// }

// Eigen::MatrixXd TreeExplainer::get_loss(
//     const Eigen::MatrixXd& x,
//     const Eigen::VectorXd& y,
//     const Eigen::MatrixXd& T0_x,
//     double y_mean_ori
// ) {
//     int max_depth = this->summary.max_depth;
//     const Eigen::MatrixXcd& store_v_invc = this->store_v_invc;
//     const Eigen::MatrixXcd& store_z = this->store_z;

//     // If y_mean_ori is not provided (NaN), calculate it
//     if (std::isnan(y_mean_ori)) {
//         y_mean_ori = y.mean();
//     }

//     // Call loss_treeshap with the parameters
//     Eigen::MatrixXd loss = loss_treeshap(x, y, this->summary, store_v_invc, store_z, T0_x);
    
//     return loss;
// }

// Eigen::VectorXd TreeExplainer::get_rsq(
//         const Eigen::MatrixXd& x,
//         const Eigen::VectorXd& y,
//         const Eigen::MatrixXd& T0_x,
//         bool loss_out,
//         int ncore,
//         int nsample,
//         double nfrac,
//         int random_state
// ) {
//     Eigen::MatrixXd x_sampled = x;
//     Eigen::VectorXd y_sampled = y;
    
//     // Handle sample size (nsample)
//     if (nsample > 0) {
//         if (nsample <= 0 || nsample >= x.rows()) {
//             // TODO: Error handling for invalid nsample
//             // "Sampling sample size must be larger than 0 and smaller than the total number of samples"
//         }
        
//         // TODO: Implement random sampling
//         // Set random seed to random_state
//         // Sample nsample rows from x and y
//         // For now we'll just use the first nsample rows
//         x_sampled = x.topRows(nsample);
//         y_sampled = y.head(nsample);
//     }
    
//     // Handle sample fraction (nfrac)
//     if (nfrac > 0.0 && nsample <= 0) {
//         if (nfrac <= 0.0 || nfrac >= 1.0) {
//             // TODO: Error handling for invalid nfrac
//             // "Sample fraction must be between (0, 1)"
//         }
        
//         // Calculate sample size from fraction
//         nsample = static_cast<int>(x.rows() * nfrac);
        
//         // TODO: Implement random sampling
//         // Set random seed to random_state
//         // Sample nsample rows from x and y
//         // For now we'll just use the first nsample rows
//         x_sampled = x.topRows(nsample);
//         y_sampled = y.head(nsample);
//     }
    
//     double y_mean_ori = y_sampled.mean();
    
    
//     double sst = (y_sampled.array() - y_mean_ori).square().sum();
    
//     Eigen::MatrixXd loss_matrix = get_loss(x_sampled, y_sampled, T0_x, y_mean_ori);
    
//     Eigen::VectorXd rsq_values = -loss_matrix.colwise().sum() / sst;
    
//     return rsq_values;
// }