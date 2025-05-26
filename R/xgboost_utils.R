#' @importFrom xgboost xgb.model.dt.tree
NULL

# Loss implementation for xgboost model
#' @keywords internal
qshap_loss_xgboost <- function(explainer, x, y, y_mean_ori = NULL) {
  # Extract necessary components from explainer
  # max_depth <- explainer$max_depth
  model <- explainer$model
  store_v_invc <- explainer$store_v_invc
  store_z <- explainer$store_z
  base_score <- explainer$base_score
  xgb_trees <- explainer$trees
  
  # Get number of trees
  num_tree <- length(xgb_trees)
  
  # Initialize loss matrix
  loss <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  
  # Convert x to matrix if it's not already
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  

  # Does this work? Not looking at individual trees, but the whole model
  shap_cur <- predict(model, x, predcontrib = TRUE)
  shap_cur <- shap_cur[, -ncol(shap_cur)] # Remove bias col

  for (i in seq_len(num_tree)) {
    # Calculate residuals for current tree
    if (i == 1) {
      res <- y - base_score
    } else {
      pred_partial <- predict(model, x, iterationrange = c(1, i))
      res <- y - pred_partial
    }
    
    # Get summary tree for current tree
    summary_tree <- summarize_tree(xgb_trees[[i]])


    # if (i == 1) {
    #   T0_x <- shap_cur
    # }
    # else {
    #   T0_x <- shap_cur - shap_prev
    # }

    # shap_prev <- shap_cur

    # print(str(T0_x))

    tree_loss <- loss_treeshap(x, res, summary_tree, store_v_invc, store_z, shap_cur, 1.0)
    
    if (i == 1) {
      loss <- tree_loss
    } else {
      loss <- loss + tree_loss
    }
  }
  
  loss
}


# Formats an xgboost model into a list of simple_tree objects
#' @keywords internal
xgb_formatter <- function(xgb_model, max_depth) {
  # Extract tree information using xgb.model.dt.tree with use_int_id=TRUE for easier indexing
  tree_dt <- xgb.model.dt.tree(model = xgb_model, use_int_id = TRUE)
  
  # Get unique tree IDs
  tree_ids <- unique(tree_dt$Tree)
  
  # Initialize list to store simple_tree objects
  xgb_trees <- list()
  
  # Process each tree
  for (tree_id in tree_ids) {
    # Filter data for current tree
    tree_data <- tree_dt[tree_dt$Tree == tree_id, ]
    
    # Sort by Node to ensure proper ordering
    tree_data <- tree_data[order(tree_data$Node), ]
    
    # Get number of nodes
    node_count <- nrow(tree_data)
    
    # Initialize vectors for simple_tree structure
    children_left <- rep(-1L, node_count)
    children_right <- rep(-1L, node_count)
    feature <- rep(-1L, node_count)
    threshold <- rep(0.0, node_count)
    n_node_samples <- rep(0.0, node_count)  # Using Cover from xgboost
    value <- rep(0.0, node_count)  # Using Quality from xgboost
    
    # Create a mapping from Node ID to array index
    node_to_index <- setNames(seq_len(node_count) - 1L, tree_data$Node)  # 0-based indexing
    
    # Process each node
    for (i in seq_len(node_count)) {
      node_data <- tree_data[i, ]
      
      # Handle feature indices - convert feature names to indices if needed
      if (!is.na(node_data$Feature) && node_data$Feature != "Leaf") {
        # If Feature is a character (feature name), we need to convert to index
        if (is.character(node_data$Feature)) {
          # Get all unique feature names in the entire model (across all trees)
          all_features <- unique(tree_dt$Feature[tree_dt$Feature != "Leaf" & !is.na(tree_dt$Feature)])
          feature[i] <- match(node_data$Feature, all_features) - 1L  # 0-based indexing
        } else {
          # If it's already numeric, use it directly (but ensure 0-based)
          feature[i] <- as.integer(node_data$Feature)
        }
        
        threshold[i] <- as.numeric(node_data$Split)
      }
      
      # Set children indices using the node mapping
      if (!is.na(node_data$Yes)) {
        yes_node_id <- as.character(node_data$Yes)
        if (yes_node_id %in% names(node_to_index)) {
          children_left[i] <- node_to_index[yes_node_id]
        }
      }
      
      if (!is.na(node_data$No)) {
        no_node_id <- as.character(node_data$No)
        if (no_node_id %in% names(node_to_index)) {
          children_right[i] <- node_to_index[no_node_id]
        }
      }
      
      # Set node samples (using Cover)
      n_node_samples[i] <- as.numeric(node_data$Cover)
      
      # Set node value (using Quality)
      value[i] <- as.numeric(node_data$Quality)
    }
    
    # Create simple_tree object
    tree_obj <- simple_tree(
      children_left = children_left,
      children_right = children_right,
      feature = feature,
      threshold = threshold,
      max_depth = max_depth,
      n_node_samples = n_node_samples,
      value = value,
      node_count = node_count
    )
    
    # Add to list
    xgb_trees[[length(xgb_trees) + 1]] <- tree_obj
  }
  
  xgb_trees
}
