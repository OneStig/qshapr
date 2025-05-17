#' @importFrom xgboost xgb.model.dt.tree
#' @importFrom data.table split.data.table
NULL

# Constructor for a simple_tree object, similar to a Python dataclass
simple_tree <- function(children_left,
                        children_right,
                        feature,
                        threshold,
                        max_depth,
                        n_node_samples,
                        value,
                        node_count) {
  structure(list(
    children_left   = children_left,
    children_right  = children_right,
    feature         = feature,
    threshold       = threshold,
    max_depth       = max_depth,
    n_node_samples  = n_node_samples,
    value           = value,
    node_count      = node_count
  ), class = "simple_tree")
}

# Formats an xgboost model into a list of simple_tree objects
#' @keywords internal
xgb_formatter <- function(xgb_model, max_depth) {
  # Extract node-level data for every tree
  dt <- xgb.model.dt.tree(model = xgb_model)
  trees <- split(dt, by = "Tree", keep.by = FALSE)

  lapply(trees, function(tree) {
    nodes <- tree$Node
    n     <- length(nodes)

    left_idx  <- match(tree$Yes, nodes) - 1L
    right_idx <- match(tree$No,  nodes) - 1L

    feat <- ifelse(tree$Feature == "Leaf",
                   -2L,
                   as.integer(sub("^f", "", tree$Feature)))

    thresh <- ifelse(tree$Feature == "Leaf",
                     NA_real_,
                     tree$Split)

    samples <- tree$Cover
    vals    <- ifelse(tree$Feature == "Leaf",
                      tree$Quality,
                      NA_real_)

    simple_tree(
      children_left   = left_idx,
      children_right  = right_idx,
      feature         = feat,
      threshold       = thresh,
      max_depth       = max_depth,
      n_node_samples  = samples,
      value           = vals,
      node_count      = n
    )
  })
}
