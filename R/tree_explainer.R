#' @importFrom xgboost xgb.parameters
NULL

#' Create a QSHAPR Tree Explainer
#' 
#' Add this later
#' 
#' @param tree_model A tree model object
#' @return A TreeExplainer object
#' @export
create_tree_explainer <- function(tree_model) {
  if (inherits(tree_model, "xgb.Booster")) {
    # Setting defauult depth to 6 to match Python package
    # max_depth_xgb <- xgb.parameters(tree_model)$max_depth

    # if (is.null(max_depth_xgb)) {
    #   max_depth_xgb <- 6
    # }

    # TODO: Add support for optinoal base score
    base_score <- 0.5
    # base_score <- xgb.parameters(tree_model)$base_score

    xgb_trees <- xgb_formatter(tree_model, max_depth_xgb)

    explainer <- create_xgboost_explainer()

    class(explainer) <- c("qshapr_tree_explainer", class(explainer))
    return(explainer)
  }
  else {
    stop("tree_model must be a supported model object")
  }
  
}

#' Get a summary of the tree structure
#' 
#' @param explainer A tree explainer object
#' @return A list containing the tree structure
#' @export
get_summary <- function(explainer) {
  if (!inherits(explainer, "qshapr_tree_explainer")) {
    stop("explainer must be a qshapr_tree_explainer object")
  }
  
  # Get the summary from C++
  return(get_tree_summary(explainer))
} 