#' Create a Tree Explainer from an rpart model
#' 
#' This function takes an rpart decision tree model and creates a TreeExplainer
#' object that can be used for explanation.
#' 
#' @param tree_model An rpart model object
#' @return A TreeExplainer object
#' @export
create_tree_explainer <- function(tree_model) {
  if (!inherits(tree_model, "rpart")) {
    stop("tree_model must be an rpart object")
  }
  
  # Create the TreeExplainer
  explainer <- create_tree_explainer_cpp(tree_model)
  
  # Add class for method dispatch
  class(explainer) <- c("qshapr_tree_explainer", class(explainer))
  
  return(explainer)
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