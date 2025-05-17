library(xgboost)
library(qshapr)

X <- matrix(1:10, ncol=1)
y <- rep(0, 10)
model <- xgboost(data = X, label = y, nrounds = 1, verbose = 0)
print(class(model))

explainer <- qshapr::create_tree_explainer(model)