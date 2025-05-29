library(xgboost)
library(qshapr)

url <- "https://raw.githubusercontent.com/ageron/handson-ml2/master/datasets/housing/housing.csv"
housing <- read.csv(url)

housing$AveRooms <- housing$total_rooms / housing$households
housing$AveBedrms <- housing$total_bedrooms / housing$households
housing$AveOccup <- housing$population / housing$households

X <- data.frame(
    MedInc = housing$median_income,
    AveOccup = housing$AveOccup,
    Longitude = housing$longitude,
    Latitude = housing$latitude,
    HouseAge = housing$housing_median_age,
    AveRooms = housing$AveRooms,
    AveBedrms = housing$AveBedrms,
    Population = housing$population
)

y <- housing$median_house_value

model <- xgboost(
  data = as.matrix(X),
  label = y,
  nrounds = 50,
  max_depth = 2,
  verbose = 0,
)

explainer <- qshapr::create_tree_explainer(model)
rsq = qshapr::qshap_rsq(explainer, X, y, nsample = 1024)

print(rsq)