## code to prepare `toy_model` dataset goes here

data("ec")

library(xgboost)

data <- xgb.DMatrix(as.matrix(ec[, -"y_Circular"]), label = ec$y_Circular)

toy_model <- xgboost(data,
  eta = 0.3,
  nthread = 10,
  nrounds = 100,
  objective = "binary:logistic",
  verbose = 2,
  params = list(
    scale_pos_weight = 100,
    eval_metric = "aucpr"
  )
)

predict(toy_model, as.matrix(ec[, -"y_Circular"]))

saveRDS(toy_model, file = "inst/extdata/toy_model.rds")
