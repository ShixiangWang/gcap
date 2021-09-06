data <- readRDS("../ecDNA/data/train_data_v1.rds")

set.seed(2021)
idx <- c(
  sample(which(data[, "y_Circular"] == 1), size = 20),
  sample(which(data[, "y_Circular"] == 0), size = 2000)
)
ec <- data[idx, ]

ec <- data.table::as.data.table(ec)
ec$y_nonLinear <- NULL

usethis::use_data(ec, overwrite = TRUE)

# dir.create("inst/extdata", recursive = TRUE)
# save(ec, file = "inst/extdata/ec.RData")
