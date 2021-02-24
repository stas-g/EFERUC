source("lib_utility_functions.R")

classifiers  <- c("NB", "RF", "DT", "KNN")
regressors   <- c("LS", "RR", "RF", "XGB")
CRPAIRS <- expand.grid(
  classifier = classifiers, 
  regressor = regressors,
  stringsAsFactors = FALSE
)
PMETRIC <- "rsquared"
ICVFOLDS <- 5 
MAX_BINS <- 10
DISCRETIZERS <- c("frequency")
OUTPUT_PATH <- "../output/"
DATASET_NAMES <- read_dataset_names()