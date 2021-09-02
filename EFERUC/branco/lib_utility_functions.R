library(caret)
library(e1071)
library(doMC)
library(glmnet)
library(ranger)
library(rpart)
library(quadprog)

options(warn = 1)

# General utilities -----------------------------------------------------------

generate_cv_fold_indices <- function(n, folds){
  set.seed(13579)
  split(sample(1:n), rep(1:folds, length = n))
}

generate_train_test_split <- function(n, percent_split = 0.3){
  set.seed(23452)
  nsequence <- seq_len(n)
  sample_size <- floor(percent_split * n)
  test_indices <- sample(nsequence, size = sample_size)
  train_indices <- setdiff(nsequence, test_indices)
  
  list(train_indices = train_indices, test_indices = test_indices)
}

confirm_non_intersecting_indices <- function(xnames, ynames, dataset_name){
  if (length(intersect(xnames, ynames)) != 0){
    notification <- paste0("Shared samples between X and Y entities in ",
                           dataset_name, " dataset!")
    stop(notification)
  }
}

estimate_regressor_performance <- function(y_actual, y_predicted, spec = NULL, 
                                           dec = 3){

  mse       <- sum((y_actual - y_predicted)^2) / length(y_actual)
  rmse      <- sqrt(sum((y_actual - y_predicted)^2) / length(y_actual))
  rsquared  <- 1 - (sum((y_actual - y_predicted)^2) / 
                    sum((y_actual - mean(y_actual))^2))

  error_estimates <- c(mse = mse, rmse = rmse, rsquared = rsquared)
  error_estimates <- round(error_estimates, dec)

  # Function has been specialised for RMSE for these experiments
  if (!is.null(spec)){
    if (spec == "rmse"){
      round(rmse, dec)
    } else if (spec == "rsquared"){
      round(rsquared, dec)
    } else {
      stop("Specified metric may be defined but not directly returnable!")
    }
  } else {
    error_estimates
  }
}

# Statistical significance testing is performed on all unique pairs
# of column vectors in dframe.
perform_significance_testing <- function(dframe, dec = 3){
  entities <- colnames(dframe)
  entity_pairs <- combn(entites, 2)

  pvalue_estimates <- matrix(nrow = ncol(entity_pairs), ncol = 3)
  rnames <- vector(length = ncol(entity_pairs))
  for (i in seq_len(ncol(entity_pairs))){
    pair <- entity_pairs[, i]
    x_vect <- dframe[, pair[1]]
    y_vect <- dframe[, pair[2]]

    dcount <- length(which(x_vect > y_vect))
    signtest_pv <- binom.test(dcount, length(y_vect))$p.value
    pairedtt_pv <- t.test(x_vect, y_vect, paired = TRUE)$p.value
    wilcoxon_pv <- wilcox.test(x_vect, y_vect, paired = TRUE)$p.value

    pvalues <- c(signtest = signtest_pv, pairedtt = pairedtt_pv,
                 wilcoxon = wilcoxon_pv)
    pvalues <- round(pvalues, dec)
    pvalue_estimates[i] <- pvalues
    rnames[i] <- paste0(pair[1], "-", pair[2])
  }
  rownames(pvalue_estimates) <- rnames
  pvalue_estimates
}

# File I/0 --------------------------------------------------------------------

read_dataset_names <- function(){
  read.table(file = "../without_nominal_datasets/datasets.txt", 
    header = F, stringsAsFactors = F)[, 1]
}

read_experiment_data <- function(dataset){
  df <- data.matrix(read.csv(
    file = paste0("../without_nominal_datasets/", dataset, ".csv"),
    header = TRUE
  ))

  train_idx <- read.table(
    file = paste0("../without_nominal_datasets/splits/", dataset, "_train.txt"),
    header = FALSE,
    stringsAsFactors = FALSE)[,1]
  
  test_idx <- read.table(
    file = paste0("../without_nominal_datasets/splits/", dataset, "_test.txt"),
    header = FALSE,
    stringsAsFactors = FALSE)[,1]

  df_learning <- df[train_idx, ]
  df_test <- df[test_idx, ]

  x_learning  <- df_learning[, 2:ncol(df_learning)]
  ys_learning <- df_learning[, 1, drop = FALSE]

  x_test  <- df_test[, 2:ncol(df_test)]
  ys_test <- df_test[, 1, drop = FALSE]

  rownames(x_learning) <- train_idx
  rownames(ys_learning) <- train_idx
  rownames(x_test)  <- test_idx
  rownames(ys_test) <- test_idx

  # Sanity check
  confirm_non_intersecting_indices(
    rownames(x_learning), 
    rownames(x_test),
    dataset
  )
  confirm_non_intersecting_indices(
    rownames(ys_learning), 
    rownames(ys_test),
    dataset
  )

  list(dataset_name = dataset,
       x_learning   = x_learning, 
       x_test       = x_test, 
       ys_learning  = ys_learning,
       ys_test      = ys_test)
}

check_output_directories <- function(output_path, dataset_name){
  subdirs <- c(
    "weight_details/",
    "best_bin_size/", 
    "performance/regression/",
    "performance/regression/all_bin_performance/", 
    "performance/classification/",
    "performance/classification/RDS/"
  )

  dataset_path <- paste0(output_path, dataset_name, "/")
  for (subdir in subdirs){
    subdir_path <- paste0(dataset_path, subdir)
    if (!dir.exists(subdir_path)){
      dir.create(subdir_path, recursive = TRUE)
    }
  }
}

save_performance_results <- function(dataset_name,
                                     response,
                                     performance_file,
                                     new_matrix){
  if (file.exists(performance_file)){
    current_file <- read.csv(
      file = performance_file, 
      header = TRUE, 
      row.names = 1
    )

    if (response %in% rownames(current_file)){
      note <- paste0("Modifying existing error estimates for response: ", 
                      response, " in dataset: ", dataset_name, ".")
      message(note)
      updated_file <- current_file[!rownames(current_file) %in% response, ]
      updated_file <- rbind(updated_file, new_matrix)
      write.csv(updated_file, file = performance_file)
    } else {
      updated_file <- rbind(current_file, new_matrix)
      write.csv(updated_file, file = performance_file)
    }

  } else {
    write.csv(new_matrix, file = performance_file)
  }
}

log_regression_performance <- function(dataset_name, 
                                       response, 
                                       crpair,
                                       discretizer,
                                       crpair_weight_details,
                                       regression_performance){
  fname <- paste0(crpair[1], "-", crpair[2], "_", discretizer)
  weights_file <- paste0(
    OUTPUT_PATH, 
    dataset_name, 
    "/weight_details/", 
    fname, 
    ".txt"
  )
  weight_line <- c(response, crpair_weight_details$weights)
  write(
    x = weight_line, 
    ncolumns = length(weight_line), 
    append = TRUE, 
    file = weights_file
  )

  best_bin_size_file <- paste0(
    OUTPUT_PATH, 
    dataset_name, 
    "/best_bin_size/", 
    fname, 
    ".txt"
  )
  best_bin_line <- c(response, crpair_weight_details$best_single)
  write(
    x = best_bin_line, 
    ncolumns = 2, 
    append = TRUE, 
    file = best_bin_size_file
  )

  bins_performance_file <- paste0(
    OUTPUT_PATH, 
    dataset_name, 
    "/performance/regression/all_bin_performance/", 
    fname, 
    "_", response, ".csv"
  )
  bins_reg_matrix <- matrix(
    regression_performance$all_bin_performance, 
    nrow = 1, 
    dimnames = list(
      response = response, 
      names(regression_performance$all_bin_performance)
    )
  )
  write.csv(bins_reg_matrix, file = bins_performance_file)

  gen_performance_file <- paste0(
    OUTPUT_PATH, 
    dataset_name, 
    "/performance/regression/", 
    fname, 
    ".csv"
  )
  gen_reg_matrix <- matrix(
    regression_performance$general_performance, 
    nrow = 1, 
    dimnames = list(
      response = response, 
      names(regression_performance$general_performance)
    )
  )
  save_performance_results(
    dataset_name,
    response, 
    gen_performance_file,
    gen_reg_matrix
  )
}

log_classification_performance <- function(dataset_name,
                                           response,
                                           crpair,
                                           discretizer,
                                           classification_performance){
  fname <- paste0(crpair[1], "-", crpair[2], "_", discretizer)

  # Store full confusion matrix as RDS file
  base_cmpath <-  paste0(
    OUTPUT_PATH, 
    dataset_name, 
    "/performance/classification/RDS/", 
    fname
  )
  cms <- classification_performance$confusion_matrices
  for (i in 1:length(cms)){
    saveRDS(
      cms[[i]], 
      file = paste0(base_cmpath, "_", response, "_", i+1, ".RDS")
    )
  }

  # Save the accuracies to file by response
  performance_file <- paste0(
    OUTPUT_PATH, 
    dataset_name, 
    "/performance/classification/", 
    fname, 
    ".csv"
  )

  accuracies <- classification_performance$accuracy
  cls_matrix <- matrix(
    accuracies, 
    nrow = 1, 
    dimnames = list(
      response = response, 
      names(accuracies)
    )
  )
  save_performance_results(
    dataset_name,
    response, 
    performance_file,
    cls_matrix
  )
}

# model_type -> {classifier, regressor, agent, ...}
log_built_model <- function(response, model, model_type, fdir){
  filepath <- paste0(fdir, response, "_", model_type, ".rds")
  saveRDS(model, filepath)
}

# Learning Algorithms ---------------------------------------------------------
# Regressors

build_rf_regressor <- function(x, y){
  set.seed(45322)
  ranger(
    x = x, 
    y = y, 
    num.trees = 1000, 
    mtry = floor(ncol(x)/3), 
    verbose = FALSE
  )
}

build_xgb_regressor <- function(x, y){
  ctrl <- trainControl(
    method = "LGOCV", 
    p = 0.7,
    number = 1,
    allowParallel = FALSE
  )

  grid <- expand.grid(
    nrounds = 1500,
    max_depth = 6,
    eta = c(0.001, 0.01, 0.1, 0.2, 0.3),
    gamma = 0.05,
    subsample = 0.5,
    colsample_bytree = 1,
    min_child_weight = 1
  )

  # Having objective = "reg:squarederror" below causes the following
  # warning without any detriment to actual model building. 
  # "The following parameters were provided multiple times:
	#    objective
  #  Only the last value for each of them will be used."
  # If objective = "reg:squarederror" is not passed to train, the 
  # following warning is received:
  # "reg:linear is now deprecated in favor of reg:squarederror."
  # Therefore, the only option is to suppress local warnings within 
  # this function. - OIO, 31-Oct-2020 

  set.seed(47567)
  suppressWarnings(
    xgb_fit <- train(
      x = x, 
      y = y,
      method = "xgbTree",
      eval_metric = "rmse",
      verbose = 0, 
      trControl = ctrl,
      tuneGrid = grid,
      objective = "reg:squarederror"
    )
  )
  xgb_fit
}

build_lasso_regressor <- function(x, y){
  registerDoMC(cores = detectCores())
  set.seed(35342)
  lasso_model <- cv.glmnet(
    x, y,
    standardize = FALSE,
    alpha = 1,
    lambda = 10^seq(10, -3, length=100),
    parallel = TRUE
  )
  registerDoSEQ() # Disable parallelism

  lasso_model
}

build_ridge_regressor <- function(x, y){
  registerDoMC(cores = detectCores())
  set.seed(65462)
  ridge_model <- cv.glmnet(
    x, y,
    standardize = FALSE,
    alpha = 0,
    lambda = 10^seq(10, -3, length=100),
    parallel = TRUE
  )
  registerDoSEQ() # Disable parallelism

  ridge_model
}

# Classifiers

build_knn_classifier <- function(x, y){
  registerDoMC(cores = detectCores()) # Enable parallelism
  ctrl <- trainControl(
    method = "LGOCV", 
    p = 0.7,
    number = 1,
    allowParallel = TRUE
  )

  grid <- expand.grid(k = c(1, 3, 5, 7, 11, 13, 15))

  set.seed(74562)
  knn_fit <- train(
    x = x, 
    y = y,
    method = "knn",
    tuneGrid = grid,
    trControl = ctrl
  )
  registerDoSEQ() # Disable parallelism

  knn_fit
}

build_dt_classifier <- function(x, y){
  edata <- as.data.frame(cbind(x, y = y))
  rpart(y ~., data = edata, method = "class")
}

build_rf_classifier <- function(x, y){
  set.seed(45322)
  ranger(
    x = x, 
    y = y, 
    num.trees = 1000, 
    probability = TRUE,
    verbose = FALSE
  )
}

build_nb_classifier <- function(x, y){
  naiveBayes(x, y)
}

convex_lr <- function(x, y){
  Rinv <- solve(chol(t(x) %*% x))
  C <- cbind(rep(1,ncol(Rinv)), diag(ncol(Rinv)))
  b <- c(1, rep(0, ncol(Rinv)))
  d <- t(y) %*% x
  qsolve <- solve.QP(
    Dmat = Rinv, 
    factorized = TRUE, 
    dvec = d, 
    Amat = C, 
    bvec = b, 
    meq = 1
  )
  qsolve$solution
}

vscale <- function(y){(y-min(y))/(max(y)-min(y))}

sconvex_lr <- function(x, y){
  y <- vscale(y)
  x <- apply(x, 2, FUN = function(ny) vscale(ny))

  Rinv <- solve(chol(t(x) %*% x))
  C <- cbind(rep(1,ncol(Rinv)), diag(ncol(Rinv)))
  b <- c(1, rep(0, ncol(Rinv)))
  d <- t(y) %*% x
  qsolve <- solve.QP(
    Dmat = Rinv, 
    factorized = TRUE, 
    dvec = d, 
    Amat = C, 
    bvec = b, 
    meq  = 1
  )
  qsolve$solution
}