library(arules) # For discretize(), used in unsupervised_discretization()
source("lib_experiment_parameters.R")

# Global constants ------------------------------------------------------------
# Mostly read in from lib_experiment_parameters.R

regressor_functions  <- list("LS"  = build_lasso_regressor,
                             "RR"  = build_ridge_regressor,
                             "RF"  = build_rf_regressor,
                             "XGB" = build_xgb_regressor)
classifier_functions <- list("RF"  = build_rf_classifier,
                             "NB"  = build_nb_classifier,
                             "DT"  = build_dt_classifier,
                             "KNN" = build_knn_classifier)

# Core functions --------------------------------------------------------------

unsupervised_discretization <- function(y, bin_size, type){
  # Will produce a warning if the method cannot discretize the y vector into
  # the specified bin size because of how it works.
  # In this situation, the maximum number of bins are returned.
  idx <- discretize(y, method = type, breaks = bin_size)
  ybins <- split(y, idx)
  names(ybins) <- 0:(length(ybins)-1)
  for (ybin in ybins){
    if (length(ybin) < 50){
      return(NA)
    }
    if (length(unique(ybin)) < 2) {
      return(NA) 
    }
  }
  ybins 
}

get_layered_bins <- function(y, discretizer, max_bins = NULL){
  layered_bins <- list()
  if (is.null(max_bins)){
    for (bin_size in 1:MAX_BINS){
      temp_bins <- tryCatch(unsupervised_discretization(
                              y        = y, 
                              bin_size = bin_size, 
                              type     = discretizer
                            ), 
                            error = function(e) e,
                            warning = function(w) w        
                  )

      if (class(temp_bins)[1] == "simpleWarning"){
        break
      
      } else if (is.na(temp_bins[1])){
        break

      } else {
        layered_bins[[bin_size]] <- temp_bins
      }
    }
  
  } else {
    for (bin_size in 1:max_bins){
      temp_bins <- unsupervised_discretization(
        y        = y, 
        bin_size = bin_size, 
        type     = discretizer
      )
      if (is.na(temp_bins[1])) {
        break     
      } else {
        layered_bins[[bin_size]] <- temp_bins
      }
    }
  }
  
  layered_bins
}

# Learning functions ----------------------------------------------------------

build_regressors <- function(x, ybinned, bin_size, regressor){
  learner <- regressor_functions[[regressor]]
  models  <- vector(mode = "list", length = bin_size)
  for (i in 1:bin_size){
    y <- ybinned[[i]]
    x_train <- x[names(y), ]
    models[[i]] <- learner(x_train, y)
  }
  models
}

build_classifier <- function(x, y, classifier){
  learner <- classifier_functions[[classifier]]
  learner(x, y)
}

upsample_dataset <- function(x, ybinned){
  yclass <- NULL
  for (i in 1:length(ybinned)){
    ybin <- ybinned[[i]]
    ybinclass <- rep(i, length(ybin))
    names(ybinclass) <- names(ybin)
    yclass <- c(yclass, ybinclass)
  }
  yclass <- as.factor(yclass)
  xclass <- x[names(yclass), ]

  set.seed(98328)
  upSample(x = xclass, y = yclass, list = TRUE)
}

build_models <- function(x, layered_ybins, max_bins, crpair){
  crpair_models <- vector(mode = "list", length = max_bins)
  for (bin_size in 1:max_bins){
    ybinned <- layered_ybins[[bin_size]]
    if (bin_size == 1){ # For base model.
      # Build only base regressor as there is nothing to classify.
      # regressors should have exactly one model.
      regressors <- build_regressors(x, ybinned, bin_size, crpair$regressor)
      crpair_models[[bin_size]]$regressor <- regressors[[1]]

    } else {
      # Upsample training set. Then build single classifier.
      udf <- upsample_dataset(x, ybinned)
      classifier <- build_classifier(udf$x, udf$y, crpair$classifier)
      # Build regressors for each bin.
      regressors <- build_regressors(x, ybinned, bin_size, crpair$regressor)
    
      crpair_models[[bin_size]]$classifier <- classifier
      crpair_models[[bin_size]]$regressors <- regressors
    }
  }

  crpair_models
}

get_regression_predictions <- function(model, regressor, x_new){
  # There's a need to add special checks because of how the different
  # learning algorithms output predictions. Some output a vector
  # while others output a matrix.      
  if (regressor == "RF"){
    predict(model, x_new)$predictions
  } else if (regressor == "RR" || regressor == "LS"){
    predict(model, s = model$lambda.min, newx = x_new)[,1]
  } else { # For caret based models.
    # Note that models created using caret might either return a vector
    # or a matrix when model$finalModel is used, as is the case with
    # XGBoost. However, if the direct model is used, as in 
    # "model" in the predict function, this simply returns a vector,
    # which is what we want.
    predict(model, x_new)
  }
}

get_classification_probabilities <- function(model, classifier, x_new){
  if (classifier == "RF"){
    predict(model, x_new)$predictions
  
  } else if (classifier == "NB"){
    predict(model, x_new, type = "raw")
  
  } else if (classifier == "DT"){
    predict(model, as.data.frame(x_new), type = "prob")
  
  } else if (classifier == "KNN") {
    predict(model, x_new, type = "prob")
  
  } else {
    notification <- paste0("Non-existent classifier type has been passed ",
                           "to get_classifier_probabilities().")
    stop(notification)
  }
}

weight_regs_with_probs <- function(xreg_preds, xclass_probs){
  rowSums(xreg_preds*xclass_probs)
}

make_federated_predictions <- function(crpair_models, x_new, max_bins, crpair){
  crpair_predictions <- matrix(nrow = nrow(x_new), ncol = max_bins)
  colnames(crpair_predictions) <- 1:max_bins
  for (bin_size in 1:max_bins){
    bin_size_models <- crpair_models[[bin_size]]

    if (bin_size == 1){ # Base case, so there is only a regressor.
      xreg_pred <- get_regression_predictions(
                      model = bin_size_models$regressor, 
                      regressor = crpair$regressor, 
                      x_new = x_new
                    )
      crpair_predictions[, bin_size] <- xreg_pred

    } else {
      xclass_probs <- get_classification_probabilities(
                        model = bin_size_models$classifier,
                        classifier = crpair$classifier,
                        x_new = x_new    
                      )
      xreg_preds <- sapply(
                      bin_size_models$regressors, 
                      FUN = function(x){
                              get_regression_predictions(
                                model = x,
                                regressor = crpair$regressor,
                                x_new = x_new
                              )
                            }
                    )
      
      xreg_pred <- weight_regs_with_probs(xreg_preds, xclass_probs)
      crpair_predictions[, bin_size] <- xreg_pred
    }
  }
  crpair_predictions  
}

learn_bin_size_weights <- function(x_learning, 
                                   y_learning, 
                                   discretizer,
                                   max_bins, 
                                   crpair, 
                                   fold_indices){
  # The way the x and y are constructed for crpair could be slow
  # for large enough datasets, but this should not be a problem in 
  # our evaluation.
  x_crpair <- NULL
  y_crpair <- NULL
  nfolds <- length(fold_indices)
  for (fold_index in 1:nfolds){
    message(paste0("     Fold ", fold_index, " of ", nfolds))
    validation_idx <- fold_indices[[fold_index]]
    x_train <- x_learning[-validation_idx, ]
    x_val   <- x_learning[validation_idx, ]
    y_train <- y_learning[-validation_idx]
    y_val   <- y_learning[validation_idx]
 
    layered_ytrain_bins <- get_layered_bins(y_train, discretizer, max_bins) 

    # Model building might fail because there are too few samples for a 
    # given bin size.
    built_models  <- tryCatch(
      build_models(
        x_train, 
        layered_ytrain_bins, 
        max_bins, 
        crpair
      ),
      error = function(e) NA,
      warning = function(w) NA
    )
    
    if (is.na(built_models[1])){ 
      message("       Model building in weight learning failed!")
      next
    }
    
    crpair_predictions  <- make_federated_predictions(
      crpair_models = built_models, 
      x_new = x_val,
      max_bins = max_bins, 
      crpair = crpair
    )
    rownames(crpair_predictions) <- rownames(x_val)

    x_crpair <- rbind(x_crpair, crpair_predictions)
    y_crpair <- c(y_crpair, y_val)
  }

  if (is.null(x_crpair)){
    message("       No input dataset(X) for internal CV generated.")
    return(NA) 
  }
  # Additional sanity check
  x_crpair <- as.matrix(x_crpair[names(y_crpair), ])

  # Learn weights using convex linear regression
  crpair_weights <- tryCatch(
      round(convex_lr(x_crpair, y_crpair), 2),
      error = function(e) NA,
      warning = function(w) NA
  )
  if (is.na(crpair_weights[1])){
    crpair_weights <- round(sconvex_lr(x_crpair, y_crpair), 2)
  }

  # Get individual best performing bin size.
  bin_size_rsq <- apply(x_crpair, 2, FUN = function(y){
                                             estimate_regressor_performance(
                                             y_crpair, y, "rsquared"
                                             )
                                           })
  best_bin_size <- as.vector(which(bin_size_rsq == max(bin_size_rsq)))[1]

  list(weights = crpair_weights, best_single = best_bin_size)
}

class_test_samples <- function(y_test, layered_ybins, max_bins){
  y_classifications <- vector(
    mode = "list", 
    length = length(layered_ybins)
  )

  # Base case sample IDs. Used in get_classification_performance().
  y_classifications[[1]] <- names(y_test) 
  for (bin_size in 2:max_bins){
    ybinned <- layered_ybins[[bin_size]]
    bslist <- vector(mode = "list", length = bin_size)
    for (bin in 1:bin_size){
      bin_values <- ybinned[[bin]]
      imin <- min(bin_values)
      imax <- max(bin_values)
      
      if (bin == bin_size){
        csamps <- names(which((y_test >= imin) & (y_test <= imax)))
        bslist[[bin]] <- csamps
      } else {
        csamps <- names(which((y_test >= imin) & (y_test < imax)))
        bslist[[bin]] <- csamps
      }
    }
    y_classifications[[bin_size]] <- bslist
  }

  y_classifications
}

# Final performance estimators ------------------------------------------------

get_final_regression_performance <- function(crpair_weight_details, 
                                             crpair_predictions, 
                                             y_test){
  # 1. Get weighted prediction
  weighted_prediction <- rowSums(crpair_predictions %*% 
                                 diag(crpair_weight_details$weights)
                                )
  # 2. Get best prediction  
  best_bin_prediction <- crpair_predictions[, crpair_weight_details$best_single]

  # 3. Estimate prediction error  
  base_performance <- estimate_regressor_performance(
    y_actual = y_test,
    y_predicted = crpair_predictions[, 1], # 1 because that's the base case.
    spec = PMETRIC
  )
  best_bin_performance <- estimate_regressor_performance(
    y_actual = y_test, 
    y_predicted = best_bin_prediction,
    spec = PMETRIC
  )
  lr_weighted_performance <- estimate_regressor_performance(
    y_actual = y_test, 
    y_predicted = weighted_prediction,
    spec = PMETRIC
  )
  # Get performance for all bins
  all_bin_performance <- apply(
    crpair_predictions, 
    2, 
    FUN = function(y){
            estimate_regressor_performance(
              y_test, y, "rsquared"
            )
          }
  )
  
  general_performance <- c(
    base = base_performance, 
    best_bin = best_bin_performance,
    lr_weighted = lr_weighted_performance
  )
  list(
    general_performance = general_performance,
    all_bin_performance = all_bin_performance
  )                                          
}

get_classification_performance <- function(built_models,
                                           max_bins,
                                           crpair,
                                           x_test, 
                                           y_test_classes){

  accuracies <- vector(length = max_bins - 1)
  confusion_matrices <- vector(mode = "list", length = max_bins - 1)
  for (bin_size in 2:max_bins){
    y_actual <- vector()
    bs_ytest_classes <- y_test_classes[[bin_size]] 
    for (bin in 1:bin_size){
      for (sample_name in bs_ytest_classes[[bin]]){
        y_actual[sample_name] <- bin
      }
    }

    x_test_bs <- x_test[names(y_actual), ]
    yclass_probs <- get_classification_probabilities(
                      model = built_models[[bin_size]]$classifier,
                      classifier = crpair$classifier,
                      x_new = x_test_bs    
                    )
    y_predicted <- apply(yclass_probs, 1, which.max)

    y_actual    <- as.factor(y_actual)
    y_predicted <- factor(y_predicted, levels = levels(y_actual))
  
    cf <- confusionMatrix(y_predicted, y_actual)
    accuracy <- cf$overall["Accuracy"]
    
    accuracies[bin_size - 1] <- accuracy
    confusion_matrices[[bin_size - 1]] <- cf
  }
  names(accuracies) <- paste0("BS", 2:max_bins)
  accuracies <- round(accuracies, 2)

  list(accuracy = accuracies, confusion_matrices = confusion_matrices)
}

get_stable_max_bins <- function(y_learning, discretizer, fold_indices){
  stable_max_bins <- MAX_BINS
  nfolds <- length(fold_indices)
  for (fold_index in 1:nfolds){
    message(paste0("  Validating Fold ", fold_index, " of ", nfolds))
    validation_idx <- fold_indices[[fold_index]]
    y_train <- y_learning[-validation_idx]  
    layered_ytrain_bins <- get_layered_bins(y_train, discretizer, MAX_BINS) 

    nmaxbins <- length(layered_ytrain_bins) 
    if (nmaxbins < MAX_BINS){
      note <- paste0("   Internal CV max bins less than specified MAX BINS: ", 
                     nmaxbins)
      message(note) 
      if (nmaxbins < stable_max_bins){
        stable_max_bins <- nmaxbins 
      } 
    } else if (nmaxbins == MAX_BINS){ 
      message("   Bin selection is stable") 
    }
  }
  stable_max_bins
}

# Main loop -------------------------------------------------------------------

perform_feruc <- function(dataset_name, discretizer){
  expdata <- read_experiment_data(dataset_name)
  binlog_path <- paste0("../output/", dataset_name, "/bin_logs.txt")
  if (file.exists(binlog_path)){
    file.remove(binlog_path)
  }

  for(response in colnames(expdata$ys_learning)){
    message(paste0("Target: ", response))
    y_learning <- expdata$ys_learning[, response]
    names(y_learning) <- rownames(expdata$ys_learning)
    icv_fold_indices <- generate_cv_fold_indices(length(y_learning), ICVFOLDS)
    
    # Check if MAX_BINS is stable for internal cross-validation
    max_bins <- get_stable_max_bins(y_learning, discretizer, icv_fold_indices)
    if (max_bins == 1){   
      message(paste0("  Cannot discretize, max bins is: ", max_bins))  
      next
    }
    layered_ylearning_bins <- get_layered_bins(y_learning, discretizer, max_bins)
        
    message(paste0("  \nLearning Max Bins: ", max_bins))  
    for (crpair_index in 1:nrow(CRPAIRS)) {
      crpair <- CRPAIRS[crpair_index, ]
       # Log the specified MAX BINS and the actual discretizer max bins
      bin_note <- paste0(crpair$classifier, "-", crpair$regressor, 
                         "\t", discretizer, "\t", response, "\t", MAX_BINS, 
                         " - ", max_bins)
       
      write(bin_note, file = binlog_path, append = TRUE)
      message(paste0("  CR: ", crpair[1], "-", crpair[2], " | ", response))

      ## Step 1
      # Learn weights for each classifier-regressor pair (crpair)
      message("   Learning weights")
      crpair_weight_details  <- learn_bin_size_weights(
        x_learning = expdata$x_learning,
        y_learning = y_learning, 
        discretizer = discretizer,
        max_bins = max_bins,
        crpair = crpair,
        fold_indices = icv_fold_indices
      )
      
      if (is.na(crpair_weight_details[1])){ 
        message("    Weight learning failed!")
        next 
      }

      ## Build models using full learning set and make test set predictions
      message("   Building models")
      built_models <- build_models(
        x = expdata$x_learning,
        layered_ybins = layered_ylearning_bins,
        max_bins = max_bins,
        crpair = crpair
      )

      crpair_predictions  <- make_federated_predictions(
        crpair_models = built_models,
        x_new = expdata$x_test, 
        max_bins = max_bins,
        crpair = crpair
      )

      ## Step 2. Get performance results
      message("   Estimating performance")
      y_test <- expdata$ys_test[, response]
      # Get regression performance
      regression_performance <- get_final_regression_performance(
        crpair_weight_details = crpair_weight_details, 
        crpair_predictions = crpair_predictions, 
        y_test = y_test
      )

      # Get classification performance
      y_test_classes <- class_test_samples(
        y_test = y_test, 
        layered_ybins = layered_ylearning_bins, 
        max_bins = max_bins
      )
      classification_performance <- get_classification_performance(
        built_models = built_models,
        max_bins = max_bins,
        crpair = crpair,
        x_test = expdata$x_test,
        y_test_classes = y_test_classes
      )

      ## Step 3
      # Log performance.
      # Not logging models because given the scale of the experiments, storing
      # all generated models will be prohibitive space-wise. It is also not
      # clear if potential readers of the paper will be interested in using 
      # models directly. However, we will be providing everything needed for
      # an independent party to both reproduce and replicate our results in
      # copious amounts of detail.
      message("   Logging performance results")
      log_regression_performance(
        dataset_name = dataset_name, 
        response = response, 
        crpair = crpair,
        discretizer = discretizer,
        crpair_weight_details = crpair_weight_details,
        regression_performance = regression_performance
      )

      log_classification_performance(
        dataset_name = dataset_name,
        response = response,
        crpair = crpair,
        discretizer = discretizer,
        classification_performance = classification_performance
      )
      message("")
    }
  }
}

main <- function(dataset_index){
  if (!is.na(dataset_index)){
    dataset_name <- DATASET_NAMES[dataset_index]
    check_output_directories(OUTPUT_PATH, dataset_name)

    for (discretizer in DISCRETIZERS){
      message(paste0("Dataset: ", dataset_name, " - Discretizer: ", discretizer))
      perform_feruc(dataset_name, discretizer)
    }

  } else {
    for (dataset_name in DATASET_NAMES){
      check_output_directories(OUTPUT_PATH, dataset_name)
      for (discretizer in DISCRETIZERS){
        message(paste0("Dataset: ", dataset_name, " - Discretizer: ", discretizer))
        perform_feruc(dataset_name, discretizer)
      }
    }
  }
}
 
dataset_index <- strtoi(commandArgs(trailingOnly = TRUE)[1])
main(dataset_index)