library(caret)
library(ranger)
library(glmnet)
library(xgboost)
library(funfun)
library(randomFunctions) #devtools::install_github('chr1swallace/random-functions') or remotes::install_github('chr1swallace/random-functions')

## extremesResampling.R can be obtained from the Code and Data section of https://www.dcc.fc.up.pt/~ltorgo/ExpertSystems/ (Torgo et al.  Re-sampling Strategies for Regression. Expert Systems. DOI: 10.1111/exsy.1208)
source("CodeAndData/extremesResampling.R")
source("smoter_helper.R")

# choose parameters: dataset collection (dat), dataset (a), learner, resampling technique (r), undersampling level (u), and oversampling level (o).
ARGS <- getArgs()
dat <- ARGS$d
a <- ARGS$a
smoter <- as.logical(ARGS$s)
perc.under <- as.numeric(ARGS$u)

if(smoter)  perc.over <- as.numeric(ARGS$o)

# information about whether datset contains nominal variables and, if so, which ones
NOM <- readRDS("nominal_var_info.rds")[[dat]]
# bespoke relevance parameters for datasets that require them
REL <- readRDS("relevance.rds")[[dat]]

#===============================================================================================
PATH <- paste0("datasets_and_splits/datasets/", dat)
if(dat != "PaoBrancoImbalanced") {
  INFO <- readLines(file.path(PATH, "datasets.txt"), warn = FALSE) #all datasets
} else {
  INFO <- gsub("\\.csv", "", list.files(PATH, pattern = "csv$"))
}
INFO_NN <- INFO[!(INFO %in% names(NOM))] #non-nominal datasets

if(!file.exists('output')) dir.create('output')
#===============================================================================================
## reading data in

if(dat == "gene_expression") {
  x_train <- read.csv(file.path(PATH, "X_train.csv"), row.names = 1)
  x_test <- read.csv(file.path(PATH, "X_test.csv"), row.names = 1)

  ind1 <- colnames(xqc(x_train))
  ind2 <- colnames(xqc(x_test))
  ind <- intersect(ind1, ind2)
  x_train <- xtf(x_train[, ind])
  x_test <- xtf(x_test[, ind])

  y_train <- read.csv(file.path(PATH, paste0(a, "_Ys_train.csv")), row.names = 1)
  y_test <- read.csv(file.path(PATH, paste0(a, "_Ys_test.csv")), row.names = 1)[,1]

  data_train <- cbind(y_train, x_train)
  data_test <- cbind(y_test, x_test)

  method <- "extremes"
  control.pts <- NULL
}


if(dat %in% c("OpenML", "QSAR", "Yeast")) {
  x <- read.csv(sprintf("%s/%s_X.csv", PATH, a), row.names = 1)
  y <- read.csv(sprintf("%s/%s_Ys.csv", PATH, a), row.names = 1)

  if(dat != "OpenML") {
    x <- xtf(x)
  } else {
    if(a %in% names(NOM)) x <- xtf(x, cols = NOM[[a]])
  }

  test_ind <- read.csv(sprintf("datasets_and_splits/splits/%s/%s_test.txt", dat, a))[,1]
  train_ind <- read.csv(sprintf("datasets_and_splits/splits/%s/%s_train.txt", dat, a))[,1]

  x_train <- x[rownames(x) %in% train_ind, ]
  x_test <- x[rownames(x) %in% test_ind, ]

  y_test <- y[rownames(y) %in% test_ind, 1]

  ind1 <- colnames(xqc(x_train))
  ind2 <- colnames(xqc(x_test))
  ind <- intersect(ind1, ind2)
  data_train <- cbind(y[rownames(y) %in% train_ind,,drop = FALSE], x_train[, ind])
  data_test <- cbind(y[rownames(y) %in% test_ind,,drop = FALSE], x_test[, ind])

  if(a %in% names(REL)) {
    method <- "range"
    control.pts <- REL[[a]]
  } else {
      method <- "extremes"
      control.pts <- NULL
  }
}


if(dat == "PaoBrancoImbalanced") {
  data <- read.csv(sprintf("%s/%s.csv", PATH, a))

  if(a %in% names(NOM)) data <- xtf(data, cols = NOM[[a]])

  train_ind <- read.csv(sprintf("datasets_and_splits/splits/%s/%s_train.txt", dat, a))[,1]
  test_ind <- read.csv(sprintf("datasets_and_splits/splits/%s/%s_test.txt", dat, a))[,1]

  data_train <- data[train_ind, ]
  data_test <- data[test_ind, ]

  y_test <- data[test_ind, 1]

  method <- "extremes"
  control.pts <- NULL
}

#===============================================================================================
## create a SmotED dataset
o <- ifelse(smoter, sprintf("o%s_", perc.over), "")

form <- as.formula(paste(colnames(data_train)[1], "~ ."))
if(smoter) {
  tryCatch(newdataset <- my_smoter(form, data = data_train,
        thr.rel = 0.7,
        perc.over = perc.over,
        perc.under = perc.under,
        k = 5,
        method = method,
        control.pts = control.pts,
        learner = NULL), error = function(e) err <<- e)
 } else {
      tryCatch(newdataset <- my_undersampler(form = form, data = data_train,
          thr.rel = 0.7,
          perc.under = perc.under,
          method = method,
          control.pts = control.pts,
          learner = NULL), error = function(e) err <<- e)
}

if(exists('err')) stop("smoter failed! halting...")
if(is.null(newdataset)) stop("no extreme values in dataset! halting...")


##run RF
rfmod <- ranger(form, data = newdataset,
    num.trees = 1000,
    mtry = floor((ncol(newdataset) - 1)/3))

rfpred <- predict(rfmod, data_test)$predictions
rfrsq <- calc.mse(y_test, rfpred, rsq = TRUE)

saveRDS(rfmod, file = sprintf('output/%s_%s_ranger_%s_u%s_%smod.rds', dat, a, ifelse(smoter, "smoter", "undersampler"), perc.under, o))
saveRDS(rfpred, file = sprintf('output/%s_%s_ranger_%s_u%s_%spred.rds', dat, a, ifelse(smoter, "smoter", "undersampler"), perc.under, o))
saveRDS(rfrsq, file = sprintf('output/%s_%s_ranger_%s_u%s_%srsq.rds', dat, a, ifelse(smoter, "smoter", "undersampler"), perc.under, o))


## if dataset a doesn't have any nominal variables, run lasso, ridge and xgb
if(a %in% INFO_NN) {
    ##LASSO
    x <- newdataset[, -1]
    y <- newdataset[, 1]
    x <- xtn(x)
    x <- as.matrix(x)

    lmod <- cv.glmnet(x = x, y = y,
        alpha = 1,
        family = "gaussian",
        standardize = FALSE,
        lambda = 10^seq(10, -3, length=100)
      )

    xnew <- as.matrix(xtn(data_test[, -1]))
    lpred  <- predict(lmod, xnew, s = "lambda.1se")
    lrsq <- calc.mse(y_test, lpred, rsq = TRUE)

    saveRDS(lmod, file = sprintf('output/%s_%s_lasso_%s_u%s_o%s_mod.rds', dat, a, ifelse(smoter, "smoter", "undersampler"), perc.under, o))
    saveRDS(lpred, file = sprintf('output/%s_%s_lasso_%s_u%s_o%s_pred.rds', dat, a, ifelse(smoter, "smoter", "undersampler"), perc.under, o))
    saveRDS(lrsq, file = sprintf('output/%s_%s_lasso_%s_u%s_o%s_rsq.rds', dat, a, ifelse(smoter, "smoter", "undersampler"), perc.under, o))

    ##RIDGE
    x <- newdataset[, -1]
    y <- newdataset[, 1]
    x <- xtn(x)
    x <- as.matrix(x)

    rmod <- cv.glmnet(x = x, y = y,
        alpha = 0,
        family = "gaussian",
        standardize = FALSE,
        lambda = 10^seq(10, -3, length=100)
      )

    xnew <- as.matrix(xtn(data_test[, -1]))
    rpred  <- predict(rmod, xnew, s = "lambda.1se")
    rrsq <- calc.mse(y_test, rpred, rsq = TRUE)

    saveRDS(rmod, file = sprintf('output/%s_%s_ridge_%s_u%s_o%s_mod.rds', dat, a, ifelse(smoter, "smoter", "undersampler"), perc.under, o))
    saveRDS(rpred, file = sprintf('output/%s_%s_ridge_%s_u%s_o%s_pred.rds', dat, a, ifelse(smoter, "smoter", "undersampler"), perc.under, o))
    saveRDS(rrsq, file = sprintf('output/%s_%s_ridge_%s_u%s_o%s_rsq.rds', dat, a, ifelse(smoter, "smoter", "undersampler"), perc.under, o))

    ##XGBOOST
    grid <- expand.grid(
        nrounds = 1500,
        max_depth = 6,
        eta = c(0.001, 0.01, 0.1, 0.2, 0.3),
        gamma = 0.05,
        subsample = 0.5,
        colsample_bytree = 1,
        min_child_weight = 1
      )

      ctrl <- trainControl(
        method = "LGOCV",
        p = 0.7,
        number = 1,
        allowParallel = FALSE
      )

      xgbmod <- train(
        x = xtn(newdataset[, -1]),
        y = newdataset[, 1],
        method = "xgbTree",
        eval_metric = "rmse",
        verbose = 1,
        trControl = ctrl,
        tuneGrid = grid,
        objective = "reg:squarederror",
      )

      xnew <- as.matrix(xtn(data_test[, -1]))
      xgbpred <- predict(xgbmod$finalModel, xnew)
      xgbrsq <- calc.mse(y_test, xgbpred, rsq = TRUE)

      saveRDS(xgbmod, file = sprintf('output/%s_%s_xgboost_%s_u%s_o%s_mod.rds', dat, a, ifelse(smoter, "smoter", "undersampler"), perc.under, o))
      saveRDS(xgbpred, file = sprintf('output/%s_%s_xgboost_%s_u%s_o%s_pred.rds', dat, a, ifelse(smoter, "smoter", "undersampler"), perc.under, o))
      saveRDS(xgbrsq, file = sprintf('output/%s_%s_xgboost_%s_u%s_o%s_rsq.rds', dat, a, ifelse(smoter, "smoter", "undersampler"), perc.under, o))
}
