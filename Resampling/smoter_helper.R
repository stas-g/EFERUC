#===============================================================================================
###HELPER FUNCTIONS for SMOTE/undersampler analysis
#===============================================================================================
#modification of SMOTER.R function from Torgo et al. to: 1) utomatically calculate all the necessary inputs for the phi function given method and exr.type which default to 'extremes' and 'both', resp. 2) handle errors when no extreme values are detected
my_smoter <- function(form, data,
                    thr.rel = 0.5,
                    perc.over = 200, k = 5,
                    perc.under = 200,
                    method = "extremes",
                    extr.type = "both",
                    control.pts = NULL,
                    learner = NULL,
                    tgt = 1, ... #tgt = location of the target variable
                    )
{

  y <- resp(form,data)
  pc <- list()
  if(method == "extremes") {
  phi_params <- phi.control(y, method = method, extr.type = extr.type)
  } else {
      phi_params <- phi.control(y, method = method, control.pts = control.pts)
}

  pc$method <- method
  pc$npts <- phi_params$npts
  pc$control.pts <- phi_params$control.pts

  both <- all(pc$control.pts[c(2,8)] == c(1,1))
  y.relev <- phi(y, pc)

  if (both) {  # we have both low and high extrs
    rare.low <- which(y.relev > thr.rel & y < pc$control.pts[4])
    rare.high <- which(y.relev > thr.rel & y > pc$control.pts[4])
    if(length(rare.low) == 0 & length(rare.high) == 0) return(NULL)
    if(length(rare.low) == 0) smote.exsL <- NULL else smote.exsL <- smote.exs(data[rare.low,],tgt,perc.over,k)
    if(length(rare.high) == 0) smote.exsH <- NULL else smote.exsH <- smote.exs(data[rare.high,],tgt,perc.over,k)
    rare.cases <- c(rare.low,rare.high)
    smoted.exs <- rbind(smote.exsL,smote.exsH)
  } else {

    # the indexes of the cases with rare target variable values
    rare.cases <- if (pc$control.pts[2] == 1)  which(y.relev > thr.rel & y < pc$control.pts[4]) else which(y.relev > thr.rel & y > pc$control.pts[4])
    if(length(rare.cases) == 0) return(NULL)
    # Get the smoted examples associated with these rare cases
    smoted.exs <- smote.exs(data[rare.cases,],tgt,perc.over,k)
  }

  # get the undersample of the "majority class" examples
  sel.maj <- sample((1:NROW(data))[-rare.cases],
#                    as.integer((perc.under/100)*nrow(smoted.exs)),
                    as.integer((perc.under/100)*(nrow(smoted.exs)+length(rare.cases))),
                    replace = TRUE)

  # the final data set (the undersample+the rare cases+the smoted exs)
  newdataset <- rbind(data[sel.maj,],data[rare.cases,],smoted.exs)

  # learn a model if required
  if (is.null(learner)) return(newdataset)
  else do.call(learner,list(form,newdataset,...))
}


#modification of random.underSample function from Torgo et al. to: 1) utomatically calculate all the necessary inputs for the phi function given method and exr.type which default to 'extremes' and 'both', resp. 2) handle errors when no extreme values are detected
my_undersampler <- function(form, data,
                               thr.rel = 0.5,
                               perc.under = 200,
                               method = "extremes",
                               extr.type = "both",
                               control.pts = NULL,
                               learner = NULL,
                               tgt = 1, ...) #tgt = location of the target variable
{

  y <- resp(form, data)
  pc <- list()
  if(method == "extremes") {
  phi_params <- phi.control(y, method = method, extr.type = extr.type)
  } else {
      phi_params <- phi.control(y, method = method, control.pts = control.pts)
}

  pc$method <- method
  pc$npts <- phi_params$npts
  pc$control.pts <- phi_params$control.pts

  both <- all(pc$control.pts[c(2,8)] == c(1,1))
  y.relev <- phi(y, pc)

  # the column where the target variable is
#  tgt <- target.col(form,data)

  # get the relevance values of the target variable values
#  x <- data[,tgt]
#  if (is.null(relev.pars)) x.relev <- phi(x)
#  else x.relev <- phi(x,phi.param=relev.pars)

  # the indexes of the cases with rare target variable values
  rare.cases <- which(y.relev > thr.rel)
  if(length(rare.cases) == 0) return(NULL)


  # start by adding all the cases of the minority class
  newdataset <- data[rare.cases,]

  # get the undersample of the "majority class" examples
  sel.maj <- sample((1:NROW(data))[-rare.cases],
                    as.integer((perc.under/100)*nrow(newdataset)),
                    replace=T)

  newdataset <- rbind(newdataset,data[sel.maj,])

  if (is.null(learner))
   {
     return(newdataset)
   } else {
      if(learner == "ranger") mod <- do.call(learner, list(form, newdataset,...))

      if(learner == "cv.glmnet") {
        x <- newdataset[, -tgt]
        y <- newdataset[, tgt]
        x <- xtn(x)
        x <- as.matrix(x)
        mod <- do.call(learner, list(x, y, ...))
  }
      if(learner == "xgboost") {
        x <- newdataset[, -tgt]
        y <- newdataset[, tgt]
        x <- xtn(x)
        dat <- xgb.DMatrix(as.matrix(x), label = y)
        mod <- do.call(learner, list(data = dat, ...))
        }

      return(mod)
    }
}


#===============================================================================================
### data manipulation functions
#===============================================================================================

#recode columns of x as factors; cols = vector of column names to be recoded. if col = NULL, then all columns are recoded (and are assumed to be non-character). if cols are provided, then only indicated columns are recoded; if a column is non-numeric, then levels of the resulting factor are coded as {1,...,number_unique_values}.
xtf <- function(x, cols = NULL) {
  if(is.matrix(x)) x <- as.data.frame(x)
  if(is.null(cols)) for(i in 1 : ncol(x)) x[, i] <- as.factor(x[, i]) else {
    for(a in cols) {
      z <- x[, a]
      if(is.numeric(z)) x[, a] <- as.factor(z) else {
        levs <- unique(z)
        x[, a] <- as.factor(match(z, unique(z)))
      }
    }
  }
  x
}

#QC for weeding out constant or constant variables out
xqc <- function(x) {
  ind <- diag(var(x))
  x[, ind > 0]
}

#convert factor SNP vars in x to numeric (i.e. columns are assumbed to be factors with levels {0, 1})
xtn <- function(x) {
  ii <- which(sapply(1 : ncol(x), FUN = function(i) class(x[, i]) != "numeric"))
  if(length(ii) > 0) for(k in ii) x[, k] <- (as.numeric(x[, k]) - 1)
  x
}


#calculate mse (rsq = FALSE) or r-squared (rsq = TRUE) for a vector (matrix, columnwise) of observed values and a vector (matrix, columnwise) of predicted values
calc.mse
function (obs, pred, rsq = FALSE)
{
    if (is.vector(obs))
        obs <- as.matrix(obs)
    if (is.vector(pred))
        pred <- as.matrix(pred)
    n <- nrow(obs)
    rss <- colSums((obs - pred)^2, na.rm = TRUE)
    if (rsq == FALSE)
        rss/n
    else {
        tss <- diag(var(obs, na.rm = TRUE)) * (n - 1)
        1 - rss/tss
    }
}
