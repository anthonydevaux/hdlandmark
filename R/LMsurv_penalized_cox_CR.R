#' Title
#'
#' @param data.surv
#' @param coxnet.submodels
#' @param cause
#'
#' @return
#' @export
#'
#' @import glmnet
#' @importFrom survival Surv coxph
#' @import riskRegression
#'
#' @examples
LMsurv.coxnet.CR <- function(data.surv, coxnet.submodels, cause = 1){

  model.coxnet <- list()

  best.lambda <- NULL
  best.alpha <- NULL

  data.surv.omit <- na.omit(data.surv[,!(names(data.surv) %in% "subject")])

  nb.cause <- sort(unique(data.surv.omit$event))[-1] # unique cause without censoring indicator
  nb.cause <- c(nb.cause[cause],nb.cause[-cause]) # order with interest cause in first

  lambda.min.lasso <- lambda.min.ridge <- lambda.min.EN <- list()

  ########################
  # tune lambda parameter for each cause for each penalty

  for (ind.cause in 1:length(nb.cause)){

    data.surv.cause <- data.surv.omit
    data.surv.cause$event <- ifelse(data.surv.cause$event==nb.cause[ind.cause],1,0) # censure for other cause

    data.surv.X <- model.matrix( ~ ., data.surv.cause[,!(names(data.surv.cause) %in% c("time.event","event"))])[,-1]
    data.surv.Y <- data.surv.cause[,c("time.event","event")]

    if (any(coxnet.submodels %in% c("lasso"))){

      coxnet.fit.lasso <- glmnet::cv.glmnet(data.surv.X,
                                            Surv(data.surv.Y$time.event, data.surv.Y$event),
                                            family = "cox", alpha = 1,
                                            nfolds = 10)

      lambda.min.lasso[[ind.cause]] <- coxnet.fit.lasso$lambda.min

    }

    if (any(coxnet.submodels %in% c("ridge"))){

      coxnet.fit.ridge <- glmnet::cv.glmnet(data.surv.X,
                                            Surv(data.surv.Y$time.event, data.surv.Y$event),
                                            family = "cox", alpha = 0,
                                            nfolds = 10)

      lambda.min.ridge[[ind.cause]] <- coxnet.fit.ridge$lambda.min

    }

  }

  ###### CSC fit #####

  allVar <- colnames(data.surv.X)

  newformula <- reformulate(termlabels = allVar,
                            response = "Hist(time.event,event)")

  if (any(coxnet.submodels %in% c("lasso"))){

    CSC.lasso.fit <- riskRegression::CSC(newformula,
                                         data = data.surv.omit,
                                         cause = cause, fitter = "glmnet",
                                         lambda = lambda.min.lasso,
                                         alpha = 1)


    model.coxnet[["lasso-CR"]] <- CSC.lasso.fit

  }

  if (any(coxnet.submodels %in% c("ridge"))){

    CSC.ridge.fit <- riskRegression::CSC(newformula,
                                         data = data.surv.omit,
                                         cause = cause, fitter = "glmnet",
                                         lambda = lambda.min.ridge,
                                         alpha = 0)


    model.coxnet[["ridge-CR"]] <- CSC.ridge.fit

  }

  return(list(model = model.coxnet, lambda.opt = best.lambda, alpha.opt = best.alpha,
              surv.name = "penalized-cox"))

}
