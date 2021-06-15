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
#'
#' @examples
LMsurv.coxnet.CR <- function(data.surv, coxnet.submodels, cause = 1){

  model.coxnet <- list()

  best.lambda <- NULL
  best.alpha <- NULL

  data.surv.omit <- na.omit(data.surv[,!(names(data.surv) %in% "subject")])

  nb.cause <- sort(unique(data.surv.omit$event))[-1] # unique cause without censoring indicator

  lambda.min.lasso <- lambda.min.ridge <- list()
  lasso.coef <- ridge.coef <- c()

  ########################
  # tune lambda parameter for each cause for each penalty

  for (ind.cause in 1:length(nb.cause)){

    data.surv.cause <- data.surv.omit
    data.surv.cause$event <- ifelse(data.surv.cause$event==nb.cause[ind.cause],1,0) # censure for other cause

    data.surv.X <- model.matrix( ~ ., data.surv.cause[,!(names(data.surv.cause) %in% c("time.event","event"))])[,-1]
    data.surv.Y <- data.surv.cause[,c("time.event","event")]

    if (any(coxnet.submodels %in% c("lasso"))){

      nzero <- 0

      while (nzero==0){

        coxnet.fit.lasso <- glmnet::cv.glmnet(data.surv.X,
                                              Surv(data.surv.Y$time.event, data.surv.Y$event),
                                              family = "cox", alpha = 1,
                                              nfolds = 10)

        # no var issue
        best.cvm <- coxnet.fit.lasso$cvm[which(coxnet.fit.lasso$lambda==coxnet.fit.lasso$lambda.min)]
        nzero <- coxnet.fit.lasso$nzero[which(coxnet.fit.lasso$cvm==best.cvm)]

      }

      lambda.min.lasso[[ind.cause]] <- coxnet.fit.lasso$lambda.min

      coxnet.fit.lasso <- glmnet::glmnet(data.surv.X,
                                         Surv(data.surv.Y$time.event, data.surv.Y$event),
                                         family = "cox", lambda = lambda.min.lasso[[ind.cause]], alpha = 1)

      lasso.coef <- c(lasso.coef, as.numeric(coef(coxnet.fit.lasso)))

    }

    if (any(coxnet.submodels %in% c("ridge"))){

      coxnet.fit.ridge <- glmnet::cv.glmnet(data.surv.X,
                                            Surv(data.surv.Y$time.event, data.surv.Y$event),
                                            family = "cox", alpha = 0,
                                            nfolds = 10)

      lambda.min.ridge[[ind.cause]] <- coxnet.fit.ridge$lambda.min

      coxnet.fit.ridge <- glmnet::glmnet(data.surv.X,
                                         Surv(data.surv.Y$time.event, data.surv.Y$event),
                                         family = "cox", lambda = lambda.min.ridge[[ind.cause]], alpha = 0)

      ridge.coef <- c(ridge.coef, as.numeric(coef(coxnet.fit.ridge)))

    }

  }

  ######Cause Specific coxph fit #####

  time.event <- data.surv$time.event
  event <- factor(data.surv$event, 0:max(nb.cause))
  data.surv.X <- subset(data.surv, select = -c(time.event, event))

  allVar <- colnames(data.surv.X)[-which(colnames(data.surv.X)=="subject")]

  newformula <- reformulate(termlabels = allVar,
                            response = "Surv(time.event,event)")

  if (any(coxnet.submodels %in% c("lasso"))){

    coxph.lasso.fit <- survival::coxph(formula = newformula, data = data.surv.X, id = subject,
                                       init = lasso.coef, iter = 0, x = TRUE)


    model.coxnet[["lasso-CR"]] <- coxph.lasso.fit

  }

  if (any(coxnet.submodels %in% c("ridge"))){

    coxph.ridge.fit <- survival::coxph(formula = newformula, data = data.surv.X, id = subject,
                                       init = ridge.coef, iter = 0, x = TRUE)


    model.coxnet[["ridge-CR"]] <- coxph.ridge.fit

  }

  return(list(model = model.coxnet, lambda.opt = best.lambda, alpha.opt = best.alpha,
              surv.name = "penalized-cox"))

}
