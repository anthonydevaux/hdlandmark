#' Title
#'
#' @param data.surv
#' @param coxnet.opt
#' @param coxnet.lasso
#' @param coxnet.ridge
#'
#' @return
#' @export
#'
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom survival Surv coxph
#'
#' @examples
LMsurv.coxnet <- function(data.surv, coxnet.opt = TRUE, coxnet.lasso = coxnet.opt, coxnet.ridge = coxnet.opt){

  model.coxnet <- list()

  best.lambda <- NULL
  best.alpha <- NULL

  data.surv.omit <- na.omit(data.surv[,!(names(data.surv) %in% "subject")])
  data.surv.X <- model.matrix( ~ ., data.surv.omit[,!(names(data.surv.omit) %in% c("time.event","event"))])[,-1]
  data.surv.Y <- data.surv.omit[,c("time.event","event")]

  #######################################################

  # elastic net tuning

  if (coxnet.opt){

    best.cvm <- Inf

    for (current.alpha in seq(0,1,0.1)){

      coxnet.fit <- cv.glmnet(data.surv.X,
                                 Surv(data.surv.Y$time.event, data.surv.Y$event),
                                 family = "cox", alpha = current.alpha,
                                 nfolds = 10)

      current.cvm <- min(coxnet.fit$cvm)
      current.lambda <- coxnet.fit$lambda.min

      if (current.cvm < best.cvm){

        best.cvm <- current.cvm
        best.lambda <- current.lambda
        best.alpha <- current.alpha

      }

    }

    coxnet.fit <- glmnet(data.surv.X,
                          Surv(data.surv.Y$time.event, data.surv.Y$event),
                          family = "cox", lambda = best.lambda, alpha = best.alpha)

    coxnetFit.coef <- coef(coxnet.fit)
    coxnetFit.activVar <- coxnetFit.coef@Dimnames[[1]][coxnetFit.coef@i+1] # keep variables with non zero coefficients

    coxnetFit.coxph <- coxph(Surv(data.surv.Y$time.event, data.surv.Y$event) ~ .,
                               data = as.data.frame(data.surv.X[,coxnetFit.activVar]),
                               init = coxnetFit.coef@x, iter = 0, x = TRUE)

    model.coxnet[["opt"]] <- coxnetFit.coxph

  }

  browser()

  # lasso

  if (coxnet.lasso){

    coxnet.fit <- cv.glmnet(data.surv.X,
                            Surv(data.surv.Y$time.event, data.surv.Y$event),
                            family = "cox", alpha = 1,
                            nfolds = 10)

    coxnet.fit <- glmnet(data.surv.X,
                         Surv(data.surv.Y$time.event, data.surv.Y$event),
                         family = "cox", lambda = coxnet.fit$lambda.min, alpha = 1)

    coxnetFit.coef <- coef(coxnet.fit)
    coxnetFit.activVar <- coxnetFit.coef@Dimnames[[1]][coxnetFit.coef@i+1] # keep variables with non zero coefficients

    coxnetFit.coxph <- coxph(Surv(data.surv.Y$time.event, data.surv.Y$event) ~ .,
                             data = as.data.frame(data.surv.X[,coxnetFit.activVar]),
                             init = coxnetFit.coef@x, iter = 0, x = TRUE)

    model.coxnet[["lasso"]] <- coxnetFit.coxph

  }

  # ridge

  if (coxnet.ridge){

    coxnet.fit <- cv.glmnet(data.surv.X,
                            Surv(data.surv.Y$time.event, data.surv.Y$event),
                            family = "cox", alpha = 0,
                            nfolds = 10)

    coxnet.fit <- glmnet(data.surv.X,
                         Surv(data.surv.Y$time.event, data.surv.Y$event),
                         family = "cox", lambda = coxnet.fit$lambda.min, alpha = 0)

    coxnetFit.coef <- coef(coxnet.fit)
    coxnetFit.activVar <- coxnetFit.coef@Dimnames[[1]][coxnetFit.coef@i+1] # keep variables with non zero coefficients

    coxnetFit.coxph <- coxph(Surv(data.surv.Y$time.event, data.surv.Y$event) ~ .,
                             data = as.data.frame(data.surv.X[,coxnetFit.activVar]),
                             init = coxnetFit.coef@x, iter = 0, x = TRUE)

    model.coxnet[["ridge"]] <- coxnetFit.coxph

  }


  return(list(model = model.coxnet, lambda.opt = best.lambda, alpha.opt = best.alpha,
              surv.name = "penalized-cox"))

}
