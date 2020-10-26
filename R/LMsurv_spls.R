#' Title
#'
#' @param data.surv
#' @param spls.submodels
#'
#' @return
#' @export
#'
#' @importFrom plsRcox cv.coxsplsDR coxsplsDR
#'
#' @examples
LMsurv.spls <- function(data.surv, spls.submodels){

  model.spls <- list()

  data.surv.omit <- na.omit(data.surv[,!(names(data.surv) %in% "subject")])
  data.surv.X <- model.matrix( ~ ., data.surv.omit[,!(names(data.surv.omit) %in% c("time.event","event"))])[,-1]
  data.surv.Y <- data.surv.omit[,c("time.event","event")]

  best.eta <- best.ncomp <- NULL

  ##############################################
  # drop variables with null standard deviation

  var.nullsd <- which(apply(data.surv.X, 2, sd)==0)

  if (length(var.nullsd)>0){

    data.surv.X <- data.surv.X[,-var.nullsd]

  }

  ##############################################

  if (any(spls.submodels %in% c("opt"))){

    best.auc <- 0.5

    for (eta in seq(0,0.9,0.1)){ # loop for eta tuning (Chun et Keles, 2010)

      cat(eta,"\n")

      error.flag <- "error"
      error.ind <- 0

      while(error.flag=="error"&error.ind<10){ # loop for fold issues

        error.flag <- tryCatch(cv.coxsplsDR(list(x = data.surv.X, time = data.surv.Y$time.event,
                                                 status = data.surv.Y$event),
                                            eta = eta, nt = 10, nfold = 5, plot.it = FALSE,
                                            allCVcrit = FALSE, details = TRUE),
                               error = function(e){return("error")})

        error.ind <- error.ind + 1

      }

      cv.splsdrFit <- error.flag

      temp.auc <- cv.splsdrFit$cv.error10[cv.splsdrFit$lambda.min10+1]

      if (temp.auc > best.auc){
        best.auc <- temp.auc
        best.eta <- eta
        best.ncomp <- cv.splsdrFit$lambda.min10
      }
    }

    res.spls <- coxsplsDR(Xplan = data.surv.X, time = data.surv.Y$time.event,
                          time2 = data.surv.Y$event,
                          ncomp = best.ncomp, eta = best.eta,
                          trace = TRUE, allres = TRUE)

    model.spls[["opt"]] <- res.spls

  }

  if (any(spls.submodels %in% c("nosparse"))){

    error.flag <- "error"
    error.ind <- 0

    while(error.flag=="error"&error.ind<10){

      error.flag <- tryCatch(cv.coxsplsDR(list(x = data.surv.X, time = data.surv.Y$time.event,
                                               status = data.surv.Y$event),
                                          eta = 0, nt = 10, nfold = 5, plot.it = FALSE,
                                          allCVcrit = FALSE, details = TRUE),
                             error = function(e){return("error")})

      error.ind <- error.ind + 1

    }

    cv.splsdrFit <- error.flag

    res.spls <- coxsplsDR(Xplan = data.surv.X, time = data.surv.Y$time.event,
                          time2 = data.surv.Y$event,
                          ncomp = cv.splsdrFit$lambda.min10, eta = 0,
                          trace = TRUE, allres = TRUE)

    model.spls[["nosparse"]] <- res.spls

  }

  if (any(spls.submodels %in% c("maxsparse"))){

    error.flag <- "error"
    error.ind <- 0

    while(error.flag=="error"&error.ind<10){

      error.flag <- tryCatch(cv.coxsplsDR(list(x = data.surv.X, time = data.surv.Y$time.event,
                                               status = data.surv.Y$event),
                                          eta = 0.9, nt = 10, nfold = 5, plot.it = FALSE,
                                          allCVcrit = FALSE, details = TRUE),
                             error = function(e){return("error")})

      error.ind <- error.ind + 1

    }

    cv.splsdrFit <- error.flag

    res.spls <- coxsplsDR(Xplan = data.surv.X, time = data.surv.Y$time.event,
                          time2 = data.surv.Y$event,
                          ncomp = cv.splsdrFit$lambda.min10, eta = 0.9,
                          trace = TRUE, allres = TRUE)

    model.spls[["maxsparse"]] <- res.spls

  }

  return(list(model = model.spls, eta.opt = best.eta, ncomp.opt = best.ncomp,
              surv.name = "spls"))
}
