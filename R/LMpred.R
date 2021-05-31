#' Title
#'
#' @param data.surv
#' @param model.surv
#' @param long.method
#' @param surv.methods
#' @param tHor
#' @param cause
#' @param CR
#'
#' @return
#' @export
#'
#' @importFrom survival survfit
#' @importFrom randomForestSRC predict.rfsrc
#' @import ranger
#' @import riskRegression
#'
#' @examples
LMpred <- function(data.surv, model.surv, long.method, surv.methods, tHor, cause = 1, CR = FALSE){

  cat("Survival prediction on data test","\n")
  cat("----------------------------------","\n")

  models <- unlist(lapply(model.surv, FUN = function(x) length(x$model)))
  models.nb <- sum(models)

  models.ind <- 1

  pred.surv <- matrix(NA, nrow = nrow(data.surv), ncol = models.nb,
                      dimnames = list(data.surv$subject,
                                      paste0("V", 1:models.nb)))

  for (surv.method in surv.methods){

    # Cox, Penalized-cox

    if (surv.method == "cox"){

      cat("Cox...")

      sub.methods <- names(model.surv[[surv.method]]$model)

      for (sub.method in sub.methods){

        model.current <- model.surv[[surv.method]]$model[[sub.method]]
        method.name <- paste(long.method, surv.method, sub.method, sep = "-")
        res.survfit <- tryCatch(survfit(model.current, data.surv), error = function(e){return(NULL)})
        id.time <- sum(res.survfit$time <= tHor)
        pred.surv[colnames(res.survfit$surv), models.ind] <- res.survfit$surv[id.time,]
        colnames(pred.surv)[models.ind] <- method.name
        models.ind <- models.ind + 1

      }

    }

    # penalized-cox

    if (surv.method == "penalized-cox"){

      cat("Penalized-Cox...")

      sub.methods <- names(model.surv[[surv.method]]$model)

      for (sub.method in sub.methods){

        model.current <- model.surv[[surv.method]]$model[[sub.method]]

        method.name <- paste(long.method, surv.method, sub.method, sep = "-")

        if (any(sub.method %in% c("opt","lasso","ridge"))){

          data.surv.coxnet <- as.data.frame(model.matrix( ~ ., na.omit(data.surv))[,-1])
          res.survfit <- tryCatch(survfit(model.current, data.surv.coxnet), error = function(e){return(NULL)})
          id.time <- sum(res.survfit$time <= tHor)
          pred.surv[colnames(res.survfit$surv), models.ind] <- res.survfit$surv[id.time,]
          colnames(pred.surv)[models.ind] <- method.name


        }

        if (any(sub.method %in% c("opt-CR","lasso-CR","ridge-CR"))){

          data.surv.coxnet <- as.data.frame(model.matrix( ~ ., na.omit(data.surv))[,-1])
          res.survfit <- riskRegression::predictCauseSpecificCox(model.current,
                                                                 newdata = data.surv.coxnet,
                                                                 times = tHor,
                                                                 cause = cause,
                                                                 product.limit = FALSE,
                                                                 #type = "survival")
                                                                 type = "absRisk")

          #pred.surv[rownames(na.omit(data.surv)), models.ind] <- res.survfit$survival[,1]
          pred.surv[rownames(na.omit(data.surv)), models.ind] <- res.survfit$absRisk[,1]
          colnames(pred.surv)[models.ind] <- method.name

        }

        models.ind <- models.ind + 1

      }
    }

    # Penalized-FG

    if (surv.method == "penalized-FG"){

      cat("penalized-FG...")

      sub.methods <- names(model.surv[[surv.method]]$model)

      for (sub.method in sub.methods){

        model.current <- model.surv[[surv.method]]$model[[sub.method]]
        method.name <- paste(long.method, surv.method, sub.method, sep = "-")

        data.surv.penaFG <- as.data.frame(model.matrix( ~ ., na.omit(data.surv))[,-1])

        pred.penaFG.fit <- predict(object = model.current,
                                   newdata = data.surv.penaFG,
                                   times = tHor)

        pred.surv[rownames(na.omit(data.surv)), models.ind] <- pred.penaFG.fit[,1]
        colnames(pred.surv)[models.ind] <- method.name
        models.ind <- models.ind + 1

      }

    }

    # sPLS

    if (surv.method == "spls"){

      cat("sPLS...")

      sub.methods <- names(model.surv[[surv.method]]$model)

      for (sub.method in sub.methods){

        model.current <- model.surv[[surv.method]]$model[[sub.method]]

        method.name <- paste(long.method, surv.method, sub.method, sep = "-")

        Xnames <- rownames(model.current$splsDR_modplsr$loadings$X)

        data.surv.spls <- as.data.frame(model.matrix( ~ ., na.omit(data.surv))[,-1])

        # centre reduit la matrice des nouveaux individus à partir du mean/sd du train
        Xh.scale <- t((t(data.surv.spls[,Xnames])-model.current$XplanCent[Xnames])/model.current$XplanScal[Xnames])

        X.spls <- matrix(NA, nrow = nrow(Xh.scale), ncol = ncol(model.current$tt_splsDR),
                         dimnames = list(rownames(Xh.scale), colnames(model.current$tt_splsDR)))

        u <- model.current$splsDR_modplsr$loadings$X

        X.spls[,1] <- Xh.scale%*%u[,1]

        if (ncol(X.spls) > 1){

          for (h in 2:ncol(X.spls)){

            th <- Xh.scale%*%u[,h-1]

            proj.num <- th%*%t(th)
            proj.den <- as.numeric(t(th)%*%th)
            proj <- proj.num / proj.den

            Xh.scale <- Xh.scale - proj%*%Xh.scale

            Xh <- Xh.scale%*%u[,h]

            X.spls[,h] <- Xh

          }

        }

        X.spls.df <- as.data.frame(X.spls)
        rownames(X.spls.df) <- rownames(data.surv.spls)

        res.survfit <- survfit(model.current$cox_splsDR, X.spls.df)

        id.time <- sum(res.survfit$time <= tHor)

        pred.surv[colnames(res.survfit$surv), models.ind] <- res.survfit$surv[id.time,]

        colnames(pred.surv)[models.ind] <- method.name

        models.ind <- models.ind + 1

      }

    }

    # rsf

    if (surv.method == "rsf"){

      cat("RSF...")

      sub.methods <- names(model.surv[[surv.method]]$model)

      for (sub.method in sub.methods){

        model.current <- model.surv[[surv.method]]$model[[sub.method]]
        method.name <- paste(long.method, surv.method, sub.method, sep = "-")

        if (any(sub.method %in% c("logrank-opt","logrank-noVS","logrank-default",
                                  "bs.gradient-opt","bs.gradient-noVS","bs.gradient-default"))){

          res.survfit <- predict.rfsrc(model.current, data.surv)
          id.time <- sum(res.survfit$time.interest <= tHor)
          formula.xvar <- as.formula(as.character(model.current$call$formula)[c(1,3)])
          id.noNA <- rownames(model.frame(formula.xvar,
                                          data.surv[,model.current$xvar.names, drop = FALSE])) # id without NA
          pred.surv[id.noNA, models.ind] <- res.survfit$survival[,id.time]

        }

        if (any(sub.method %in% c("logrank-opt-CR","logrank-noVS-CR","logrank-default-CR"))){

          res.survfit <- predict.rfsrc(model.current, data.surv)
          id.time <- sum(res.survfit$time.interest <= tHor)
          formula.xvar <- as.formula(as.character(model.current$call$formula)[c(1,3)])
          id.noNA <- rownames(model.frame(formula.xvar,
                                          data.surv[,model.current$xvar.names, drop = FALSE])) # id without NA
          pred.surv[id.noNA, models.ind] <- res.survfit$cif[,id.time, cause]

        }

        if (any(sub.method %in% c("ranger"))){

          data.surv.omit <- na.omit(data.surv)
          id.noNA <- id.noNA <- rownames(data.surv.omit)
          res.survfit <- predict(model.current, data.surv.omit)
          id.time <- sum(res.survfit$unique.death.times <= tHor)
          pred.surv[id.noNA, models.ind] <- res.survfit$survival[,id.time]

        }

        colnames(pred.surv)[models.ind] <- method.name
        models.ind <- models.ind + 1

      }

    }

  }

  cat("--", "\n")

  return(list(pred.surv = pred.surv, models.name = names(models)))

}
