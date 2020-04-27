#' Title
#'
#' @param data.surv
#' @param model.surv
#' @param long.methods
#' @param surv.methods
#' @param tHor
#'
#' @return
#' @export
#'
#' @importFrom survival survfit
#' @importFrom randomForestSRC predict.rfsrc
#'
#' @examples
LMpred <- function(data.surv, model.surv, long.methods, surv.methods, tHor){

  models <- unlist(lapply(model.surv, FUN = function(x) {lapply(x, FUN = function(y) length(y$model))}))
  models.nb <- sum(models)

  pred.surv <- matrix(NA, nrow = nrow(data.surv[[1]]), ncol = models.nb,
                      dimnames = list(data.surv[[1]]$subject,
                                      paste0("V", 1:models.nb)))

  models.ind <- 1

  for (long.method in long.methods){

    for (surv.method in surv.methods){

      # Cox, Penalized-cox

      if (any(surv.method %in% c("cox", "penalized-cox"))){

        sub.methods <- names(model.surv[[long.method]][[surv.method]]$model)

        for (sub.method in sub.methods){

          model.current <- model.surv[[long.method]][[surv.method]]$model[[sub.method]]

          method.name <- paste(long.method, surv.method, sub.method, sep = "-")

          if (surv.method == "cox"){

            res.survfit <- survfit(model.current, data.surv[[long.method]])

          }

          if (surv.method == "penalized-cox"){

            data.surv.coxnet <- as.data.frame(model.matrix( ~ ., na.omit(data.surv[[long.method]]))[,-1])
            res.survfit <- survfit(model.current, data.surv.coxnet)

          }

          id.time <- sum(res.survfit$time <= tHor)

          pred.surv[colnames(res.survfit$surv), models.ind] <- res.survfit$surv[id.time,]

          colnames(pred.surv)[models.ind] <- method.name

          models.ind <- models.ind + 1

        }

      }

      # sPLS

      if (surv.method == "sPLS"){

        sub.methods <- names(model.surv[[long.method]][[surv.method]]$model)

        for (sub.method in sub.methods){

          model.current <- model.surv[[long.method]][[surv.method]]$model[[sub.method]]

          method.name <- paste(long.method, surv.method, sub.method, sep = "-")

          Xnames <- rownames(model.current$splsDR_modplsr$loadings$X)

          data.surv.spls <- as.data.frame(model.matrix( ~ ., na.omit(data.surv[[long.method]]))[,-1])

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

        sub.methods <- names(model.surv[[long.method]][[surv.method]]$model)

        for (sub.method in sub.methods){

          browser()

          model.current <- model.surv[[long.method]][[surv.method]]$model[[sub.method]]
          method.name <- paste(long.method, surv.method, sub.method, sep = "-")

          res.survfit <- predict.rfsrc(model.current, data.surv[[long.method]])
          id.time <- sum(res.survfit$time.interest <= tHor)
          formula.xvar <- as.formula(as.character(model.current$call$formula)[c(1,3)])
          id.noNA <- rownames(model.frame(formula.xvar,
                                           data.surv[[long.method]][,model.current$xvar.names, drop = FALSE])) # id without NA
          pred.surv[id.noNA, models.ind] <- res.survfit$survival[,id.time]

          colnames(pred.surv)[models.ind] <- method.name

          models.ind <- models.ind + 1

        }

      }

    }

  }

}
