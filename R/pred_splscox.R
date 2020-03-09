#' Title
#'
#' @param model
#' @param method
#' @param newdata
#' @param tHor
#' @param pred_surv
#'
#' @return
#' @export
#'
#' @examples
pred.splscox <- function(model, method, newdata, tHor, pred_surv){

  Xnames <- rownames(model$splsDR_modplsr$loadings$X)

  # centre reduit la matrice des nouveaux individus à partir du mean/sd du train
  Xh.scale <- t((t(newdata[,Xnames])-model$XplanCent[Xnames])/model$XplanScal[Xnames])

  X.spls <- matrix(NA, nrow = nrow(Xh.scale), ncol = ncol(model$tt_splsDR),
                   dimnames = list(rownames(Xh.scale), colnames(model$tt_splsDR)))

  u <- model$splsDR_modplsr$loadings$X

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

  pred_surv <- pred.phm(model = model$cox_splsDR, method = method,
                        newdata = as.data.frame(X.spls), var_list = NULL,
                        tHor = tHor, pred_surv = pred_surv)

  return(pred_surv)

}
