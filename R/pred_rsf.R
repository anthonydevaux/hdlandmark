#' Title
#'
#' @param model
#' @param method
#' @param newdata
#' @param var_list
#' @param tHor
#' @param pred_surv
#'
#' @return
#' @export
#'
#' @examples
pred.rsf <- function(model, method, newdata, var_list, tHor, pred_surv){

  n <- nrow(newdata)

  if (!is.null(var_list)){

    newdata <- newdata[,var_list]

  }

  if (n > 1){

    res_survfit <- randomForestSRC::predict.rfsrc(model, newdata)
    id_time <- sum(res_survfit$time.interest <= tHor)
    formula.xvar <- as.formula(as.character(model$call$formula)[c(1,3)])
    id_no_na <- rownames(model.frame(formula.xvar,
                                     newdata[,model$xvar.names, drop = FALSE])) # id without NA
    pred_surv[id_no_na, method] <- res_survfit$survival[,id_time]

  }else{

    if (any(is.na(newdata[,model$xvar.names]))){

      pred_surv[1, method] <- NA

    }else{

      res_survfit <- randomForestSRC::predict.rfsrc(model, newdata)
      id_time <- sum(res_survfit$time.interest <= tHor)
      pred_surv[1, method] <- res_survfit$survival[id_time]

    }

  }

  return(pred_surv)

}
