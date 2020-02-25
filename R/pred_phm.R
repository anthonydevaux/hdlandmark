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
pred.phm <- function(model, method, newdata, var_list, tHor, pred_surv){

  n <- nrow(newdata)

  if (!is.null(var_list)){

    newdata <- newdata[,var_list]

  }

  res_survfit <- survival::survfit(model, newdata)
  id_time <- sum(res_survfit$time <= tHor)

  if (n > 1){

    pred_surv[colnames(res_survfit$surv), method] <- res_survfit$surv[id_time,]

  }else{

    pred_surv[1, method] <- res_survfit$surv[id_time]

  }

  return(pred_surv)

}
