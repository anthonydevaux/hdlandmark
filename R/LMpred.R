#' Title
#'
#' @param method
#' @param cox_object
#' @param data
#' @param time
#' @param event
#' @param tHor
#'
#' @return
#' @export
#'
#' @examples
#'
LMpred <- function(method = c("cox","rf"), cox_object, data, time, event, tHor){

  if (method == "cox"){
    if (is.null(cox_object)){
      stop("cox_object is needed for method = cox", "\n")
    }

    res_survfit <- survfit(cox_object, data)
    id_time <- sum(res_survfit$time <= tHor)

    res_surv <- data.frame("S_tHor" = res_survfit$surv[id_time,],
                           "ICinf" = res_survfit$lower[id_time,],
                           "ICsup" = res_survfit$upper[id_time,])

  }

  if (method == "rf"){

  }

}
