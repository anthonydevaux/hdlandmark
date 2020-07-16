#' Title
#'
#' @param pred.surv
#' @param data.surv
#' @param tHor
#' @param simu
#'
#' @return
#' @export
#'
#' @importFrom pec pec
#' @importFrom timeROC timeROC
#' @importFrom survival Surv
#' @importFrom prodlim Hist
#'
#' @examples
LMassess <- function(pred.surv, data.surv, tHor){

  methods <- colnames(pred.surv)

  res <- list()

  res[["BS"]] <- res[["AUC"]] <- rep(NA, length = length(methods))
  names(res[["BS"]]) <- names(res[["AUC"]]) <- methods

  for (method in methods){

    # Brier score

    surv.method <- na.omit(cbind(rep(1,length(pred.surv[,method])), pred.surv[,method]))
    data.method <- data.surv[which(!is.na(pred.surv[,method])),, drop = FALSE]

    pec.method <- tryCatch(pec(object = surv.method, formula = Surv(time.event, event) ~ 1,
                               data = data.method, cens.model = "marginal",
                               exact = FALSE, times = c(0,tHor), ptime = 0),
                           error = function(e){return(NULL)})

    res[["BS"]][method] <- ifelse(!is.null(pec.method), pec.method$AppErr$matrix[2], NA)

    # AUC

    tHor.AUC <- tHor - 0.001

    ROC.method <- tryCatch(timeROC(T = data.method$time.event, delta = data.method$event,
                                   marker = 1 - surv.method[,2],
                                   cause = 1, iid = TRUE,
                                   times = tHor.AUC),
                           error = function(e){return(NULL)})

    res[["AUC"]][method] <- ifelse(!is.null(ROC.method), ROC.method$AUC[2], NA)

  }

  return(res)

}


