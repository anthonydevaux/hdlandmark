#' Title
#'
#' @param data.surv
#' @param cox.autoVar
#' @param cox.allVar
#'
#' @return
#' @export
#'
#' @importFrom survival Surv coxph
#'
#' @examples
LMsurv.cox <- function(data.surv, cox.autoVar = TRUE, cox.allVar = !cox.autoVar){

  model.cox <- list()

  data.surv <- data.surv[,!(names(data.surv) %in% "subject")]

  #### Cox with variables selection

  if (cox.autoVar){

    data.surv.omit <- na.omit(data.surv) # omit values to compute step function

    best.coxFit <- step(coxph(Surv(time.event, event) ~ ., data = data.surv.omit))

    coxFit <- coxph(best.coxFit$formula, data = data.surv, x = TRUE)

    model.cox[["autoVar"]] <- coxFit

  }

  #### Cox without variables selection (all variables)

  if (cox.allVar){

    coxFit <- coxph(Surv(time.event, event) ~ ., data = data.surv, x = TRUE)

    model.cox[["allVar"]] <- coxFit

  }

  return(list(model = model.cox, surv.name = "cox"))

}
