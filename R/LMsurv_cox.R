#' Title
#'
#' @param data.surv
#' @param cox.submodels
#'
#' @return
#' @export
#'
#' @importFrom survival Surv coxph
#'
#' @examples
LMsurv.cox <- function(data.surv, cox.submodels){

  model.cox <- list()

  data.surv <- data.surv[,!(names(data.surv) %in% "subject")]

  #### Cox with variables selection

  if (any(cox.submodels %in% c("autoVar"))){

    data.surv.omit <- na.omit(data.surv) # omit values to compute step function

    best.coxFit <- tryCatch(step(coxph(Surv(time.event, event) ~ ., data = data.surv.omit), trace = 0),
                            error = function(e) {return(NULL)})

    if (!is.null(best.coxFit)){

      coxFit <- coxph(best.coxFit$formula, data = data.surv, x = TRUE)

    }else{

      coxFit <- "error"

    }

    model.cox[["autoVar"]] <- coxFit

  }

  #### Cox without variables selection (all variables)

  if (any(cox.submodels %in% c("allVar"))){

    coxFit <- coxph(Surv(time.event, event) ~ ., data = data.surv, x = TRUE)

    model.cox[["allVar"]] <- coxFit

  }

  return(list(model = model.cox, surv.name = "cox"))

}
