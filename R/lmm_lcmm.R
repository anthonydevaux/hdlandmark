#' Title
#'
#' @param model
#' @param formul
#' @param data
#' @param link
#' @param range
#' @param subject
#'
#' @return
#' @export
#'
#' @importFrom lcmm lcmm
#'
#' @examples
lmm.lcmm <- function(model, formul, data, link = "linear", range = NULL, subject){

  if (class(model)!="lcmm"){

    if (!is.null(model$range)){
      range <- model$range
    }

    if (!is.null(model$link)){
      link <- model$link
    }

    cat("lcmm modelling...")

    if (is.null(model$B)){
      model.output <- lcmm(fixed = model$fixed, random = model$random,
                           subject = subject, data = data,
                           link = link, range = range, verbose = FALSE)
    }else{
      model.output <- lcmm(fixed = model$fixed, random = model$random,
                           subject = subject, data = data,
                           link = link, range = range, verbose = FALSE,
                           B = model$B, posfix = model$posfix)
    }

    pred.RE <- predRE.lcmm(model = model.output, formul = model, data = data, subject = subject)

    return(list(model.output = model.output, pred.RE = pred.RE))

  }else{

    pred.RE <- predRE.lcmm(model = model, formul = formul, data = data, subject = subject)

    return(list(model.output = model, pred.RE = pred.RE))

  }

}
