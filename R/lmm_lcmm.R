#' Title
#'
#' @param model
#' @param formul
#' @param data
#' @param link
#' @param range
#'
#' @return
#' @export
#'
#' @importFrom lcmm lcmm
#'
#' @examples
lmm.lcmm <- function(model, formul, data, link = "linear", range = NULL){

  if (class(model)!="lcmm"){

    if (!is.null(model$range)){
      range <- model$range
    }

    if (!is.null(model$link)){
      link <- model$link
    }

    cat("lcmm modelling...")

    model.output <- lcmm(fixed = model$fixed, random = model$random,
                         subject = model$subject, data = data,
                         link = link, range = range, verbose = FALSE)

    pred.RE <- predRE.lcmm(model = model.output, formul = model, data = data)

    return(list(model.output = model.output, pred.RE = pred.RE))

  }else{

    pred.RE <- predRE.lcmm(model = model, formul = formul, data = data)

    return(list(model.output = model, pred.RE = pred.RE))

  }

}
