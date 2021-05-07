#' Title
#'
#' @param model
#' @param formul
#' @param data
#'
#' @return
#' @export
#'
#' @importFrom lcmm lcmm
#'
#' @examples
lmm.lcmm <- function(model, formul, data){

  if (class(model)!="lcmm"){

    model.output <- lcmm(fixed = model$fixed, random = model$random,
                         subject = model$subject, data = data,
                         link = ifelse(!is.null(model$link), model$link, "linear"))

    pred.RE <- predRE.lcmm(model = model.output, formul = model, data = data)

    return(list(model.output = model.output, pred.RE = pred.RE))

  }else{

    pred.RE <- predRE.lcmm(model = model, formul = formul, data = data)

    return(list(model.output = model, pred.RE = pred.RE))

  }

}
