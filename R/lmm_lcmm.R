#' Title
#'
#' @param model
#' @param data
#'
#' @return
#' @export
#'
#' @importFrom lcmm hlme
#'
#' @examples
lmm.lcmm <- function(model, formul, data){

  if (class(model)!="hlme"){

    model.output <- hlme(fixed = model$fixed, random = model$random, subject = model$subject,
                  data = data)

    pred.RE <- predRE.lcmm(model = model.output, formul = model, data = data)

    return(list(model.output = model.output, pred.RE = pred.RE))

  }else{

    pred.RE <- predRE.lcmm(model = model, formul = formul, data = data)

    return(list(model.output = model, pred.RE = pred.RE))

  }

}
