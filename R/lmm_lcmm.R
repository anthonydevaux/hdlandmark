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
lmm.lcmm <- function(model, data){

  if (class(model)!="hlme"){

    model.output <- hlme(fixed = model$fixed, random = model$random, subject = model$subject,
                         data = data)
  }

  pred.RE <- predRE.lcmm(model.output, data)

  return(list(model.output = model.output, pred.RE = pred.RE))

}
