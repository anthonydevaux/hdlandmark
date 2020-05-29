#' Title
#'
#' @param model
#' @param data
#'
#' @return
#' @export
#'
#' @importFrom lme4 lmer
#'
#' @examples
lmm.lme4 <- function(model, data){

  if (class(model)!="lmerMod"){

    model.output <- lmer(model,
                         data = data, REML = FALSE,
                         control = lmerControl(optimizer ="Nelder_Mead"))
  }

  pred.RE <- predRE.lme4(model.output, data)

  return(list(model.output = model.output, pred.RE = pred.RE))

}
