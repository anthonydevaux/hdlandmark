#' Function to predict the value of longitudinal outcome for a specific time
#'
#' @param predRE A list object from \code{predRE} function
#' @param data A dataframe where each row containing some predictive variables for a specific subject
#'
#' @return A matrix containing the prediction value of the longitudinal outcome for each subject
#' @export
#'
#' @examples
#'
predY_ind <- function(predRE, data, time, tLM, threshold = NULL){

  if (is.null(threshold)){

    pred_Y <- unlist(lapply(tLM, FUN = function(itLM){

      return(predY(predRE, data, time, itLM))

    }))

  }else{

    pred_Y <- unlist(lapply(tLM, FUN = function(itLM){

      return(max(0, predY(predRE, data, time, itLM) - threshold))

    }))

  }

  return(pred_Y)

}
