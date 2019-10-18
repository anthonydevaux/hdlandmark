#' Function to predict the value of longitudinal outcome for a specific time
#'
#' @param predRE A list object from \code{predRE} function
#' @param data A dataframe where each row containing some predictive variables for a specific subject
#'
#' @return A matrix containing the prediction value of the longitudinal outcome for each subject
#'
#' @examples
#'
predY_ind <- function(predRE, data, time, tLM){

  pred_Y <- c()

  for (i in 1:length(tLM)){
    pred_Y <- c(pred_Y, predY(predRE, data, time, tLM[i]))
  }

  return(pred_Y)

}
