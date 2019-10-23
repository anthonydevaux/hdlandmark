#' Function to predict the value of longitudinal outcome for a specific time
#'
#' @param predRE A list object from \code{predRE} function
#' @param data A dataframe where each row containing some predictive variables for a specific subject
#'
#' @return A matrix containing the prediction value of the longitudinal outcome for each subject
#'
#' @examples
#'
predY_ind <- function(predRE, data, time, tLM, threshold = NULL){

  pred_Y <- c()

  if (is.null(threshold)){

    for (i in 1:length(tLM)){

      predY_value <- predY(predRE, data, time, tLM[i])

      pred_Y <- c(pred_Y, predY_value)

    }

  }else{

    for (i in 1:length(tLM)){

      predY_value <- max(0, predY(predRE, data, time, tLM[i]) - threshold)

      pred_Y <- c(pred_Y, predY_value)

    }

  }

  return(pred_Y)

}
