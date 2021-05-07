#' Function to compute the cumulative prediction value of longitudinal outcome using an history window
#'
#' @param predRE
#' @param data
#' @param time
#' @param tLM
#' @param HW
#'
#' @return
#' @export
#'
#' @examples
#'
cumulY <- function(predRE, data, time, tLM, HW){

  subject <- predRE$group

  ind_subjects <- as.integer(rownames(predRE$b_i))

  Y_cumul <- matrix(NA, nrow = length(ind_subjects), ncol = 1,
                    dimnames = list(ind_subjects, "Y_cumul"))

  Y_cumul_row <- 1

  for (ind_subject in ind_subjects){

    newdata_id <- data[data$id==ind_subject,]

    predRE_temp <- predRE

    predRE_temp$b_i <- predRE_temp$b_i[which(rownames(predRE_temp$b_i)==ind_subject),, drop = FALSE]

    predRE_temp$sigmae <- NULL # pas besoin en simulation pour approximer le cumul

    ###########################################

    predY.fct <- function(t){
      #as.numeric(predY_ind(tLM = t, predRE_temp, data, time))

      pred_Y <- unlist(lapply(t, FUN = function(itLM){

        return(predY(predRE_temp, data, time, itLM))

      }))

      return(as.numeric(pred_Y))

    }

    ##########################################

    Y_cumul[Y_cumul_row,] <- integrate(predY.fct, lower = tLM - HW, upper = tLM)$value

    Y_cumul_row <- Y_cumul_row + 1

  }

  return(Y_cumul)

}
