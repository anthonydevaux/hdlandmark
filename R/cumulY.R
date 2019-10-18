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

  subject <- predRE$formul$subject

  Y_cumul <- matrix(NA, nrow = length(unique(data[,subject])), ncol = 1,
                    dimnames = list(unique(data[,subject]), "Y_cumul"))

  Y_cumul_row <- 1

  for (ind_subject in unique(data[,subject])){

    newdata_id <- data[data$id==ind_subject,]

    predRE_temp <- predRE

    predRE_temp$b_i <- predRE_temp$b_i[which(rownames(predRE_temp$b_i)==ind_subject),, drop = FALSE]

    ###########################################

    F2 <- function(predRE_temp, data, time){
      f2 <- function(t){
        as.numeric(predY_ind(tLM = t, predRE_temp, data, time))
      }
      return(f2)
    }

    ##########################################

    Y_cumul[Y_cumul_row,] <- integrate(F2(predRE_temp, data, time), lower = tLM - HW, upper = tLM)$value

    Y_cumul_row <- Y_cumul_row + 1

  }

  return(Y_cumul)

}


# test <- function(tLM){
#   sapply(tLM, function(predRE, data, time){
#     predY(predRE, data, time)
#   })
# }
