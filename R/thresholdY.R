#' Title
#'
#' @return
#' @export
#'
#' @examples
#'
thresholdY <- function(predRE, data, time, marker_name, tLM, threshold = NULL){

  subject <- predRE$formul$subject

  ind_subjects <- unique(data[!is.na(data[,marker_name]),subject])

  Y_threshold <- matrix(NA, nrow = length(ind_subjects), ncol = 1,
                    dimnames = list(ind_subjects, "Y_threshold"))

  Y_threshold_row <- 1

  for (ind_subject in ind_subjects){

    newdata_id <- data[data$id==ind_subject,]

    predRE_temp <- predRE

    predRE_temp$b_i <- predRE_temp$b_i[which(rownames(predRE_temp$b_i)==ind_subject),, drop = FALSE]

    ###########################################

    F2 <- function(predRE_temp, data, time){
      f2 <- function(t){
        as.numeric(predY_ind(tLM = t, predRE_temp, data, time, threshold))
      }
      return(f2)
    }

    ##########################################

    Y_threshold[Y_threshold_row,] <- integrate(F2(predRE_temp, data, time), lower = 0, upper = tLM)$value

    Y_threshold_row <- Y_threshold_row + 1

  }

  return(Y_threshold)

}
