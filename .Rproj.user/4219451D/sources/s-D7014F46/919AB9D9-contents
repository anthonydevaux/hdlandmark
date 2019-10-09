#' Function to compute several summaries of longitudinal outcomes
#'
#' @param lmm_objects A list containing as many \code{hlme} output as longitudinal outcome
#' @param data A dataframe containing longitudinal data
#' @param tLM An integer indicating the landmark time
#' @param subject A character indicating the subject variable
#' @param time A character indicating the time-dependent variable
#' @param derivForm_objects A list containing as many derivation form as longitudinal outcome
#'
#' @return The original dataframe including the summaries compute from longitudinal outcomes
#' @export
#'
#' @examples
#'
landmark <- function(lmm_objects, data, tLM, subject, time, derivForm_objects){

  data_landmark <- data[which(data[,time]<tLM),]

  # build survival data

  data_surv <- data.frame()

  for (ind_subject in unique(data_landmark[,subject])){
    temp_subject <- data_landmark[which(data_landmark[,subject]==ind_subject),]
    temp_subject <- temp_subject[which.min(temp_subject[,time]),]

    data_surv <- rbind(data_surv, temp_subject)
  }

  data_surv[,time] <- tLM

  marker_ind <- 1

  for (lmm_object in lmm_objects){

    marker_name <- as.character(lmm_object$call$fixed)[2]

    pred_RE <- predRE(lmm_object, data_landmark)

    pred_Y <- predY(pred_RE, data_surv)
    data_surv[which(data_surv[,subject]%in%rownames(pred_Y)),marker_name] <- predY

    deriv_Y <- derivY(pred_RE, data_surv, derivForm_objects[[marker_ind]])
    data_surv[which(data_surv[,subject]%in%rownames(deriv_Y)),ncol(data_surv) + 1] <- deriv_Y
    colnames(data_surv)[ncol(data_surv)] <- paste(marker_name, "slope", sep = "_")

    marker_ind <- marker_ind + 1

  }

  return(data_surv)
}
