#' Function to compute several summaries of longitudinal outcomes
#'
#' @param lmm_objects A list containing as many \code{hlme} output as longitudinal outcome
#' @param data A dataframe containing longitudinal data
#' @param tLM An integer indicating the landmark time
#' @param subject A character indicating the subject variable
#' @param time A character indicating the time-dependent variable
#' @param derivForm_objects A list containing as many derivation form as longitudinal outcome
#' @param HW An integer meaning the range of history window for computing the cumulative values.
#' Default is \code{HW = tLM} resulting of the entire history.
#' @param threshold A list containing a numeric threshold
#'
#' @return The original dataframe including the summaries compute from longitudinal outcomes
#' @export
#'
#' @importFrom stringr str_detect regex
#'
#' @examples
#'
LMsum <- function(lmm_objects, data, tLM, subject, time, derivForm_objects, HW = tLM,
                  threshold = NULL){

  browser()

  data_landmark <- data[which(data[,time]<tLM),]

  # build survival data

  data_surv <- data.frame()

  for (ind_subject in unique(data_landmark[,subject])){
    temp_subject <- data_landmark[which(data_landmark[,subject]==ind_subject),]
    temp_subject <- temp_subject[which.min(temp_subject[,time]),]

    data_surv <- rbind(data_surv, temp_subject)
  }

  marker_ind <- 1

  for (lmm_object in lmm_objects){

    if (class(lmm_object)=="hlme"){

      marker_name <- as.character(lmm_object$call$fixed)[2]

    }else if (class(lmm_object)=="glmerMod"){

      marker_name <- as.character(formula(lmm_object))[2]

    }

    cat(paste0("Marker : ",marker_name), "\n")

    # random effects
    cat("random effect...")

    # marqueur continu
    if (class(data_landmark[,marker_name])=="numeric"){

      pred_RE <- predRE(lmm_object, data_landmark)

      if (nrow(pred_RE$b_i)==0){
        marker_ind <- marker_ind + 1

        # RE
        data_surv[,(ncol(data_surv) + 1):(ncol(data_surv) + ncol(pred_RE$b_i))] <- NA
        b_i_var <- colnames(pred_RE$b_i)
        b_i_var_issue <- stringr::str_detect(b_i_var, regex("(?=\\().*?(?<=\\))")) # colnames contain parenthesis ?

        if (any(b_i_var_issue)){
          b_i_var[b_i_var_issue] <-
            regmatches(b_i_var[b_i_var_issue], gregexpr("(?<=\\().*?(?=\\))", b_i_var[b_i_var_issue], perl=T))[[1]]
        }

        colnames(data_surv)[(ncol(data_surv)-ncol(pred_RE$b_i)+1):(ncol(data_surv))] <-
          paste(marker_name, "RE", b_i_var, sep = "_")

        # pred
        data_surv[,ncol(data_surv) + 1] <- NA
        colnames(data_surv)[ncol(data_surv)] <- paste(marker_name, "pred", sep = "_")

        # slope
        data_surv[,ncol(data_surv) + 1] <- NA
        colnames(data_surv)[ncol(data_surv)] <- paste(marker_name, "slope", sep = "_")

        # cumul
        data_surv[,ncol(data_surv) + 1] <- NA
        colnames(data_surv)[ncol(data_surv)] <- paste(marker_name, "cumul", sep = "_")

        # threshold
        if (!is.null(threshold[[marker_name]])){
          data_surv[,ncol(data_surv) + 1] <- NA
          colnames(data_surv)[ncol(data_surv)] <- paste(marker_name, "threshold", threshold_value, sep = "_")
        }

        next(paste("Unable to compute random effects for marker", marker_name))

      }else{

        data_surv[which(data_surv[,subject]%in%rownames(pred_RE$b_i)),
                  (ncol(data_surv) + 1):(ncol(data_surv) + ncol(pred_RE$b_i))] <- pred_RE$b_i

        b_i_var <- colnames(pred_RE$b_i)
        b_i_var_issue <- stringr::str_detect(b_i_var, regex("(?=\\().*?(?<=\\))")) # colnames contain parenthesis ?

        if (any(b_i_var_issue)){
          b_i_var[b_i_var_issue] <-
            regmatches(b_i_var[b_i_var_issue], gregexpr("(?<=\\().*?(?=\\))", b_i_var[b_i_var_issue], perl=T))[[1]]
        }

        colnames(data_surv)[(ncol(data_surv)-ncol(pred_RE$b_i)+1):(ncol(data_surv))] <-
          paste(marker_name, "RE", b_i_var, sep = "_")

        # prediction at landmark time
        cat("prediction...")
        pred_Y <- predY(pred_RE, data_surv, time, tLM)
        data_surv[which(data_surv[,subject]%in%rownames(pred_Y)),ncol(data_surv) + 1] <- pred_Y
        colnames(data_surv)[ncol(data_surv)] <- paste(marker_name, "pred", sep = "_")

        # slope at landmark time
        cat("slope...")
        deriv_Y <- derivY(pred_RE, data_surv, derivForm_objects[[marker_ind]])
        data_surv[which(data_surv[,subject]%in%rownames(deriv_Y)),ncol(data_surv) + 1] <- deriv_Y
        colnames(data_surv)[ncol(data_surv)] <- paste(marker_name, "slope", sep = "_")

        # cumulative value at landmark time
        cat("cumulative...")
        cumul_Y <- cumulY(pred_RE, data_surv, time, marker_name, tLM, HW)
        data_surv[which(data_surv[,subject]%in%rownames(cumul_Y)),ncol(data_surv) + 1] <- cumul_Y
        colnames(data_surv)[ncol(data_surv)] <- paste(marker_name, "cumul", sep = "_")

        # threshold cumulative value at landmark time
        if (!is.null(threshold[[marker_name]])){

          cat("threshold...")
          threshold_value <- threshold[[marker_name]]

          threshold_Y <- thresholdY(pred_RE, data_surv, time, marker_name, tLM, threshold = threshold_value)
          data_surv[which(data_surv[,subject]%in%rownames(threshold_Y)),ncol(data_surv) + 1] <- threshold_Y
          colnames(data_surv)[ncol(data_surv)] <- paste(marker_name, "threshold", threshold_value, sep = "_")

        }
      }

    # marqueur binaire
    }else if (class(data_landmark[,marker_name])=="factor"){

      pred_RE <- predRE.binary(lmm_object, data_landmark)

      if (nrow(pred_RE$b_i)==0){
        marker_ind <- marker_ind + 1

        # RE
        data_surv[,(ncol(data_surv) + 1):(ncol(data_surv) + ncol(pred_RE$b_i))] <- NA
        b_i_var <- colnames(pred_RE$b_i)
        b_i_var_issue <- stringr::str_detect(b_i_var, regex("(?=\\().*?(?<=\\))")) # colnames contain parenthesis ?

        if (any(b_i_var_issue)){
          b_i_var[b_i_var_issue] <-
            regmatches(b_i_var[b_i_var_issue], gregexpr("(?<=\\().*?(?=\\))", b_i_var[b_i_var_issue], perl=T))[[1]]
        }

        colnames(data_surv)[(ncol(data_surv)-ncol(pred_RE$b_i)+1):(ncol(data_surv))] <-
          paste(marker_name, "RE", b_i_var, sep = "_")

        # pred
        data_surv[,ncol(data_surv) + 1] <- NA
        colnames(data_surv)[ncol(data_surv)] <- paste(marker_name, "pred", sep = "_")

        # slope
        data_surv[,ncol(data_surv) + 1] <- NA
        colnames(data_surv)[ncol(data_surv)] <- paste(marker_name, "slope", sep = "_")

        # cumul
        data_surv[,ncol(data_surv) + 1] <- NA
        colnames(data_surv)[ncol(data_surv)] <- paste(marker_name, "cumul", sep = "_")

        # threshold
        if (!is.null(threshold[[marker_name]])){
          data_surv[,ncol(data_surv) + 1] <- NA
          colnames(data_surv)[ncol(data_surv)] <- paste(marker_name, "threshold", threshold_value, sep = "_")
        }

        next(paste("Unable to compute random effects for marker", marker_name))

      }else{

        data_surv[which(data_surv[,subject]%in%rownames(pred_RE$b_i)),
                  (ncol(data_surv) + 1):(ncol(data_surv) + ncol(pred_RE$b_i))] <- pred_RE$b_i

        b_i_var <- colnames(pred_RE$b_i)
        b_i_var_issue <- stringr::str_detect(b_i_var, regex("(?=\\().*?(?<=\\))")) # colnames contain parenthesis ?

        if (any(b_i_var_issue)){
          b_i_var[b_i_var_issue] <-
            regmatches(b_i_var[b_i_var_issue], gregexpr("(?<=\\().*?(?=\\))", b_i_var[b_i_var_issue], perl=T))[[1]]
        }

        colnames(data_surv)[(ncol(data_surv)-ncol(pred_RE$b_i)+1):(ncol(data_surv))] <-
          paste(marker_name, "RE", b_i_var, sep = "_")

      }

    }

    marker_ind <- marker_ind + 1
    cat("\n")

  }

  data_surv[,time] <- tLM

  cat("DONE!!!", "\n")

  return(data_surv)
}
