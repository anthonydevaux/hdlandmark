#' Function to compute several summaries of longitudinal outcomes
#'
#' @param params.true
#' @param data A dataframe containing longitudinal data
#' @param tLM An integer indicating the landmark time
#' @param subject A character indicating the subject variable
#' @param time A character indicating the time-dependent variable
#' @param RE
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
LMsum_simu <- function(markers.true, data, tLM, subject, time, RE, HW = tLM,
                       threshold = NULL){

  data_landmark <- data[which(data[,time]<tLM),]

  # build survival data

  data_surv <- data.frame()

  for (ind_subject in unique(data_landmark[,subject])){
    temp_subject <- data_landmark[which(data_landmark[,subject]==ind_subject),]
    temp_subject <- temp_subject[which.min(temp_subject[,time]),]

    data_surv <- rbind(data_surv, temp_subject)
  }

  marker_ind <- 1

  marker_list <- c()

  for (marker.true in markers.true){

    marker_name <- as.character(marker.true$model)[2]

    marker_list <- c(marker_list, marker_name)

    cat(paste0("Marker : ",marker_name), "\n")

    # random effects
    cat("random effect...")

    RE.marker <- list(b_i = RE[[marker_ind]],
               beta = marker.true$params$beta,
               sigmae = marker.true$params$sigmae,
               call = marker.true$model,
               group = subject)

    if (nrow(RE.marker$b_i)==0){
      marker_ind <- marker_ind + 1

      # RE
      data_surv[,(ncol(data_surv) + 1):(ncol(data_surv) + ncol(RE.marker$b_i))] <- NA
      b_i_var <- colnames(RE.marker$b_i)
      b_i_var_issue <- stringr::str_detect(b_i_var, regex("(?=\\().*?(?<=\\))")) # colnames contain parenthesis ?

      if (any(b_i_var_issue)){
        b_i_var[b_i_var_issue] <-
          regmatches(b_i_var[b_i_var_issue], gregexpr("(?<=\\().*?(?=\\))", b_i_var[b_i_var_issue], perl=T))[[1]]
      }

      colnames(data_surv)[(ncol(data_surv)-ncol(RE.marker$b_i)+1):(ncol(data_surv))] <-
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

      data_surv[,(ncol(data_surv) + 1):(ncol(data_surv) + ncol(RE.marker$b_i))] <- RE.marker$b_i

      b_i_var <- colnames(RE.marker$b_i)

      if (!is.null(b_i_var)){

        b_i_var_issue <- stringr::str_detect(b_i_var, regex("(?=\\().*?(?<=\\))")) # colnames contain parenthesis ?

        if (any(b_i_var_issue)){
          b_i_var[b_i_var_issue] <-
            regmatches(b_i_var[b_i_var_issue], gregexpr("(?<=\\().*?(?=\\))", b_i_var[b_i_var_issue], perl=T))[[1]]
        }

      }else{

        b_i_var <- seq(ncol(RE.marker$b_i)) - 1

      }

      colnames(data_surv)[(ncol(data_surv)-ncol(RE.marker$b_i)+1):(ncol(data_surv))] <-
        paste(marker_name, "RE", b_i_var, sep = "_")

      # prediction at landmark time
      cat("prediction...")
      pred_Y <- predY(RE.marker, data_surv, time, tLM)
      data_surv[which(data_surv[,subject]%in%rownames(pred_Y)),ncol(data_surv) + 1] <- pred_Y
      colnames(data_surv)[ncol(data_surv)] <- paste(marker_name, "pred", sep = "_")

      # slope at landmark time
      cat("slope...")
      deriv_Y <- derivY(RE.marker, data_surv, marker.true$deriv, time, tLM)
      data_surv[which(data_surv[,subject]%in%rownames(deriv_Y)),ncol(data_surv) + 1] <- deriv_Y
      colnames(data_surv)[ncol(data_surv)] <- paste(marker_name, "slope", sep = "_")

      # cumulative value at landmark time
      cat("cumulative...")
      cumul_Y <- cumulY(RE.marker, data_surv, time, marker_name, tLM, HW)
      data_surv[which(data_surv[,subject]%in%rownames(cumul_Y)),ncol(data_surv) + 1] <- cumul_Y
      colnames(data_surv)[ncol(data_surv)] <- paste(marker_name, "cumul", sep = "_")

      # threshold cumulative value at landmark time
      if (!is.null(threshold[[marker_name]])){

        cat("threshold...")
        threshold_value <- threshold[[marker_name]]

        threshold_Y <- thresholdY(RE.marker, data_surv, time, marker_name, tLM, threshold = threshold_value)
        data_surv[which(data_surv[,subject]%in%rownames(threshold_Y)),ncol(data_surv) + 1] <- threshold_Y
        colnames(data_surv)[ncol(data_surv)] <- paste(marker_name, "threshold", threshold_value, sep = "_")

      }
    }

    marker_ind <- marker_ind + 1
    cat("\n")

  }

  data_surv[,time] <- tLM

  attr(data_surv, "marker") <- marker_list

  cat("DONE!!!", "\n")

  return(data_surv)
}
