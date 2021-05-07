#' Title
#'
#' @param data
#' @param data.surv
#' @param markers
#' @param tLM
#' @param subject
#' @param time
#' @param summaries
#' @param HW
#' @param threshold
#'
#' @return
#' @export
#'
#' @import lme4
#'
#' @examples
LMsum.glmm <- function(data, data.surv, markers, tLM, subject, time, summaries, HW, threshold = NULL){

  markers.names <- names(markers)

  for (marker in markers.names){ # loop marker

    cat(paste0("Marker : ",marker), "\n")

    ####### Random effects #######

    cat("random effect...")

    # numeric marker
    if (class(data[,marker])%in%c("numeric","integer")){

      res.lmm <- lmm.lcmm(markers[[marker]]$model, markers[[marker]]$formul, data)
      markers[[marker]]$formul <- markers[[marker]]$model
      markers[[marker]]$model <- res.lmm$model.output
      pred.RE <- res.lmm$pred.RE

    }

    # binary marker
    if (class(data[,marker])=="factor"){

      if (class(markers[[marker]]$model)!="glmerMod"){

        markers[[marker]]$model <- glmer(markers[[marker]]$model,
                                         data = data, family = "binomial",
                                         glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

      }

      pred.RE <- predRE.binary(markers[[marker]]$model, data)

    }

    if (any(summaries=="RE")){

      data.surv[which(data.surv[,subject]%in%rownames(pred.RE$b_i)),
                (ncol(data.surv) + 1):(ncol(data.surv) + ncol(pred.RE$b_i))] <- pred.RE$b_i

      colnames(data.surv)[(ncol(data.surv)-ncol(pred.RE$b_i)+1):(ncol(data.surv))] <-
        paste(marker, "RE", seq(ncol(pred.RE$b_i))-1, sep = "_")

    }

    ####### Current value at landmark time #######

    if (any(summaries=="pred")){

      cat("prediction...")

      pred_Y <- predY(pred.RE, data.surv, time, tLM)
      data.surv[which(data.surv[,subject]%in%rownames(pred_Y)),ncol(data.surv) + 1] <- pred_Y
      colnames(data.surv)[ncol(data.surv)] <- paste(marker, "pred", sep = "_")

    }

    ####### Current slope at landmark time #######

    if (any(summaries=="slope")){

      cat("slope...")

      deriv_Y <- derivY(pred.RE, data.surv, markers[[marker]]$deriv, time, tLM)
      data.surv[which(data.surv[,subject]%in%rownames(deriv_Y)),ncol(data.surv) + 1] <- deriv_Y
      colnames(data.surv)[ncol(data.surv)] <- paste(marker, "slope", sep = "_")

    }

    # Cumulative value at landmark time

    if (any(summaries=="cumulative")){

      cat("cumulative...")

      cumul_Y <- cumulY(predRE = pred.RE, data = data.surv, time = time, tLM = tLM, HW = HW)
      data.surv[which(data.surv[,subject]%in%rownames(cumul_Y)),ncol(data.surv) + 1] <- cumul_Y
      colnames(data.surv)[ncol(data.surv)] <- paste(marker, "cumul", sep = "_")

    }

    cat("\n")

  }

  return(list(data.surv = data.surv, markers.model = markers))

}
