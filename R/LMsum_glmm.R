#' Title
#'
#' @param data
#' @param data.surv
#' @param markers
#' @param tLM
#' @param subject
#' @param time
#' @param HW
#' @param threshold
#'
#' @return
#' @export
#'
#' @import lme4
#'
#' @examples
LMsum.glmm <- function(data, data.surv, markers, tLM, subject, time, threshold = NULL, lmm.package){

  markers.names <- names(markers)

  for (marker in markers.names){ # loop marker

    cat(paste0("Marker : ",marker), "\n")

    ####### Random effects #######

    cat("random effect...")

    browser()

    # numeric marker
    if (class(data[,marker])%in%c("numeric","integer")){

      if (lmm.package=="lme4"){

        res.lmm <- lmm.lme4(markers[[marker]]$model, data)
        markers[[marker]]$model <- res.lmm$model.output
        markers[[marker]]$model$package <- "lme4"

      }

      if (lmm.package=="lcmm"){

        res.lmm <- lmm.lcmm(markers[[marker]]$model, data)
        markers[[marker]]$model <- res.lmm$model.output
        markers[[marker]]$model$package <- "lcmm"

      }

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

    data.surv[which(data.surv[,subject]%in%rownames(pred.RE$b_i)),
              (ncol(data.surv) + 1):(ncol(data.surv) + ncol(pred.RE$b_i))] <- pred.RE$b_i

    # b_i_var <- colnames(pred.RE$b_i)
    # b_i_var_issue <- stringr::str_detect(b_i_var, regex("(?=\\().*?(?<=\\))")) # colnames contain parenthesis ?
    #
    # if (any(b_i_var_issue)){
    #   b_i_var[b_i_var_issue] <-
    #     regmatches(b_i_var[b_i_var_issue], gregexpr("(?<=\\().*?(?=\\))", b_i_var[b_i_var_issue], perl=T))[[1]]
    # }
    #
    # colnames(data.surv)[(ncol(data.surv)-ncol(pred.RE$b_i)+1):(ncol(data.surv))] <-
    #   paste(marker, "RE", b_i_var, sep = "_")

    colnames(data.surv)[(ncol(data.surv)-ncol(pred.RE$b_i)+1):(ncol(data.surv))] <-
      paste(marker, "RE", seq(ncol(pred.RE$b_i))-1, sep = "_")

    ####### Current value at landmark time #######

    cat("prediction...")

    pred_Y <- predY(pred.RE, data.surv, time, tLM)
    data.surv[which(data.surv[,subject]%in%rownames(pred_Y)),ncol(data.surv) + 1] <- pred_Y
    colnames(data.surv)[ncol(data.surv)] <- paste(marker, "pred", sep = "_")


    ####### Current slope at landmark time #######

    cat("slope...")

    deriv_Y <- derivY(pred.RE, data.surv, markers[[marker]]$deriv, time, tLM)
    data.surv[which(data.surv[,subject]%in%rownames(deriv_Y)),ncol(data.surv) + 1] <- deriv_Y
    colnames(data.surv)[ncol(data.surv)] <- paste(marker, "slope", sep = "_")

    # Cumulative value at landmark time

    cat("cumulative...", "\n")

    cumul_Y <- cumulY(pred.RE, data.surv, time, marker, tLM, HW = tLM)
    data.surv[which(data.surv[,subject]%in%rownames(cumul_Y)),ncol(data.surv) + 1] <- cumul_Y
    colnames(data.surv)[ncol(data.surv)] <- paste(marker, "cumul", sep = "_")

  }

  return(list(data.surv = data.surv, markers.model = markers))

}
