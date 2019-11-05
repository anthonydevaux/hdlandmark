#' Title
#'
#' @param method
#' @param models
#' @param newdata
#' @param time
#' @param subject
#' @param tHor
#'
#' @return
#' @export
#'
#' @examples
#'
LMsurv <- function(method = c("LM-cox","LM-rsf","cox"), models,
                   newdata, time, subject,
                   tHor){

  if (!all(method%in%c("LM-cox","LM-rsf","cox"))){
    stop("Only options available for method are 'LM-cox', 'LM-rsf' and 'cox'")
  }

  nb_method <- length(method)

  n <- nrow(newdata)

  pred_surv <- matrix(NA, nrow = n, ncol = nb_method,
                      dimnames = list(newdata[,subject], method))

  if (any(method == "LM-cox")){
    if (is.null(models[["LM-cox"]])){
      stop("No model found in models for method = LM-cox", "\n")
    }

    newdata.LMcox <- subset(newdata, select = -c(serBilir2,serChol2,albumin,alkaline2,
                                                 SGOT2,platelets2,prothrombin2))

    res_survfit <- survfit(models[["LM-cox"]], newdata.LMcox)
    id_time <- sum(res_survfit$time <= tHor)

    if (n > 1){

      pred_surv[which(colnames(res_survfit$surv)%in%rownames(pred_surv)),"LM-cox"] <- res_survfit$surv[id_time,]

    }else{

      pred_surv[1,"LM-cox"] <- res_survfit$surv[id_time]

    }

  }

  if (any(method == "cox")){
    if (is.null(models[["cox"]])){
      stop("No model found in models for method = cox", "\n")
    }

    newdata.cox <- subset(newdata, select = c(drug,age,sex,serBilir2,serChol2,
                                              albumin,alkaline2,SGOT2,platelets2,prothrombin2))

    res_survfit <- survfit(models[["cox"]], newdata.cox)
    id_time <- sum(res_survfit$time <= tHor)

    if (n > 1){

      pred_surv[which(colnames(res_survfit$surv)%in%rownames(pred_surv)),"cox"] <- res_survfit$surv[id_time,]

    }else{

      pred_surv[1,"cox"] <- res_survfit$surv[id_time]

    }

  }

  if (any(method == "LM-rsf")){
    if (is.null(models[["LM-rsf"]])){
      stop("No model found in models for method = LM-rsf", "\n")
    }

    newdata.LMrsf <- subset(newdata, select = -c(serBilir2,serChol2,albumin,alkaline2,
                                                 SGOT2,platelets2,prothrombin2))

      if (n > 1){

        res_survfit <- predict.rfsrc(models[["LM-rsf"]], newdata.LMrsf)
        id_time <- sum(res_survfit$time.interest <= tHor)
        formula.xvar <- as.formula(as.character(model_obj$`LM-rsf`$call$formula)[c(1,3)])
        id_no_na <- rownames(model.frame(formula.xvar,
                                         newdata.LMrsf[,models[["LM-rsf"]]$xvar.names])) # id without NA
        pred_surv[id_no_na,"LM-rsf"] <- res_survfit$survival[,id_time]

      }else{

        if (any(is.na(newdata.LMrsf[,models[["LM-rsf"]]$xvar.names]))){

          pred_surv[1,"LM-rsf"] <- NA

        }else{

          res_survfit <- predict.rfsrc(models[["LM-rsf"]], newdata.LMrsf)
          id_time <- sum(res_survfit$time.interest <= tHor)
          pred_surv[1,"LM-rsf"] <- res_survfit$survival[id_time]

        }

      }

  }

  return(pred_surv)
}
