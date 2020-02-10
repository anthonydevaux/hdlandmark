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
LMsurv <- function(method = c("LM-cox","LM-rsf","cox","LM-coxnet","LM-splsDR"), models,
                   newdata, time, subject, marker_list, covar_list,
                   tHor){

  if (!all(method%in%c("LM-cox","LM-rsf","cox","LM-coxnet","LM-splsDR"))){
    stop("Only options available for method are 'LM-cox', 'LM-rsf', 'LM-coxnet', 'LM-splsDR' and 'cox'")
  }

  nb_method <- length(method)

  n <- nrow(newdata)

  pred_surv <- matrix(NA, nrow = n, ncol = nb_method,
                      dimnames = list(newdata[,subject], method))

  if (any(method == "LM-cox")){
    if (is.null(models[["LM-cox"]])){
      stop("No model found in models for method = LM-cox", "\n")
    }

    names.use <- names(newdata)[!(names(newdata) %in% marker_list)]

    newdata.LMcox <- newdata[,names.use]

    res_survfit <- survival::survfit(models[["LM-cox"]], newdata.LMcox)
    id_time <- sum(res_survfit$time <= tHor)

    if (n > 1){

      pred_surv[colnames(res_survfit$surv),"LM-cox"] <- res_survfit$surv[id_time,]

    }else{

      pred_surv[1,"LM-cox"] <- res_survfit$surv[id_time]

    }

  }

  if (any(method == "cox")){
    if (is.null(models[["cox"]])){
      stop("No model found in models for method = cox", "\n")
    }

    newdata.cox <- newdata[,c(covar_list, marker_list)]

    res_survfit <- survival::survfit(models[["cox"]], newdata.cox)
    id_time <- sum(res_survfit$time <= tHor)

    if (n > 1){

      pred_surv[colnames(res_survfit$surv),"cox"] <- res_survfit$surv[id_time,]

    }else{

      pred_surv[1,"cox"] <- res_survfit$surv[id_time]

    }

  }

  if (any(method == "LM-rsf")){
    if (is.null(models[["LM-rsf"]])){
      stop("No model found in models for method = LM-rsf", "\n")
    }

    names.use <- names(newdata)[!(names(newdata) %in% marker_list)]

    newdata.LMrsf <- newdata[,names.use]

    if (n > 1){

      res_survfit <- randomForestSRC::predict.rfsrc(models[["LM-rsf"]], newdata.LMrsf)
      id_time <- sum(res_survfit$time.interest <= tHor)
      formula.xvar <- as.formula(as.character(models[["LM-rsf"]]$call$formula)[c(1,3)])
      id_no_na <- rownames(model.frame(formula.xvar,
                                       newdata.LMrsf[,models[["LM-rsf"]]$xvar.names])) # id without NA
      pred_surv[id_no_na,"LM-rsf"] <- res_survfit$survival[,id_time]

    }else{

      if (any(is.na(newdata.LMrsf[,models[["LM-rsf"]]$xvar.names]))){

        pred_surv[1,"LM-rsf"] <- NA

      }else{

        res_survfit <- randomForestSRC::predict.rfsrc(models[["LM-rsf"]], newdata.LMrsf)
        id_time <- sum(res_survfit$time.interest <= tHor)
        pred_surv[1,"LM-rsf"] <- res_survfit$survival[id_time]

      }

    }

  }

  if (any(method == "LM-coxnet")){
    if (is.null(models[["LM-coxnet"]])){
      stop("No model found in models for method = LM-coxnet", "\n")
    }

    names.use <- names(newdata)[!(names(newdata) %in% c(subject, time, marker_list))]

    newdata.LMcoxnet <- na.omit(newdata[,names.use])

    newdata.LMcoxnet <- as.data.frame(model.matrix( ~ ., newdata.LMcoxnet)[,-1])

    res_survfit <- survival::survfit(models[["LM-coxnet"]], newdata = newdata.LMcoxnet)
    id_time <- sum(res_survfit$time <= tHor)

    if (n > 1){

      pred_surv[colnames(res_survfit$surv),"LM-coxnet"] <- res_survfit$surv[id_time,]

    }else{

      pred_surv[1,"LM-coxnet"] <- res_survfit$surv[id_time]

    }
  }

  if (any(method == "LM-splsDR")){
    if (is.null(models[["LM-splsDR"]])){
      stop("No model found in models for method = LM-splsDR", "\n")
    }

    names.use <- names(newdata)[!(names(newdata) %in% c(subject, time, marker_list))]

    newdata.LMsplsDR <- na.omit(newdata[,names.use])

    newdata.LMsplsDR <- as.data.frame(model.matrix( ~ ., newdata.LMsplsDR)[,-1])

    Xnames <- rownames(models[["LM-splsDR"]]$splsDR_modplsr$loadings$X)

    # centre reduit la matrice des nouveaux individus à partir du mean/sd du train
    Xh.scale <- t((t(newdata.LMsplsDR[,Xnames])-models[["LM-splsDR"]]$XplanCent[Xnames])/models[["LM-splsDR"]]$XplanScal[Xnames])

    X.spls <- matrix(NA, nrow = nrow(Xh.scale), ncol = ncol(models[["LM-splsDR"]]$tt_splsDR),
                     dimnames = list(rownames(Xh.scale), colnames(models[["LM-splsDR"]]$tt_splsDR)))

    u <- models[["LM-splsDR"]]$splsDR_modplsr$loadings$X

    X.spls[,1] <- Xh.scale%*%u[,1]

    if (ncol(X.spls) > 1){

      for (h in 2:ncol(X.spls)){

        th <- Xh.scale%*%u[,h-1]

        proj.num <- th%*%t(th)
        proj.den <- as.numeric(t(th)%*%th)
        proj <- proj.num / proj.den

        Xh.scale <- Xh.scale - proj%*%Xh.scale

        Xh <- Xh.scale%*%u[,h]

        X.spls[,h] <- Xh

      }

    }

    # prediction de la matrice X test dans le nouvel espace par spls
    res_survfit <- survival::survfit(models[["LM-splsDR"]]$cox_splsDR, newdata = as.data.frame(X.spls))
    id_time <- sum(res_survfit$time <= tHor)

    if (n > 1){

      pred_surv[colnames(res_survfit$surv),"LM-splsDR"] <- res_survfit$surv[id_time,]

    }else{

      pred_surv[1,"LM-splsDR"] <- res_survfit$surv[id_time]

    }

  }

  return(pred_surv)
}
