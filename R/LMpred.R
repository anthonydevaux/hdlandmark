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
LMpred <- function(method = c("LM-cox","LM-rsf","cox","LM-coxnet","LM-splsDR"), models,
                   newdata, time, subject, marker_list, covar_list,
                   tHor){

  if (!all(method%in%c("LM-cox", "LM-cox-noVS",
                       "LM-rsf", "LM-rsf-default", "LM-rsf-noVS",
                       "LM-rsf-BS", "LM-rsf-BS-default", "LM-rsf-BS-noVS",
                       "cox", "cox-noVS",
                       "LM-coxnet", "LM-coxnet-ridge", "LM-coxnet-lasso",
                       "LM-splsDR", "LM-splsDR-ridge", "LM-splsDR-lasso"))){
    stop("Only options available for method are 'LM-cox', 'LM-rsf', 'LM-coxnet', 'LM-splsDR' and 'cox'")
  }

  nb_method <- length(method)

  n <- nrow(newdata)

  pred_surv <- matrix(NA, nrow = n, ncol = nb_method,
                      dimnames = list(newdata[,subject], method))

  var_list_noLM <- c(covar_list, marker_list)
  var_list_LM <- names(newdata)[!(names(newdata) %in% marker_list)]

  # Cox

  if (any(method == "cox")){
    if (is.null(models[["cox"]])){
      stop("No model found in models for method = cox", "\n")
    }

    pred_surv <- pred.phm(model = models[["cox"]], method = "cox", newdata = newdata,
                          var_list = var_list_noLM, tHor = tHor, pred_surv = pred_surv)

  }

  # Cox sans selection de variables

  if (any(method == "cox-noVS")){
    if (is.null(models[["cox-noVS"]])){
      stop("No model found in models for method = cox-noVS", "\n")
    }

    pred_surv <- pred.phm(model = models[["cox-noVS"]], method = "cox-noVS", newdata = newdata,
                          var_list = var_list_noLM, tHor = tHor, pred_surv = pred_surv)

  }

  # LM-cox

  if (any(method == "LM-cox")){
    if (is.null(models[["LM-cox"]])){
      stop("No model found in models for method = LM-cox", "\n")
    }

    pred_surv <- pred.phm(model = models[["LM-cox"]], method = "LM-cox", newdata = newdata,
                          var_list = var_list_LM, tHor = tHor, pred_surv = pred_surv)

  }

  # LM-cox sans selection de variables

  if (any(method == "LM-cox-noVS")){
    if (is.null(models[["LM-cox-noVS"]])){
      stop("No model found in models for method = LM-cox-noVS", "\n")
    }

    pred_surv <- pred.phm(model = models[["LM-cox-noVS"]], method = "LM-cox-noVS", newdata = newdata,
                          var_list = var_list_LM, tHor = tHor, pred_surv = pred_surv)

  }

  # LM-rsf

  if (any(method == "LM-rsf")){
    if (is.null(models[["LM-rsf"]])){
      stop("No model found in models for method = LM-rsf", "\n")
    }

    pred_surv <- pred.rsf(model = models[["LM-rsf"]], method = "LM-rsf", newdata = newdata,
                          var_list = var_list_LM, tHor = tHor, pred_surv = pred_surv)

  }

  # LM-rsf default

  if (any(method == "LM-rsf-default")){
    if (is.null(models[["LM-rsf-default"]])){
      stop("No model found in models for method = LM-rsf-default", "\n")
    }

    pred_surv <- pred.rsf(model = models[["LM-rsf-default"]], method = "LM-rsf-default", newdata = newdata,
                          var_list = var_list_LM, tHor = tHor, pred_surv = pred_surv)

  }

  # LM-rsf sans selection de variables

  if (any(method == "LM-rsf-noVS")){
    if (is.null(models[["LM-rsf-noVS"]])){
      stop("No model found in models for method = LM-rsf-noVS", "\n")
    }

    pred_surv <- pred.rsf(model = models[["LM-rsf-noVS"]], method = "LM-rsf-noVS", newdata = newdata,
                          var_list = var_list_LM, tHor = tHor, pred_surv = pred_surv)

  }

  # LM-rsf - split bs.gradient

  if (any(method == "LM-rsf-BS")){
    if (is.null(models[["LM-rsf-BS"]])){
      stop("No model found in models for method = LM-rsf-BS", "\n")
    }

    pred_surv <- pred.rsf(model = models[["LM-rsf-BS"]], method = "LM-rsf-BS", newdata = newdata,
                          var_list = var_list_LM, tHor = tHor, pred_surv = pred_surv)

  }

  # LM-rsf default - split bs.gradient

  if (any(method == "LM-rsf-BS-default")){
    if (is.null(models[["LM-rsf-BS-default"]])){
      stop("No model found in models for method = LM-rsf-BS-default", "\n")
    }

    pred_surv <- pred.rsf(model = models[["LM-rsf-BS-default"]], method = "LM-rsf-BS-default", newdata = newdata,
                          var_list = var_list_LM, tHor = tHor, pred_surv = pred_surv)

  }

  # LM-rsf sans selection de variables - split bs.gradient

  if (any(method == "LM-rsf-BS-noVS")){
    if (is.null(models[["LM-rsf-BS-noVS"]])){
      stop("No model found in models for method = LM-rsf-BS-noVS", "\n")
    }

    pred_surv <- pred.rsf(model = models[["LM-rsf-BS-noVS"]], method = "LM-rsf-BS-noVS", newdata = newdata,
                          var_list = var_list_LM, tHor = tHor, pred_surv = pred_surv)

  }

  ##

  # LM-coxnet

  if (any(method == "LM-coxnet")){
    if (is.null(models[["LM-coxnet"]])){
      stop("No model found in models for method = LM-coxnet", "\n")
    }

    # names.use <- names(newdata)[!(names(newdata) %in% c(subject, time, marker_list))]
    # newdata.LMcoxnet <- na.omit(newdata[,names.use])
    # newdata.LMcoxnet <- as.data.frame(model.matrix( ~ ., newdata.LMcoxnet)[,-1])

    newdata.coxnet <- as.data.frame(model.matrix( ~ ., na.omit(newdata[,var_list_LM]))[,-1])

    pred_surv <- pred.phm(model = models[["LM-coxnet"]], method = "LM-coxnet",
                          newdata = newdata.coxnet, var_list = NULL, tHor = tHor,
                          pred_surv = pred_surv)

  }

  # LM-coxnet-ridge

  if (any(method == "LM-coxnet-ridge")){
    if (is.null(models[["LM-coxnet-ridge"]])){
      stop("No model found in models for method = LM-coxnet-ridge", "\n")
    }

    # names.use <- names(newdata)[!(names(newdata) %in% c(subject, time, marker_list))]
    # newdata.LMcoxnet <- na.omit(newdata[,names.use])
    # newdata.LMcoxnet <- as.data.frame(model.matrix( ~ ., newdata.LMcoxnet)[,-1])

    newdata.coxnet <- as.data.frame(model.matrix( ~ ., na.omit(newdata[,var_list_LM]))[,-1])

    pred_surv <- pred.phm(model = models[["LM-coxnet-ridge"]], method = "LM-coxnet-ridge",
                          newdata = newdata.coxnet, var_list = NULL, tHor = tHor,
                          pred_surv = pred_surv)
  }

  # LM-coxnet-lasso

  if (any(method == "LM-coxnet-lasso")){
    if (is.null(models[["LM-coxnet-lasso"]])){
      stop("No model found in models for method = LM-coxnet-lasso", "\n")
    }

    # names.use <- names(newdata)[!(names(newdata) %in% c(subject, time, marker_list))]
    # newdata.LMcoxnet <- na.omit(newdata[,names.use])
    # newdata.LMcoxnet <- as.data.frame(model.matrix( ~ ., newdata.LMcoxnet)[,-1])

    newdata.coxnet <- as.data.frame(model.matrix( ~ ., na.omit(newdata[,var_list_LM]))[,-1])

    pred_surv <- pred.phm(model = models[["LM-coxnet-lasso"]], method = "LM-coxnet-lasso",
                          newdata = newdata.coxnet, var_list = NULL, tHor = tHor,
                          pred_surv = pred_surv)
  }

  # LM-splsDR

  if (any(method == "LM-splsDR")){
    if (is.null(models[["LM-splsDR"]])){
      stop("No model found in models for method = LM-splsDR", "\n")
    }

    # names.use <- names(newdata)[!(names(newdata) %in% c(subject, time, marker_list))]
    # newdata.LMsplsDR <- na.omit(newdata[,names.use])
    # newdata.LMsplsDR <- as.data.frame(model.matrix( ~ ., newdata.LMsplsDR)[,-1])

    newdata.splscox <- as.data.frame(model.matrix( ~ ., na.omit(newdata[,var_list_LM]))[,-1])

    pred_surv <- pred.splscox(model = models[["LM-splsDR"]], method = "LM-splsDR",
                              newdata = newdata.splscox, tHor = tHor, pred_surv = pred_surv)

  }

  # LM-splsDR-ridge

  if (any(method == "LM-splsDR-ridge")){
    if (is.null(models[["LM-splsDR-ridge"]])){
      stop("No model found in models for method = LM-splsDR-ridge", "\n")
    }

    # names.use <- names(newdata)[!(names(newdata) %in% c(subject, time, marker_list))]
    # newdata.LMsplsDR <- na.omit(newdata[,names.use])
    # newdata.LMsplsDR <- as.data.frame(model.matrix( ~ ., newdata.LMsplsDR)[,-1])

    newdata.splscox <- as.data.frame(model.matrix( ~ ., na.omit(newdata[,var_list_LM]))[,-1])

    pred_surv <- pred.splscox(model = models[["LM-splsDR-ridge"]], method = "LM-splsDR-ridge",
                              newdata = newdata.splscox, tHor = tHor, pred_surv = pred_surv)

  }

  # LM-splsDR-lasso

  if (any(method == "LM-splsDR-lasso")){
    if (is.null(models[["LM-splsDR-lasso"]])){
      stop("No model found in models for method = LM-splsDR-lasso", "\n")
    }

    # names.use <- names(newdata)[!(names(newdata) %in% c(subject, time, marker_list))]
    # newdata.LMsplsDR <- na.omit(newdata[,names.use])
    # newdata.LMsplsDR <- as.data.frame(model.matrix( ~ ., newdata.LMsplsDR)[,-1])

    newdata.splscox <- as.data.frame(model.matrix( ~ ., na.omit(newdata[,var_list_LM]))[,-1])

    pred_surv <- pred.splscox(model = models[["LM-splsDR-lasso"]], method = "LM-splsDR-lasso",
                              newdata = newdata.splscox, tHor = tHor, pred_surv = pred_surv)

  }

  return(pred_surv)
}
