#' Title
#'
#' @param data
#' @param data.pred
#' @param markers
#' @param tLMs
#' @param tHors
#' @param subject
#' @param time
#' @param time.event
#' @param event
#' @param long.methods
#' @param surv.methods
#'
#' @return
#' @export
#'
#' @examples
hdlandmark <- function(data, data.pred = NULL, markers, tLMs, tHors,
                       subject, time, time.event, event,
                       long.methods = c("combine", "both", "GLMM", "MFPC"),
                       surv.methods = c("cox", "penalized-cox", "sPLS", "rsf")){

  ####### Check #######

  if (!class(data)%in%c("data.frame","matrix")){
    stop("data should be class of data.frame or matrix")
  }
  if (!class(markers)=="list"){
    stop("markers should be class of list")
  }
  if (!all(names(marker)%in%colnames(data))){
    stop("At least one marker variable is missing in data")
  }
  if (!class(tLMs)=="numeric"){
    stop("tLM should be class of numeric")
  }
  if (!all(tLMs>0)){
    stop("tLMs should be positive")
  }
  if (!class(tHors)=="numeric"){
    stop("tHor should be class of list")
  }
  if (!all(tHors>0)){
    stop("tHors should be positive")
  }
  if (!class(subject)=="character"){
    stop("subject should be class of character")
  }
  if (!subject%in%colnames(data)){
    stop("subject variable is missing in data")
  }
  if (!class(data[,subject])=="integer"){
    data[,subject] <- as.integer(data[,subject])
  }
  if (!class(time)=="character"){
    stop("time should be class of character")
  }
  if (!time%in%colnames(data)){
    stop("time variable is missing in data")
  }
  if (!class(time.event)=="character"){
    stop("time.event should be class of character")
  }
  if (!time.event%in%colnames(data)){
    stop("time.event variable is missing in data")
  }
  if (!class(event)=="character"){
    stop("event should be class of character")
  }
  if (!event%in%colnames(data)){
    stop("event variable is missing in data")
  }

  if (!is.null(data.pred)){
    if (!class(data.pred)%in%c("data.frame","matrix")){
      stop("data.pred should be class of data.frame or matrix")
    }
    if (!all(names(marker)%in%colnames(data.pred))){
      stop("At least one marker variable is missing in data.pred")
    }
    if (!subject%in%colnames(data.pred)){
      stop("subject variable is missing in data.pred")
    }
    if (!class(data.pred[,subject])=="integer"){
      data.pred[,subject] <- as.integer(data.pred[,subject])
    }
    if (!time%in%colnames(data.pred)){
      stop("time variable is missing in data.pred")
    }
  }

  ##############################################

  for (tLM in tLMs){ # tLM loop

    data.tLM <- data[which(data[,time]<tLM&data[,time.event]>tLM),]

    if (!is.null(data.pred)){ # different data training and test

      data.pred.tLM <- data.pred[which(data.pred[,time]<tLM&data.pred[,time.event]>tLM),]

    }else{ # same data estimation/prediction

      data.pred.tLM <- data.tLM

    }

    for (tHor in tHors){ # tHor loop

      res <- .hdlandmark(data = data.tLM, data.pred = data.pred.tLM, markers = markers, tLM = tLM, tHor = tHor,
                         subject = subject, time = time, time.event = time.event, event = event,
                         long.methods = long.methods,
                         surv.methods = surv.methods)

    }

  }

  return(res)

}

.hdlandmark <- function(data, data.pred, markers, tLM, tHor, subject, time, time.event, event,
                        long.methods, surv.methods){

  # estimation of summaries on training data and test data

  res.LMsum <- LMsummaries(data = data, data.pred = data.pred, markers = marker, tLM = tLM,
                           subject = subject, time = time, time.event = time.event, event = event,
                           long.methods = long.methods)

  # survival model on training data

  res.LMsurv <- LMsurv(data.surv = res.LMsum$data.surv, long.methods = long.methods, surv.methods = surv.methods)

  # survival model on test data

  res.LMpred <- LMpred(data.surv = res.LMsum$data.surv.pred, model.surv = res.LMsurv$model.surv,
                       long.methods = long.methods, surv.methods = surv.methods,
                       tHor = tHor)

}
