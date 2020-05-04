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
#' @param long.method
#' @param surv.methods
#' @param cox.autoVar
#' @param cox.allVar
#' @param coxnet.opt
#' @param coxnet.lasso
#' @param coxnet.ridge
#' @param spls.opt
#' @param spls.nosparse
#' @param spls.maxsparse
#' @param rsf.split
#' @param rsf.opt
#' @param rsf.noVS
#' @param rsf.default
#'
#' @return
#' @export
#'
#' @examples
hdlandmark <- function(data, data.pred = NULL, markers, tLMs, tHors,
                       subject, time, time.event, event,
                       long.method = c("combine", "GLMM", "MFPC"),
                       surv.methods = c("cox", "penalized-cox", "sPLS", "rsf"),
                       kfolds = 10, seed = 1234, parallel = FALSE, Ncpus = 1,
                       cox.autoVar = TRUE, cox.allVar = !cox.autoVar,
                       coxnet.opt = TRUE, coxnet.lasso = !coxnet.opt, coxnet.ridge = !coxnet.opt,
                       spls.opt = TRUE, spls.nosparse = !spls.opt, spls.maxsparse = !spls.opt,
                       rsf.split = c("logrank", "bs.gradient"), rsf.opt = TRUE, rsf.noVS = !rsf.opt, rsf.default = !rsf.opt){

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

  resu <- list()

  for (tLM in tLMs){ # landmark time loop

    data.tLM <- data[which(data[,time]<tLM&data[,time.event]>tLM),]

    if (!is.null(data.pred)){ # different data training and test

      data.pred.tLM <- data.pred[which(data.pred[,time]<tLM&data.pred[,time.event]>tLM),]
      kfolds <- 1

    }else{ # same data estimation/prediction

      data.pred.tLM <- data.tLM

    }

    if (parallel){



    }else{

      resu[[as.character(tLM)]] <- .hdlandmark(data = data.tLM, data.pred = data.pred.tLM, markers = markers, tLM = tLM, tHors = tHors,
                                               kfolds = kfolds, seed = seed,
                                               subject = subject, time = time, time.event = time.event, event = event,
                                               long.method = long.method, surv.methods = surv.methods,
                                               cox.autoVar = cox.autoVar, cox.allVar = cox.allVar,
                                               coxnet.opt = coxnet.opt, coxnet.lasso = coxnet.lasso, coxnet.ridge = coxnet.ridge,
                                               spls.opt = spls.opt, spls.nosparse = spls.nosparse, spls.maxsparse = spls.maxsparse,
                                               rsf.split = rsf.split, rsf.opt = rsf.opt, rsf.noVS = rsf.noVS, rsf.default = rsf.default)

    }

  }

  cat("DONE !", "\n")

  return(list(tLMs = tLMs, tHors = tHors, resu = resu,
              long.method = long.method, surv.methods = surv.methods,
              cox.autoVar = cox.autoVar, cox.allVar = cox.allVar,
              coxnet.opt = coxnet.opt, coxnet.lasso = coxnet.lasso,
              coxnet.ridge = coxnet.ridge,
              spls.opt = spls.opt, spls.nosparse = spls.nosparse,
              spls.maxsparse = spls.maxsparse,
              rsf.split = rsf.split, rsf.opt = rsf.opt, rsf.noVS = rsf.noVS,
              rsf.default = rsf.default))

}



.hdlandmark <- function(data, data.pred, markers, tLM, tHors,
                        kfolds, seed,
                        subject, time, time.event, event,
                        long.method, surv.methods,
                        cox.autoVar, cox.allVar,
                        coxnet.opt, coxnet.lasso, coxnet.ridge,
                        spls.opt, spls.nosparse, spls.maxsparse,
                        rsf.split, rsf.opt, rsf.noVS, rsf.default){

  ids <- unique(data[,subject])
  n <- length(ids)
  set.seed(seed)
  fold <- sample(rep(1:kfolds, length.out = n))

  pred.surv <- AUC <- BS <- MSEP <- list()

  for (k in 1:kfolds){

    cat(paste0("Fold : ",k,"/",kfolds),"\n")

    ids.test <- ids[which(fold==k)]

    if (kfolds==1){

      ids.train <- ids.test

    }else{

      ids.train <- ids[which(fold!=k)]

    }

    data.k <- data[which(data[,subject]%in%ids.train),]
    data.pred.k <- data.pred[which(data.pred[,subject]%in%ids.test),]

    # estimation of summaries on training data and test data

    res.LMsum <- LMsummaries(data = data.k, data.pred = data.pred.k, markers = marker, tLM = tLM,
                             subject = subject, time = time, time.event = time.event, event = event,
                             long.method = long.method)

    # survival model on training data

    res.LMsurv <- LMsurv(data.surv = res.LMsum$data.surv, long.method = long.method, surv.methods = surv.methods,
                         cox.autoVar = cox.autoVar, cox.allVar = cox.allVar,
                         coxnet.opt = coxnet.opt, coxnet.lasso = coxnet.lasso, coxnet.ridge = coxnet.ridge,
                         spls.opt = spls.opt, spls.nosparse = spls.nosparse, spls.maxsparse = spls.maxsparse,
                         rsf.split = rsf.split, rsf.opt = rsf.opt, rsf.noVS = rsf.noVS, rsf.default = rsf.default)

    for (tHor in tHors){ # tHor loop

      # survival model on test data

      pred.surv.tHor <- LMpred(data.surv = res.LMsum$data.surv.pred, model.surv = res.LMsurv$model.surv,
                               long.method = long.method, surv.methods = surv.methods,
                               tHor = tHor)

      res.LMassess <- LMassess(pred.surv = pred.surv.tHor, data.surv = res.LMsum$data.surv.pred, tHor = tHor)

      AUC[[as.character(tHor)]] <- rbind(AUC[[as.character(tHor)]], res.LMassess$AUC)
      BS[[as.character(tHor)]] <- rbind(BS[[as.character(tHor)]], res.LMassess$BS)
      MSEP[[as.character(tHor)]] <- rbind(MSEP[[as.character(tHor)]], res.LMassess$MSEP)

      pred.surv[[as.character(tHor)]] <- rbind(pred.surv[[as.character(tHor)]], pred.surv.tHor)
      pred.surv[[as.character(tHor)]] <- pred.surv[[as.character(tHor)]][order(as.integer(rownames(pred.surv[[as.character(tHor)]]))), ]

    }

  }

  return(list(data.surv = res.LMsum$data.surv, data.surv.pred = res.LMsum$data.surv.pred,
              model.surv = res.LMsurv$model.surv, pred.surv = pred.surv,
              AUC = AUC, BS = BS, MSEP = MSEP))

}
