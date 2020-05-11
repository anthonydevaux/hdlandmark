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
#' @param cox.submodels
#' @param coxnet.submodels
#' @param spls.submodels
#' @param rsf.submodels
#' @param rsf.split
#' @param kfolds
#' @param seed
#' @param parallel
#' @param Ncpus
#'
#' @return
#' @export
#'
#' @examples
hdlandmark <- function(data, data.pred = NULL, markers, tLMs, tHors,
                       subject, time, time.event, event,
                       long.method = c("combine", "GLMM", "MFPC"),
                       cox.submodels = c("autoVar","allVar"), coxnet.submodels = c("opt","lasso","ridge"),
                       spls.submodels = c("opt","nosparse","maxsparse"), rsf.submodels = c("opt","noVS","default"),
                       rsf.split = c("logrank", "bs.gradient"),
                       kfolds = 10, seed = 1234, parallel = FALSE, Ncpus = 1){

  ####### Check #######

  if (!class(data)%in%c("data.frame","matrix")){
    stop("data should be class of data.frame or matrix")
  }
  if (!class(markers)=="list"){
    stop("markers should be class of list")
  }
  if (!all(names(markers)%in%colnames(data))){
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
  if (length(long.method)>1){
    stop("long.method should be length of 1")
  }else{
    if (!all(long.method%in%c("combine","GLMM","MFPC"))){
      stop("Only combine, GLMM or MFPC are allowed for long.method")
    }
  }
  if (!all(cox.submodels%in%c("autoVar","allVar"))){
    stop("Only autoVar or allVar are allowed for cox.submodels")
  }
  if (!all(coxnet.submodels%in%c("opt","lasso","ridge"))){
    stop("Only opt, lasso or ridge are allowed for coxnet.submodels")
  }
  if (!all(spls.submodels%in%c("opt","nosparse","maxsparse"))){
    stop("Only opt, nosparse or maxsparse are allowed for spls.submodels")
  }
  if (!all(rsf.submodels%in%c("opt","noVS","default","ranger"))){
    stop("Only opt, noVS, default and ranger are allowed for rsf.submodels")
  }
  if (!all(rsf.split%in%c("logrank","bs.gradient"))){
    stop("Only logrank or bs.gradient are allowed for rsf.split")
  }

  if (!is.null(data.pred)){
    if (!class(data.pred)%in%c("data.frame","matrix")){
      stop("data.pred should be class of data.frame or matrix")
    }
    if (!all(names(markers)%in%colnames(data.pred))){
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

  surv.methods <- NULL

  if (length(cox.submodels)>0){
    surv.methods <- c(surv.methods, "cox")
  }
  if (length(coxnet.submodels)>0){
    surv.methods <- c(surv.methods, "penalized-cox")
  }
  if (length(spls.submodels)>0){
    surv.methods <- c(surv.methods, "spls")
  }
  if (length(rsf.submodels)>0){
    surv.methods <- c(surv.methods, "rsf")
  }

  if (length(surv.methods)==0){
    stop("No method has been selected ! Please select at least a method of cox.submodels, coxnet.submodels,
         spls.submodels or rsf.submodels")
  }

  ##############################################

  models <- list()

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

      models[[as.character(tLM)]] <- .hdlandmark(data = data.tLM, data.pred = data.pred.tLM, markers = markers, tLM = tLM, tHors = tHors,
                                                 kfolds = kfolds, seed = seed,
                                                 subject = subject, time = time, time.event = time.event, event = event,
                                                 long.method = long.method, surv.methods = surv.methods,
                                                 cox.submodels = cox.submodels, coxnet.submodels = coxnet.submodels,
                                                 spls.submodels = spls.submodels, rsf.submodels = rsf.submodels,
                                                 rsf.split = rsf.split)

    }

  }

  cat("DONE !", "\n")

  resu <- list(tLMs = tLMs, tHors = tHors, models = models,
               long.method = long.method, surv.methods = surv.methods,
               models.name = colnames(models[[1]]$pred.surv[[1]]), kfolds = kfolds)

  class(resu) <- "hdlandmark"

  return(resu)

}



.hdlandmark <- function(data, data.pred, markers, tLM, tHors,
                        kfolds, seed,
                        subject, time, time.event, event,
                        long.method, surv.methods,
                        cox.submodels, coxnet.submodels,
                        spls.submodels, rsf.submodels,
                        rsf.split){

  ids <- unique(data[,subject])
  n <- length(ids)
  set.seed(seed)
  fold <- sample(rep(1:kfolds, length.out = n))

  pred.surv <- AUC <- BS <- list()

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

    res.LMsum <- LMsummaries(data = data.k, data.pred = data.pred.k, markers = markers, tLM = tLM,
                             subject = subject, time = time, time.event = time.event, event = event,
                             long.method = long.method)

    for (tHor in tHors){ # tHor loop

      # censuring to horizon time
      data.surv <- res.LMsum$data.surv
      data.surv[which(data.surv$time.event > tHor), "event"] <- 0
      data.surv$time.event <- pmin(data.surv$time.event, tHor)

      data.surv.pred <- res.LMsum$data.surv.pred
      data.surv.pred[which(data.surv.pred$time.event > tHor), "event"] <- 0
      data.surv.pred$time.event <- pmin(data.surv.pred$time.event, tHor)

      # survival model on training data
      res.LMsurv <- LMsurv(data.surv = data.surv, surv.methods = surv.methods,
                           cox.submodels = cox.submodels, coxnet.submodels = coxnet.submodels,
                           spls.submodels = spls.submodels, rsf.submodels = rsf.submodels,
                           rsf.split = rsf.split)

      # survival model on test data
      res.LMpred <- LMpred(data.surv = data.surv.pred, model.surv = res.LMsurv$model.surv,
                           long.method = long.method, surv.methods = surv.methods,
                           tHor = tHor)

      res.LMassess <- LMassess(pred.surv = res.LMpred$pred.surv, data.surv = data.surv.pred, tHor = tHor)

      AUC[[as.character(tHor)]] <- rbind(AUC[[as.character(tHor)]], res.LMassess$AUC)
      BS[[as.character(tHor)]] <- rbind(BS[[as.character(tHor)]], res.LMassess$BS)

      pred.surv[[as.character(tHor)]] <- rbind(pred.surv[[as.character(tHor)]], res.LMpred$pred.surv)
      pred.surv[[as.character(tHor)]] <- pred.surv[[as.character(tHor)]][order(as.integer(rownames(pred.surv[[as.character(tHor)]]))), ]

    }

  }

  return(list(data.surv = res.LMsum$data.surv, data.surv.pred = res.LMsum$data.surv.pred,
              model.surv = res.LMsurv$model.surv, pred.surv = pred.surv,
              AUC = AUC, BS = BS))

}
