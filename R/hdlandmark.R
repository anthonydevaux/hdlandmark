#' Compute individual dynamic prediction of clinical endpoint using large dimensional longitudinal biomarker history
#'
#' @description
#'
#' hdlandmark provides individual survival probabilities using covariates and summaries build
#' on longitudinal data from biomarkers collected over the time.
#' For each biomarker, an ensemble of predictive summaries are computed at the user-specified landmark time \code{tLM}.
#' For instance, we use random-effects, level, slope and cumulative level. Then, these summaries and covariates are used as input in several survival prediction methods including: Cox model (his extension with penalty), sparse-Partial Least Square for survival data and random survival forests
#' For each survival prediction method, we provide the individual prediction on horizon time \code{tHor}.
#'
#' @param data data.frame object containing longitudinal and survival data
#' @param data.pred (optional) data.frame object for predictions. If missing, \code{data.pred} is made using \code{kfolds} cross-validation
#' @param markers list containing the modeling of repeated measures for each marker
#' @param tLMs numeric vector of landmark times
#' @param tHors numeric vector of horizon times
#' @param subject variable name in data (and \code{data.pred}) that identifies the different subjects
#' @param time variable name in data (and \code{data.pred}) which contains time measurements
#' @param time.event variable name in data (and \code{data.pred}) which contains time-to-event
#' @param event variable name in data (and \code{data.pred}) which contains time-to-event
#' @param long.method character that specifies how to model the longitudinal data. Choices are \code{GLMM} for generalized mixed model \insertCite{laird_random-effects_1982}{hdlandmark},
#' \code{MFPC} for multivariate functional principal components \insertCite{yao_functional_2005}{hdlandmark} (works only on continuous markers) or \code{combine} for both.
#' @param surv.covar covariates measure at \code{baseline} or last observation before landmark time \code{LOtLM}
#' @param cox.submodels a character vector containing Cox submodels \insertCite{cox_regression_1972}{hdlandmark}. \code{autoVar} for Cox with backward variable selection. \code{allVar} for Cox with all variables
#' @param coxnet.submodels a character vector containing penalized Cox submodels \insertCite{simon_regularization_2011}{hdlandmark}. \code{opt} for tuning the elastic net parameter penalty, \code{lasso} for lasso penalty and \code{ridge} for ridge penalty.
#' @param penaFG.submodels
#' @param spls.submodels a character vector containing Deviance residuals sparse-Partial Least Square sub-methods \insertCite{bastien_deviance_2015}{hdlandmark}. \code{opt} for tuning sparcity parameter \eqn{\eta}, \code{nosparse} for \eqn{\eta = 0} and \code{maxsparse} for \eqn{\eta = 0.9} \insertCite{@see also @chun_sparse_2010}{hdlandmark}
#' @param rsf.submodels a character vector containing random survival forests sub-methods \insertCite{ishwaran_random_2008}{hdlandmark}.
#' @param rsf.split a character vector containing the split criterion for random survival forests sub-methods. \code{logrank} for log-rank splitting or \code{bs.gradient} for gradient-based brier score splitting.
#' @param cause
#' @param HW
#' @param summaries
#' @param kfolds number of fold in cross-validation
#' @param seed (optional) seed number
#' @param scaling boolean to scale summaries (default is \code{FALSE})
#' @param SL.weights (optional) allow to compute individual probabilities from a superlearner using numeric vector of weights for each sub-methods
#'
#' @import Rdpack
#'
#' @return
#' \item{tLMs}{landmark time(s)}
#' \item{tHors}{horizon time(s)}
#' \item{models}{a list for each landmark time(s):}
#' \itemize{
#'  \item \code{data.surv} input data in survival methods for training (only available on the last fold)
#'  \item \code{data.surv.pred} input data in survival methods for predicting (only available on the last fold)
#'  \item \code{model.surv} output object for the selected survival predictive methods (only available on the last fold)
#'  \item \code{pred.surv} for each horizon time(s), containing the individual probabilities for the selected survival predictive methods
#'  \item \code{AUC} list of horizon time(s) containing AUC for each fold for the selected survival predictive methods
#'  \item \code{BS} list of horizon time(s) containing BS for each fold for the selected survival predictive methods
#' }
#' \item{long.method}{method(s) used to modeling the biomarkers}
#' \item{surv.methods}{method(s) used to compute the individual survival prediction}
#' \item{models.name}{name of survival prediction methods}
#' \item{kfolds}{number of folds}
#'
#' @author Anthony Devaux (\email{anthony.devaux@u-bordeaux.fr}) (maintener), Robin Genuer and Cécile Proust-Lima
#'
#' @references
#' \insertAllCited{}
#' @examples
#'
#' \dontrun{
#'
#' data(pbc2)
#'
#' # Formula for the modeling of the biomarkers using splines
#' serBilir = list(model = list(fixed = serBilir ~ year,
#'                              random = ~ year,
#'                              subject = "id"),
#'                 deriv = list(fixed = ~ 1,
#'                              indFixed = 2,
#'                              random = ~ 1,
#'                              indRandom = 2))
#'
#' serChol = list(model = list(fixed = serChol ~ year + I(year^2),
#'                             random = ~ year + I(year^2),
#'                             subject = "id"),
#'                deriv = list(fixed = ~ 2*year,
#'                             indFixed = c(2,3),
#'                             random = ~ 2*year,
#'                             indRandom = c(2,3)))
#'
#' albumin = list(model = list(fixed = albumin ~ year,
#'                             random = ~ year,
#'                             subject = "id"),
#'                deriv = list(fixed = ~ 1,
#'                             indFixed = 2,
#'                             random = ~ 1,
#'                             indRandom = 2))
#'
#' alkaline = list(model = list(fixed = alkaline ~ year,
#'                              random = ~ year,
#'                              subject = "id"),
#'                 deriv = list(fixed = ~ 1,
#'                              indFixed = 2,
#'                              random = ~ 1,
#'                              indRandom = 2))
#'
#' SGOT = list(model = list(fixed = SGOT ~ year,
#'                          random = ~ year,
#'                          subject = "id"),
#'             deriv = list(fixed = ~ 1,
#'                          indFixed = 2,
#'                          random = ~ 1,
#'                          indRandom = 2))
#'
#' platelets = list(model = list(fixed = platelets ~ year + I(year^2),
#'                               random = ~ year + I(year^2),
#'                               subject = "id"),
#'                  deriv = list(fixed = ~ 2*year,
#'                               indFixed = c(2,3),
#'                               random = ~ 2*year,
#'                               indRandom = c(2,3)))
#'
#' prothrombin = list(model = list(fixed = prothrombin ~ year,
#'                                 random = ~ year,
#'                                 subject = "id"),
#'                    deriv = list(fixed = ~ 1,
#'                                 indFixed = 2,
#'                                 random = ~ 1,
#'                                 indRandom = 2))
#'
#' ascites = list(model = ascites ~ year + (1 + year|id),
#'                deriv = list(fixed = ~ 1,
#'                             indFixed = 2,
#'                             random = ~ 1,
#'                             indRandom = 2))
#'
#' hepatomegaly = list(model = hepatomegaly ~ year + (1 + year|id),
#'                     deriv = list(fixed = ~ 1,
#'                                  indFixed = 2,
#'                                  random = ~ 1,
#'                                  indRandom = 2))
#'
#' spiders = list(model = spiders ~ year + (1 + year|id),
#'                deriv = list(fixed = ~ 1,
#'                             indFixed = 2,
#'                             random = ~ 1,
#'                             indRandom = 2))
#'
#' edema2 = list(model = list(fixed = edema2 ~ year,
#'                            random = ~ year,
#'                            subject = "id"),
#'               deriv = list(fixed = ~ 1,
#'                            indFixed = 2,
#'                            random = ~ 1,
#'                            indRandom = 2))
#'
#' marker <- list(serBilir = serBilir, serChol = serChol, albumin = albumin,
#'                alkaline = alkaline, SGOT = SGOT, platelets = platelets,
#'                prothrombin = prothrombin, ascites = ascites,
#'                hepatomegaly = hepatomegaly, spiders = spiders,
#'                edema2 = edema2)
#'
#' # compute hdlandmark methodology
#' hdlandmark.res <- hdlandmark(data = pbc2, data.pred = pbc2, markers = marker,
#'                              tLMs = 4, tHors = 3,
#'                              subject = "id", time = "year", time.event = "years", event = "status2",
#'                              long.method = "GLMM", surv.covar = "baseline",
#'                              cox.submodels = "allVar",
#'                              coxnet.submodels = "lasso",
#'                              spls.submodels = "nosparse",
#'                              rsf.submodels = "default",
#'                              rsf.split = c("logrank"),
#'                              kfolds = 10)
#'
#' # get individual predictions for each method
#' hdlandmark.res$models[[`4`]]$pred.surv$`3`
#' }
#'
#' @export
hdlandmark <- function(data, data.pred = NULL, markers, tLMs, tHors,
                       subject, time, time.event, event,
                       long.method = c("combine", "GLMM", "MFPC"),
                       surv.covar = c("baseline","LOtLM"),
                       cox.submodels = c("autoVar","allVar"), coxnet.submodels = c("opt","lasso","ridge"),
                       penaFG.submodels = c("GCV","BIC"), spls.submodels = c("opt","nosparse","maxsparse"),
                       rsf.submodels = c("opt","noVS","default"), rsf.split = c("logrank", "bs.gradient"),
                       cause = 1, HW = NULL, summaries = c("RE","score","pred","slope","cumulative"),
                       kfolds = 10, seed = 1234, scaling = FALSE, SL.weights = NULL,
                       nodesize.grid = NULL, mtry.grid = NULL){

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
  if (length(unique(data[,event]))>2){ # competitive risks
    CR <- TRUE
  }else{
    CR <- FALSE
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
  if (length(surv.covar)>1){
    surv.covar <- surv.covar[1]
    if (!all(surv.covar%in%c("baseline","LOtLM"))){
      stop("Only baseline or LOtLM packages are allowed for surv.covar")
    }
  }
  if (!all(cox.submodels%in%c("autoVar","allVar"))){
    stop("Only autoVar or allVar are allowed for cox.submodels")
  }
  if (!all(coxnet.submodels%in%c("opt","lasso","ridge"))){
    stop("Only opt, lasso or ridge are allowed for coxnet.submodels")
  }
  if (CR){
    if (any(coxnet.submodels=="opt")){
      warning("opt coxnet.submodels is not available for competitive risks !")
      coxnet.submodels <- coxnet.submodels[-which(coxnet.submodels=="opt")]
    }
  }else{
    if (!is.null(penaFG.submodels)){
      warning("penaFG.submodels is only available for competitive risks !")
      penaFG.submodels <- NULL
    }
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
  if (CR){
    if (any(rsf.submodels=="ranger")){
      warning("ranger rsf.submodels is not available for competitive risks !")
      rsf.submodels <- rsf.submodels[-which(rsf.submodels=="ranger")]
    }
    if (any(rsf.split=="bs.gradient")){
      warning("bs.gradient splitting rule is not available for competitive risks !")
      rsf.split <- rsf.split[-which(rsf.split=="bs.gradient")]
    }
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

  if (!is.null(HW)){
    if (HW>=tLM){
      stop("HW value should be smaller than the landmark time !")
    }
  }else{
    HW <- min(c(data[,time], data.pred[,time]), na.rm = T)
  }

  if (!all(summaries%in%c("RE","score","pred","slope","cumulative"))){
    stop("Wrong summary chosen ! Actually, only the following summaries are available : RE, score, pred, slope and cumulative !")
  }

  surv.methods <- NULL

  if (length(cox.submodels)>0){
    surv.methods <- c(surv.methods, "cox")
  }
  if (length(coxnet.submodels)>0){
    surv.methods <- c(surv.methods, "penalized-cox")
  }
  if (length(penaFG.submodels)>0){
    surv.methods <- c(surv.methods, "penalized-FG")
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

    data.tLM <- data[which(data[,time]<=tLM&data[,time.event]>tLM),]

    if (!is.null(data.pred)){ # different data training and test

      data.pred.tLM <- data.pred[which(data.pred[,time]<=tLM&data.pred[,time.event]>tLM),]
      kfolds <- 1

    }else{ # same data estimation/prediction

      data.pred.tLM <- data.tLM

    }

    ids <- unique(data.tLM[,subject])
    n <- length(ids)
    set.seed(seed)
    fold <- sample(rep(1:kfolds, length.out = n))

    pred.surv <- AUC <- BS <- list()

    for (k in 1:kfolds){

      cat(paste0("Fold : ",k,"/",kfolds),"\n")

      ids.test <- ids[which(fold==k)]

      if (kfolds==1){

        ids.train <- ids.test

        data.k <- data.tLM
        data.pred.k <- data.pred.tLM

      }else{

        ids.train <- ids[which(fold!=k)]

        data.k <- data.tLM[which(data.tLM[,subject]%in%ids.train),]
        data.pred.k <- data.pred.tLM[which(data.pred.tLM[,subject]%in%ids.test),]

      }

      # estimation of summaries on training data and test data

      res.LMsum <- LMsummaries(data = data.k, data.pred = data.pred.k, markers = markers, tLM = tLM,
                               subject = subject, time = time, time.event = time.event, event = event,
                               long.method = long.method, surv.covar = surv.covar,
                               scaling = scaling, HW = HW, summaries = summaries)

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
                             penaFG.submodels = penaFG.submodels,
                             spls.submodels = spls.submodels, rsf.submodels = rsf.submodels,
                             rsf.split = rsf.split, cause = cause, CR = CR,
                             nodesize.grid = nodesize.grid, mtry.grid = mtry.grid)

        # survival model on test data
        res.LMpred <- LMpred(data.surv = data.surv.pred, model.surv = res.LMsurv$model.surv,
                             long.method = long.method, surv.methods = surv.methods,
                             tHor = tHor, cause = cause, CR = CR)

        if (!is.null(SL.weights)){

          SL.pred <- res.LMpred$pred.surv%*%SL.weights
          colnames(SL.pred) <- "superlearner"
          res.LMpred$pred.surv <- cbind(res.LMpred$pred.surv, SL.pred)

        }

        #res.LMassess <- LMassess(pred.surv = res.LMpred$pred.surv, data.surv = data.surv.pred, tHor = tHor)

        #AUC[[as.character(tHor)]] <- rbind(AUC[[as.character(tHor)]], res.LMassess$AUC)
        #BS[[as.character(tHor)]] <- rbind(BS[[as.character(tHor)]], res.LMassess$BS)

        pred.surv[[as.character(tHor)]] <- rbind(pred.surv[[as.character(tHor)]], res.LMpred$pred.surv)
        pred.surv[[as.character(tHor)]] <- pred.surv[[as.character(tHor)]][order(as.integer(rownames(pred.surv[[as.character(tHor)]]))), , drop = FALSE]

      }

    }

    models[[as.character(tLM)]] <- list(data.surv = res.LMsum$data.surv, data.surv.pred = res.LMsum$data.surv.pred,
                                        model.surv = res.LMsurv$model.surv, pred.surv = pred.surv,
                                        AUC = AUC, BS = BS)

  }

  cat("DONE !", "\n")

  resu <- list(tLMs = tLMs, tHors = tHors, models = models,
               long.method = long.method, surv.methods = surv.methods,
               models.name = colnames(models[[1]]$pred.surv[[1]]), kfolds = kfolds)

  class(resu) <- "hdlandmark"

  return(resu)

}
