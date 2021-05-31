#' Title
#'
#' @param data.surv
#' @param penaFG.submodels
#' @param cause
#'
#' @return
#' @export
#'
#' @importFrom crrp crrp
#' @importFrom riskRegression FGR
#' @examples
LMsurv.penaFG <- function(data.surv, penaFG.submodels, cause = 1){

  model.penaFG <- list()

  data.surv.omit <- na.omit(data.surv[,!(names(data.surv) %in% "subject")])
  data.surv.X <- model.matrix( ~ ., data.surv.omit[,!(names(data.surv.omit) %in% c("time.event","event"))])[,-1]
  data.surv.Y <- data.surv.omit[,c("time.event","event")]

  #######################################################
  # lasso penalisation

  crrp.FG.fit <- crrp::crrp(time=data.surv.Y$time.event, fstatus=data.surv.Y$event, X=data.surv.X,
                            penalty="LASSO", failcode = cause, cencode = 0,
                            lambda = 10^seq(log10(0.05), log10(0.0002), length = 50), eps = 1E-6)

  if (any(penaFG.submodels %in% c("GCV"))){

    # find the tuning parameter with GCV for a lasso penalisation
    beta.GCV <- crrp.FG.fit$beta[, which.min(crrp.FG.fit$GCV)]
    beta.GCV<- as.matrix(beta.GCV)

    FGcoef.GCV <- Matrix(beta.GCV, sparse=TRUE)#transforme en sparse pour pouvoir recuperer les noms de colonnes qui nous int?resse
    FGactivVar.GCV <- FGcoef.GCV@Dimnames[[1]][FGcoef.GCV@i+1] # keep variables with non zero coefficients

    newformula.GCV <- reformulate(termlabels = FGactivVar.GCV,
                                  response = "Hist(time.event,event)")

    # with riskRegression to be able to use the FG prediction function of riskRegression
    riskReg.FG.GCV.fit <- riskRegression::FGR(formula=newformula.GCV, data=data.surv.omit, maxiter=0, init=FGcoef.GCV@x, cause=cause)

    model.penaFG[["GCV"]] <- riskReg.FG.GCV.fit

  }

  if (any(penaFG.submodels %in% c("BIC"))){

    # find the tuning parameter with BIC for a lasso penalisation
    beta.BIC <- crrp.FG.fit$beta[, which.min(crrp.FG.fit$BIC)]
    beta.BIC <-as.matrix(beta.BIC)

    FGcoef.BIC <- Matrix(beta.BIC, sparse=TRUE)#transforme en sparse pour pouvoir recuperer les noms de colonnes qui nous int?resse
    FGactivVar.BIC <- FGcoef.BIC@Dimnames[[1]][FGcoef.BIC@i+1] # keep variables with non zero coefficients

    newformula.BIC <- reformulate(termlabels = FGactivVar.BIC,
                                  response = "Hist(time.event,event)")
    #with riskRegression to be able to use the FG prediction function of riskRegression
    riskReg.FG.BIC.fit <- riskRegression::FGR(formula=newformula.BIC, data=data.surv.omit, maxiter=0, init=FGcoef.BIC@x, cause=cause)

    model.penaFG[["BIC"]] <- riskReg.FG.BIC.fit

  }

  return(list(model = model.penaFG, surv.name = "penalized-FG"))

}
