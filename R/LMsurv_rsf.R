#' Title
#'
#' @param data.surv
#' @param rsf.split
#' @param rsf.submodels
#' @param cause
#' @param CR
#'
#' @return
#' @export
#'
#' @import randomForestSRC ranger
#' @importFrom survival Surv
#'
#' @examples
LMsurv.rsf <- function(data.surv, rsf.split, rsf.submodels, cause = 1, CR = FALSE){

  model.rsf <- list()

  data.surv <- data.surv[,!(names(data.surv) %in% "subject")]

  best.param <- NULL

  for (splitrule in rsf.split){

    if (any(rsf.submodels %in% c("opt"))){ # rsf tuning

      best.err <- 1

      mtry.max <- ncol(data.surv) - 2 # (no count time.event and event variables)
      nodesize.max <- 20

      for (nodesize in seq(1,nodesize.max,5)){
        for (mtry in seq(1,mtry.max,5)){

          cat(paste0("Nodesize : ", nodesize, " and mtry : ", mtry), "\n")

          res.rsf <- rfsrc(Surv(time.event, event) ~ ., data.surv,
                           nodesize = nodesize,
                           mtry = mtry,
                           ntree = 500,
                           nsplit = 0,
                           splitrule = splitrule,
                           cause = cause,
                           bootstrap = "by.root", samptype = "swr")

          rf.err <- tail(res.rsf$err.rate, 1)

          if (rf.err < best.err){
            best.err <- rf.err
            best.param <- c("nodesize" = res.rsf$nodesize, "mtry" = res.rsf$mtry)
            best.rsf <- res.rsf
          }
          cat(paste0("Best error : ", round(best.err,4), " with nodesize = ",
                     best.param[1], " and mtry = ", best.param[2]), "\n")
        }
      }

      if (any(rsf.submodels %in% c("noVS"))){

        model.rsf[[paste0(splitrule, "-noVS", ifelse(CR, "-CR", ""))]] <- res.rsf

      }

      VI.obj <- vimp.rfsrc(res.rsf, importance = "permute")

      # VIMP > 0.005
      VI.var.selected <- names(VI.obj$importance[VI.obj$importance>0.005])

      # Top 10% variables
      # topten <- round(length(VI.obj$importance)/10)
      # VI.var.selected <- names(sort(VI.obj$importance, decreasing = T)[1:topten])

      best.rf.VI <- rfsrc(Surv(time.event, event) ~ ., data.surv[,c(VI.var.selected,"time.event","event")],
                          ntree = 1000,
                          nsplit = 0,
                          splitrule = splitrule,
                          cause = cause,
                          bootstrap = "by.root", samptype = "swr")

      model.rsf[[paste0(splitrule, "-opt", ifelse(CR, "-CR", ""))]] <- best.rf.VI

    }

    if (any(rsf.submodels %in% c("default"))){

      res.rsf <- rfsrc(Surv(time.event, event) ~ ., data.surv,
                       ntree = 1000,
                       nsplit = 0,
                       splitrule = splitrule,
                       cause = cause,
                       bootstrap = "by.root", samptype = "swr")

      model.rsf[[paste0(splitrule, "-default", ifelse(CR, "-CR", ""))]] <- res.rsf

    }

  }

  ############################################
  ### Ranger forests

  if (any(rsf.submodels %in% c("ranger"))){ # rsf tuning

    data.surv.omit <- na.omit(data.surv)

    res.rsf <- ranger(Surv(time.event, event)~., data = data.surv.omit,
                      num.trees = 1000)
    model.rsf[["ranger"]] <- res.rsf

  }

  if (!is.null(best.param)){
    return(list(model = model.rsf, mtry.opt = best.param[2], nodesize.opt = best.param[1],
                surv.name = "rsf"))
  }else{
    return(list(model = model.rsf, mtry.opt = NULL, nodesize.opt = NULL,
                surv.name = "rsf"))
  }



}
