#' Title
#'
#' @param data
#' @param id
#' @param time
#' @param u.time
#' @param markers
#' @param pve
#' @param nbasis
#' @param pred.data
#' @param plot
#'
#' @return
#' @export
#'
#' @import funData
#'
#' @examples
LM.MFPC <- function(data, id, time, u.time, markers, pred.data = NULL, pve = 0.95, nbasis = 10, plot = FALSE){

  ids <- unique(data[,id])
  P <- length(markers)
  n.time <- length(u.time)
  u.time.scale <- u.time/max(u.time)

  multivar <- array(NA, c(length(ids), n.time, P))

  for(i in 1:length(ids)){

    visits <- which(u.time %in% (data[which(data[,id] == ids[i]), time]))

    for (marker in markers){

      multivar[i,visits, which(markers%in%marker)] <- data[which(data[,id] == ids[i]), marker]

    }
  }

  Xi <- Xi.pred <- L <- sigma2 <- NULL
  phi <- mean.fct <- fit <- list()

  if (is.null(pred.data)){

    # estimate parameters from data

    for(p in 1:P){

      fct.data <- funData(u.time.scale, multivar[,,p])
      res.PACE <- PACE(fct.data, pve = pve, nbasis = nbasis)
      Xi = cbind(Xi, res.PACE$scores) # FPC scores
      L = c(L, dim(res.PACE$scores)[2])
      phi[[p]] = t(res.PACE$functions@X) # FPC eigenfunctions
      mean.fct[[p]] = res.PACE$mu@X # estimated mean functions
      sigma2 <- c(sigma2, res.PACE$sigma2) # error measurement

    }

    # surv data with scores

    data.surv <- data.frame()

    for (ind_subject in ids){
      temp_subject <- data[which(data[,id]==ind_subject),]
      temp_subject <- temp_subject[which.min(temp_subject[,time]),]

      data.surv <- rbind(data.surv, temp_subject)
    }

    colnames(Xi) <- paste("xi.scores", seq(ncol(Xi)), sep = ".")
    data.surv <- cbind(data.surv, Xi)

  }else{

    ids.pred <- unique(pred.data[,id])

    multivar.pred <- array(NA, c(length(ids.pred), n.time, P))

    for(i in 1:length(ids.pred)){

      visits <- which(u.time %in% (pred.data[which(pred.data[,id] == ids.pred[i]), time]))

      for (marker in markers){

        multivar.pred[i,visits, which(markers%in%marker)] <- pred.data[which(pred.data[,id] == ids.pred[i]), marker]

      }
    }

    for(p in 1:P){

      fct.data <- funData(u.time.scale, multivar[,,p])
      fct.data.pred <- funData(u.time.scale, multivar.pred[,,p])
      res.PACE <- PACE(fct.data, pve = pve, nbasis = nbasis)
      res.PACE.pred <- PACE(fct.data, fct.data.pred, pve = pve, nbasis = nbasis)
      Xi = cbind(Xi, res.PACE$scores) # FPC scores
      Xi.pred <- cbind(Xi.pred, res.PACE.pred$scores) # dynamic FPC scores for test subjects
      L = c(L, dim(res.PACE$scores)[2])
      phi[[p]] = t(res.PACE$functions@X) # FPC eigenfunctions
      mean.fct[[p]] = res.PACE$mu@X # estimated mean functions
      sigma2 <- c(sigma2, res.PACE$sigma2) # error measurement
      fit[[p]] <- res.PACE.pred$fit # fit from Karhunen-Loeve expansion

    }

    # surv data with scores

    data.surv <- data.frame()

    for (ind_subject in ids){
      temp_subject <- data[which(data[,id]==ind_subject),]
      temp_subject <- temp_subject[which.min(temp_subject[,time]),]

      data.surv <- rbind(data.surv, temp_subject)
    }

    colnames(Xi) <- paste("xi.scores", seq(ncol(Xi)), sep = ".")
    data.surv <- cbind(data.surv, Xi)

    data.surv.pred <- data.frame()

    for (ind_subject in ids.pred){
      temp_subject <- pred.data[which(pred.data[,id]==ind_subject),]
      temp_subject <- temp_subject[which.min(temp_subject[,time]),]

      data.surv.pred <- rbind(data.surv.pred, temp_subject)
    }

    colnames(Xi.pred) <- paste("xi.scores", seq(ncol(Xi.pred)), sep = ".")
    data.surv.pred <- cbind(data.surv.pred, Xi.pred)

  }

  # plot mean function and eigenfunctions for each marker

  if (plot){

    for (p in 1:P){

      print(plot(u.time, mean.fct[[p]][1,], type = "l",
                 main = paste(markers[p], "/", "mean function")))

      for (ind.dim in 1:ncol(phi[[p]])){

        print(plot(u.time, phi[[p]][,ind.dim], type = "l",
                   main = paste(markers[p], "/", "eigenfunction", ind.dim)))

      }
    }

  }

  return(list(data.surv = data.surv, Xi = Xi,
              data.surv.pred = data.surv.pred,
              Xi.pred = Xi.pred, L = L, phi = phi, mean.fct = mean.fct,
              fit = fit, sigma2 = sigma2))

}
