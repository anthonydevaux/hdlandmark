#' Title
#'
#' @param data
#' @param data.surv
#' @param markers
#' @param tLM
#' @param subject
#' @param time
#' @param time.interval
#' @param pve
#' @param nbasis
#'
#' @return
#' @export
#'
#' @importFrom arules discretize
#'
#' @examples
LMsum.mfpc <- function(data, data.surv, data.pred, data.surv.pred, markers, tLM, subject, time,
                       time.interval = 0.1, pve = 0.95, nbasis = 3){

  markers.names <- names(markers)

  # only on numeric variables

  markers.names <- markers.names[unlist(lapply(markers.names, function(x) return(!is.factor(data[,x]))))]
  P <- length(markers.names)

  ###### Training data ######
  # discretization of time data

  time.discret <- discretize(data[,time], method = "fixed",
                             breaks = seq(0, tLM, time.interval), labels = FALSE)
  u.time <- attr(time.discret, "discretized:breaks")
  data[,time] <- u.time[time.discret]
  #u.time <- sort(unique(data[,time]))
  n.time <- length(u.time)
  u.time.scale <- u.time/max(u.time) # scale the time domain to [0,1]
  ids <- unique(data[,subject])

  # build markers value array

  markers.value <- array(NA, c(length(ids), n.time, P))

  for(i in 1:length(ids)){

    visits <- which(u.time %in% (data[which(data[,subject] == ids[i]), time]))

    for (marker in markers.names){

      markers.value[i,visits, which(markers.names%in%marker)] <- data[which(data[,subject] == ids[i]), marker]

    }
  }

  ####### Test data ######
  # discretization of time data

  time.discret <- discretize(data.pred[,time], method = "fixed",
                             breaks = u.time, labels = FALSE)
  u.time.pred <- attr(time.discret, "discretized:breaks")
  data.pred[,time] <- u.time.pred[time.discret]

  ids.pred <- unique(data.pred[,subject])

  # build markers value array

  markers.value.pred <- array(NA, c(length(ids.pred), n.time, P))

  for(i in 1:length(ids.pred)){

    visits <- which(u.time %in% (data.pred[which(data.pred[,subject] == ids.pred[i]), time]))

    for (marker in markers.names){

      markers.value.pred[i,visits, which(markers.names%in%marker)] <- data.pred[which(data.pred[,subject] == ids.pred[i]), marker]

    }
  }

  Xi <- Xi.pred <- L <- sigma2 <- NULL
  phi <- mean.fct <- fit <- list()

  ########## PACE #######
  # estimate parameters from data for each marker + prediction on test data

  for(p in 1:P){

    fct.data <- funData(u.time.scale, markers.value[,,p])
    fct.data.pred <- funData(u.time.scale, markers.value.pred[,,p])
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

  colnames(Xi) <- colnames(Xi.pred) <- unlist(lapply(seq_along(markers.names),
                                function(x) return(paste0(markers.names[x], ".fpc.score" , seq(L[x])))))

  data.surv <- cbind(data.surv, Xi)
  data.surv.pred <- cbind(data.surv.pred, Xi.pred)

  return(list(data.surv = data.surv, data.surv.pred = data.surv.pred, Xi = Xi,
              Xi.pred = Xi.pred, L = L, phi = phi, mean.fct = mean.fct,
              fit = fit, sigma2 = sigma2))

}
