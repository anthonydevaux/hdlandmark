#' Title
#'
#' @param data
#' @param data.pred
#' @param markers
#' @param tLM
#' @param subject
#' @param time
#' @param time.event
#' @param event
#' @param long.method
#' @param surv.covar
#' @param scaling
#' @param HW
#' @param summaries
#'
#' @return
#' @export
#'
#' @examples
LMsummaries <- function(data, data.pred, markers, tLM,
                        subject, time, time.event, event,
                        long.method, surv.covar, scaling, HW, summaries){

  ########### Survival data at landmark time ###########

  # Training data

  data.surv <- data.frame()

  if (surv.covar=="baseline"){

    for (ind_subject in unique(data[,subject])){
      temp_subject <- data[which(data[,subject]==ind_subject),]
      temp_subject <- temp_subject[which.min(temp_subject[,time]),]

      data.surv <- rbind(data.surv, temp_subject)
    }

  }else if (surv.covar=="LOtLM"){

    for (ind_subject in unique(data[,subject])){
      temp_subject <- data[which(data[,subject]==ind_subject),]
      data_subject <- temp_subject[which.max(temp_subject[,time]),]

      var.na <- is.na(data_subject[,!(names(data_subject)%in%c(names(markers),subject,time,time.event,event))])

      if (any(var.na)){

        varnames.na <- colnames(var.na)[which(var.na==TRUE)]

        for (var in varnames.na){

          data_subject[,var] <- tail(na.omit(temp_subject[,var]),1)

        }

      }

      data.surv <- rbind(data.surv, data_subject)
    }

  }

  data.surv[,time.event] <- data.surv[,time.event] - tLM # scaling time data
  colnames(data.surv)[which(colnames(data.surv)==time.event)] <- "time.event"
  colnames(data.surv)[which(colnames(data.surv)==event)] <- "event"

  # order column for scaling
  data.surv <- data.surv[,c(which(colnames(data.surv)!="event"),which(colnames(data.surv)=="event"))]

  # Prediction data

  data.surv.pred <- data.frame()

  if (surv.covar=="baseline"){

    for (ind_subject in unique(data.pred[,subject])){
      temp_subject <- data.pred[which(data.pred[,subject]==ind_subject),]
      temp_subject <- temp_subject[which.min(temp_subject[,time]),]

      data.surv.pred <- rbind(data.surv.pred, temp_subject)
    }

  }else if (surv.covar=="LOtLM"){

    for (ind_subject in unique(data.pred[,subject])){
      temp_subject <- data.pred[which(data.pred[,subject]==ind_subject),]
      data_subject <- temp_subject[which.max(temp_subject[,time]),]

      var.na <- is.na(data_subject[,!(names(data_subject)%in%c(names(markers),subject,time,time.event,event))])

      if (any(var.na)){

        varnames.na <- colnames(var.na)[which(var.na==TRUE)]

        for (var in varnames.na){

          data_subject[,var] <- tail(na.omit(temp_subject[,var]),1)

        }

      }

      data.surv.pred <- rbind(data.surv.pred, data_subject)

    }

  }

  data.surv.pred[,time.event] <- data.surv.pred[,time.event] - tLM # scaling time data
  colnames(data.surv.pred)[which(colnames(data.surv.pred)==time.event)] <- "time.event"
  colnames(data.surv.pred)[which(colnames(data.surv.pred)==event)] <- "event"

  # order column for scaling
  data.surv.pred <- data.surv.pred[,c(which(colnames(data.surv.pred)!="event"),which(colnames(data.surv.pred)=="event"))]

  ############# Summaries ############

  # GLMM

  if (long.method == "GLMM"){

    # estimation

    res.glmm <- LMsum.glmm(data = data, data.surv = data.surv, markers = markers,
                           tLM = tLM, subject = subject, time = time,
                           summaries = summaries, HW = HW)
    data.surv.long <- res.glmm$data.surv[,!(names(res.glmm$data.surv) %in%
                                              c(time, names(markers)))]

    colnames(data.surv.long)[which(colnames(data.surv.long)==subject)] <- "subject"

    rownames(data.surv.long) <- data.surv.long$subject

    # prediction

    res.glmm.pred <- LMsum.glmm(data = data.pred, data.surv = data.surv.pred, markers = res.glmm$markers.model,
                                tLM = tLM, subject = subject, time = time,
                                summaries = summaries, HW = HW)
    data.surv.long.pred <- res.glmm.pred$data.surv[,!(names(res.glmm.pred$data.surv) %in%
                                                        c(time, names(markers)))]

    colnames(data.surv.long.pred)[which(colnames(data.surv.long.pred)==subject)] <- "subject"

    rownames(data.surv.long.pred) <- data.surv.long.pred$subject

  }

  # MFPC

  if (long.method == "MFPC"){

    # estimation + prediction

    res.mfpc <- LMsum.mfpc(data, data.surv, data.pred, data.surv.pred, markers, tLM, subject, time,
                           time.interval = 0.05, pve = 0.99, nbasis = 3)

    data.surv.long <- res.mfpc$data.surv[,!(names(res.mfpc$data.surv) %in%
                                              c(time, names(markers)))]
    data.surv.long.pred <- res.mfpc$data.surv.pred[,!(names(res.mfpc$data.surv.pred) %in%
                                                        c(time, names(markers)))]

    colnames(data.surv.long)[which(colnames(data.surv.long)==subject)] <- "subject"

    rownames(data.surv.long) <- data.surv.long$subject

    colnames(data.surv.long.pred)[which(colnames(data.surv.long.pred)==subject)] <- "subject"

    rownames(data.surv.long.pred) <- data.surv.long.pred$subject

  }

  # combine

  if (long.method == "combine"){

    # estimation

    res.glmm <- LMsum.glmm(data = data, data.surv = data.surv, markers = markers,
                           tLM = tLM, subject = subject, time = time,
                           summaries = summaries, HW = HW)
    data.surv.glmm <- res.glmm$data.surv

    res.mfpc  <- LMsum.mfpc(data, data.surv, data.pred, data.surv.pred, markers, tLM, subject, time,
                            time.interval = 0.05, pve = 0.99, nbasis = 3)
    data.surv.mfpc <- res.mfpc$data.surv

    data.surv.combine <- merge(data.surv.glmm, data.surv.mfpc, sort = FALSE)

    data.surv.long <- data.surv.combine[,!(names(data.surv.combine) %in% c(time, names(markers)))]

    colnames(data.surv.long)[which(colnames(data.surv.long)==subject)] <- "subject"

    rownames(data.surv.long) <- data.surv.long$subject

    # prediction

    res.glmm.pred <- LMsum.glmm(data = data.pred, data.surv = data.surv.pred, markers = res.glmm$markers.model,
                                tLM = tLM, subject = subject, time = time,
                                summaries = summaries, HW = HW)
    data.surv.pred.glmm <- res.glmm.pred$data.surv
    data.surv.pred.mfpc <- res.mfpc$data.surv.pred

    data.surv.combine.pred <- merge(data.surv.pred.glmm, data.surv.pred.mfpc, sort = FALSE)

    data.surv.long.pred <- data.surv.combine.pred[,!(names(data.surv.combine.pred) %in%
                                                       c(time, names(markers)))]

    colnames(data.surv.long.pred)[which(colnames(data.surv.long.pred)==subject)] <- "subject"

    rownames(data.surv.long.pred) <- data.surv.long.pred$subject

  }

  if (scaling){

    data.surv.long[,(which(colnames(data.surv.long)=="event")+1):ncol(data.surv.long)] <-
      scale(data.surv.long[,(which(colnames(data.surv.long)=="event")+1):ncol(data.surv.long)])

    data.surv.long.pred[,(which(colnames(data.surv.long.pred)=="event")+1):ncol(data.surv.long.pred)] <-
      scale(data.surv.long.pred[,(which(colnames(data.surv.long.pred)=="event")+1):ncol(data.surv.long.pred)])

  }

  return(list(data.surv = data.surv.long, data.surv.pred = data.surv.long.pred))

}
