LMsummaries <- function(data, data.pred, markers, tLM,
                        subject, time, time.event, event,
                        long.methods){

  ########### Survival data at landmark time ###########

  data.surv.long <- list()
  data.surv.long.pred <- list()

  # Training data

  data.surv <- data.frame()

  for (ind_subject in unique(data[,subject])){
    temp_subject <- data[which(data[,subject]==ind_subject),]
    temp_subject <- temp_subject[which.min(temp_subject[,time]),]

    data.surv <- rbind(data.surv, temp_subject)
  }

  data.surv[,time.event] <- data.surv[,time.event] - tLM # scaling time data
  colnames(data.surv)[which(colnames(data.surv)%in%c(subject, time.event, event))] <-
    c("subject","time.event","event")

  # Prediction data

  data.surv.pred <- data.frame()

  for (ind_subject in unique(data.pred[,subject])){
    temp_subject <- data.pred[which(data.pred[,subject]==ind_subject),]
    temp_subject <- temp_subject[which.min(temp_subject[,time]),]

    data.surv.pred <- rbind(data.surv.pred, temp_subject)
  }

  data.surv.pred[,time.event] <- data.surv.pred[,time.event] - tLM # scaling time data
  colnames(data.surv.pred)[which(colnames(data.surv.pred)%in%c(subject, time.event, event))] <-
    c("subject","time.event","event")

  ############# Summaries ############

  # GLMM

  if (any(long.methods == "GLMM")){

    # estimation

    res.glmm <- LMsum.glmm(data, data.surv, markers, tLM, subject, time, threshold = NULL)
    data.surv.long[["GLMM"]] <- res.glmm$data.surv[,!(names(res.glmm$data.surv) %in%
                                                        c(time, names(markers)))]

    # prediction

    res.glmm.pred <- LMsum.glmm(data.pred, data.surv.pred, res.glmm$markers.model,
                                tLM, subject, time, threshold = NULL)
    data.surv.long.pred[["GLMM"]] <- res.glmm.pred$data.surv[,!(names(res.glmm.pred$data.surv) %in%
                                                                  c(time, names(markers)))]

  }

  # MFPC

  if (any(long.methods == "MFPC")){

    # estimation + prediction

    res.mfpc  <- LMsum.mfpc(data, data.surv, data.pred, data.surv.pred, markers, tLM, subject, time,
                            time.interval = 0.1, pve = 0.99, nbasis = 3)

    data.surv.long[["MFPC"]] <- res.mfpc$data.surv[,!(names(res.mfpc$data.surv) %in%
                                                        c(time, names(markers)))]
    data.surv.long.pred[["MFPC"]] <- res.mfpc$data.surv.pred[,!(names(res.mfpc$data.surv.pred) %in%
                                                                  c(time, names(markers)))]

  }

  # combine

  if (any(long.methods == "combine")){

    # estimation

    res.glmm <- LMsum.glmm(data, data.surv, markers, tLM, subject, time, threshold = NULL)
    data.surv.glmm <- res.glmm$data.surv

    res.mfpc  <- LMsum.mfpc(data, data.surv, data.pred, data.surv.pred, markers, tLM, subject, time,
                            time.interval = 0.1, pve = 0.99, nbasis = 3)
    data.surv.mfpc <- res.mfpc$data.surv

    data.surv.combine <- merge(data.surv.glmm, data.surv.mfpc, sort = FALSE)

    data.surv.long[["combine"]] <- data.surv.combine[,!(names(data.surv.combine) %in% c(time, names(markers)))]

    # prediction

    res.glmm.pred <- LMsum.glmm(data.pred, data.surv.pred, res.glmm$markers.model,
                                tLM, subject, time, threshold = NULL)
    data.surv.pred.glmm <- res.glmm.pred$data.surv
    data.surv.pred.mfpc <- res.mfpc$data.surv.pred

    data.surv.combine.pred <- merge(data.surv.pred.glmm, data.surv.pred.mfpc, sort = FALSE)

    data.surv.long.pred[["combine"]] <- data.surv.combine.pred[,!(names(data.surv.combine.pred) %in%
                                                                    c(time, names(markers)))]

  }

  return(list(data.surv = data.surv.long, data.surv.pred = data.surv.long.pred))

}
