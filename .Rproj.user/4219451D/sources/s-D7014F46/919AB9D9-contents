landmark <- function(lmm_objects, data, tLM, subject, time, derivForm_objects){

  data_landmark <- data[which(data[,time]<tLM),]

  # build survival data

  data_surv <- data.frame()

  for (ind_subject in unique(data_landmark[,subject])){
    temp_subject <- data_landmark[which(data_landmark[,subject]==ind_subject),]
    temp_subject <- temp_subject[which.min(temp_subject[,time]),]

    data_surv <- rbind(data_surv, temp_subject)
  }

  data_surv[,time] <- tLM

  marker_ind <- 1

  data_surv <- sapply(lmm_objects, FUN = function(lmm_objects){

    browser()

    marker_name <- as.character(lmm_objects$call$fixed)[2]

    predRE <- predRE_newdata(lmm_objects, data_landmark)

    predY <- predY_newdata(predRE, data_surv)
    data_surv[,marker_name] <- predY

    derivY <- derivY_newdata(predRE, data_surv, derivForm[[marker_ind]])
    data_surv[,ncol(data_surv) + 1] <- derivY
    colnames(data_surv)[ncol(data_surv)] <- paste(marker_name, "slope", sep = "_")

    marker_ind <- marker_ind + 1

    return(data_surv)

  })

  return(data_surv)
}
