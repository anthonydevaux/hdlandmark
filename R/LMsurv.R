LMsurv <- function(data.surv, subject, time.event, event, long.methods, surv.methods){

  model.list <- list()

  for (long.method in long.methods){

    for (surv.method in surv.methods){

      method.name <- paste(long.method, surv.method, sep = "-")

      if (surv.method=="cox"){

        model.list[[method.name]] <- LMsurv.cox(data.surv[[long.method]])

      }

      if (surv.method=="penalized-cox"){

        model.list[[method.name]] <- LMsurv.coxnet(data.surv[[long.method]])

      }

      if (surv.method=="sPLS"){

        model.list[[method.name]] <- LMsurv.spls(data.surv[[long.method]])

      }

      if (surv.method=="RSF"){

        model.list[[method.name]] <- LMsurv.rsf(data.surv[[long.method]])

      }

    }

  }

}
