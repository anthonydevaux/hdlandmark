LMsurv <- function(data.surv, time.event, event, long.methods, surv.methods){

  model.list <- list()

  for (long.method in long.methods){

    model.list[[long.method]] <- list()

    for (surv.method in surv.methods){

      if (surv.method=="cox"){

        model.list[[long.method]][[surv.method]] <- LMsurv.cox(data.surv[[long.method]])

      }

      if (surv.method=="penalized-cox"){

        model.list[[long.method]][[surv.method]] <- LMsurv.coxnet(data.surv[[long.method]])

      }

      if (surv.method=="sPLS"){

        model.list[[long.method]][[surv.method]] <- LMsurv.spls(data.surv[[long.method]])

      }

      if (surv.method=="rsf"){

        model.list[[long.method]][[surv.method]] <- LMsurv.rsf(data.surv[[long.method]])

      }

    }

  }

  return(list(model.surv = model.list))

}
