LMsurv <- function(data.surv, surv.methods, cox.submodels, coxnet.submodels,
                   spls.submodels, rsf.submodels, rsf.split){

  model.surv <- list()

  for (surv.method in surv.methods){

    if (surv.method=="cox"){

      model.surv[[surv.method]] <- LMsurv.cox(data.surv = data.surv, cox.submodels = cox.submodels)

    }

    if (surv.method=="penalized-cox"){

      model.surv[[surv.method]] <- LMsurv.coxnet(data.surv = data.surv, coxnet.submodels = coxnet.submodels)

    }

    if (surv.method=="spls"){

      model.surv[[surv.method]] <- LMsurv.spls(data.surv = data.surv, spls.submodels = spls.submodels)

    }

    if (surv.method=="rsf"){

      model.surv[[surv.method]] <- LMsurv.rsf(data.surv = data.surv, rsf.split = rsf.split,
                                              rsf.submodels = rsf.submodels)

    }

  }

  return(list(model.surv = model.surv))

}
