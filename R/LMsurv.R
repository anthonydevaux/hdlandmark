LMsurv <- function(data.surv, time.event, event, long.method, surv.methods,
                   cox.autoVar, cox.allVar,
                   coxnet.opt, coxnet.lasso, coxnet.ridge,
                   spls.opt, spls.nosparse, spls.maxsparse,
                   rsf.split, rsf.opt, rsf.noVS, rsf.default){

  model.surv <- list()

  for (surv.method in surv.methods){

    if (surv.method=="cox"){

      model.surv[[surv.method]] <- LMsurv.cox(data.surv = data.surv,
                                              cox.autoVar = cox.autoVar,
                                              cox.allVar = cox.allVar)

    }

    if (surv.method=="penalized-cox"){

      model.surv[[surv.method]] <- LMsurv.coxnet(data.surv = data.surv,
                                                 coxnet.opt = coxnet.opt,
                                                 coxnet.lasso = coxnet.lasso,
                                                 coxnet.ridge = coxnet.ridge)

    }

    if (surv.method=="sPLS"){

      model.surv[[surv.method]] <- LMsurv.spls(data.surv = data.surv,
                                               spls.opt = spls.opt,
                                               spls.nosparse = spls.nosparse,
                                               spls.maxsparse = spls.maxsparse)

    }

    if (surv.method=="rsf"){

      model.surv[[surv.method]] <- LMsurv.rsf(data.surv = data.surv,
                                              rsf.split = rsf.split,
                                              rsf.opt = rsf.opt,
                                              rsf.noVS = rsf.noVS,
                                              rsf.default = rsf.default)

    }

  }

  return(list(model.surv = model.surv))

}
