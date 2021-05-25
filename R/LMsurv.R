#' Title
#'
#' @param data.surv
#' @param surv.methods
#' @param cox.submodels
#' @param coxnet.submodels
#' @param spls.submodels
#' @param rsf.submodels
#' @param rsf.split
#' @param cause
#' @param CR
#'
#' @return
#' @export
#'
#' @examples
LMsurv <- function(data.surv, surv.methods, cox.submodels, coxnet.submodels,
                   spls.submodels, rsf.submodels, rsf.split, cause = 1, CR = FALSE){

  cat("Survival estimation on data train","\n")
  cat("----------------------------------","\n")

  model.surv <- list()

  for (surv.method in surv.methods){

    if (surv.method=="cox"){

      cat("cox...")

      model.surv[[surv.method]] <- LMsurv.cox(data.surv = data.surv, cox.submodels = cox.submodels)

    }

    if (surv.method=="FG"){

      cat("Fine&Gray...")

    }

    if (surv.method=="penalized-cox"){

      cat("Penalized-Cox...")

      if (!CR){

        model.surv[[surv.method]] <- LMsurv.coxnet(data.surv = data.surv, coxnet.submodels = coxnet.submodels)

      }else{

        model.surv[[surv.method]] <- LMsurv.coxnet.CR(data.surv = data.surv, coxnet.submodels = coxnet.submodels,
                                                      cause = cause)

      }

    }

    if (surv.method=="spls"){

      cat("sPLS...")

      model.surv[[surv.method]] <- LMsurv.spls(data.surv = data.surv, spls.submodels = spls.submodels)

    }

    if (surv.method=="rsf"){

      cat("RSF...")

      model.surv[[surv.method]] <- LMsurv.rsf(data.surv = data.surv, rsf.split = rsf.split,
                                              rsf.submodels = rsf.submodels, cause = cause, CR = CR)

    }

  }

  cat("--", "\n")

  return(list(model.surv = model.surv))

}
