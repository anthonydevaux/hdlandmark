#' Title
#'
#' @param model
#' @param data
#'
#' @return
#' @export
#'
#' @importFrom lme4 nobars findbars
#'
#' @examples
#'
predRE.binary_simu <- function(model, data){

  model.formula <- model$model

  var.group <- "id"

  X.formula <- as.formula(nobars(model.formula)[-2])
  Z.formula <- as.formula(paste("~", as.character(findbars(model.formula)[[1]])[2]))
  #Y.formula <- as.formula(paste("~", -1, "+", as.character(model.formula)[2]))

  # parametres
  beta <- model$params$beta # fixed effect

  # random effects matrice var-covar

  len.sd.re <- length(model$params$sd.re)

  if (len.sd.re==1){

    B <- as.matrix(model$params$sd.re)

  }else if (len.sd.re==2){

    covar <- model$params$sd.re[1] * model$params$sd.re[2] * model$params$cor.re

    B <- matrix(c(model$params$sd.re[1]^2, covar,
                  covar, model$params$sd.re[2]^2), ncol = 2)

  }else if (len.sd.re==3){

    covar01 <- model$params$sd.re[1] * model$params$sd.re[2] * model$params$cor.re[1]
    covar02 <- model$params$sd.re[1] * model$params$sd.re[3] * model$params$cor.re[2]
    covar12 <- model$params$sd.re[2] * model$params$sd.re[3] * model$params$cor.re[3]

    B <- matrix(c(model$params$sd.re[1]^2, covar01, covar02,
                  covar01, model$params$sd.re[2]^2, covar12,
                  covar02, covar12, model$params$sd.re[3]^2), ncol = 3)

  }

  # creation matrix estimation random effects

  bi.mat <- matrix(NA, nrow = length(unique(data[,var.group])), ncol = ncol(B),
                   dimnames = list(unique(data[,var.group]), colnames(B)))

  # boucle pour chaque sujet

  for (ind_subject in unique(data[,var.group])){

    # historique des mesures repetees le sujet selectionne
    data.subject <- na.omit(data[which(data[,var.group]==ind_subject), all.vars(terms(model.formula))])

    X.mat <- model.matrix(X.formula, data.subject)
    Z.mat <- model.matrix(Z.formula, data.subject)
    Y <- data.subject[,as.character(model.formula)[2]]
    levels(Y) <- c(0,1)

    bi.function <- function(bi, X.mat, beta, Z.mat, Y, B){

      num <- exp(X.mat%*%beta+Z.mat%*%bi)^ifelse(Y==1,1,0)
      den <- 1+exp(X.mat%*%beta+Z.mat%*%bi)
      f1 <- prod(num/den) # f(Yi|bi)

      q <- nrow(B)
      f2 <- ((2*pi)^(-q/2))*(det(B)^(-1/2))*exp(-0.5*bi%*%solve(B)%*%bi) # f(bi) loi normale multivariee

      return(f1*f2)
    }

    # optim

    neg.bi.function <- function(bi, X.mat, beta, Z.mat, Y, B){-bi.function(bi, X.mat, beta, Z.mat, Y, B)}

    maxi.opt <- optim(rep(0.1,ncol(B)), fn = neg.bi.function,
                      X.mat = X.mat, beta = beta, Z.mat = Z.mat, Y = Y, B = B)

    bi.mat[which(rownames(bi.mat)==ind_subject),] <- maxi.opt$par

    # marquart

    # maxi.opt <- marqLevAlg(m = ncol(B), fn = bi.function,
    #                        X.mat = X.mat, beta = beta, Z.mat = Z.mat, Y = Y, B = B,
    #                        minimize = FALSE)
    #
    # bi.mat[which(rownames(bi.mat)==ind_subject),] <- maxi.opt$b
  }

  return(list(b_i = bi.mat,
              beta = beta,
              call = model.formula,
              group = var.group))
}
