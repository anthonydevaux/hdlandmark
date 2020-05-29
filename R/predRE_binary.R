#' Title
#'
#' @param model
#' @param data
#'
#' @return
#' @export
#'
#' @import lme4
#' @importFrom Matrix bdiag
#'
#' @examples
predRE.binary <- function(model, data){

  model.formula <- formula(model)

  var.group <- names(ranef(model))

  X.formula <- as.formula(nobars(model.formula)[-2])
  Z.formula <- as.formula(paste("~", as.character(findbars(model.formula)[[1]])[2]))
  #Y.formula <- as.formula(paste("~", -1, "+", as.character(model.formula)[2]))

  # parametres
  beta <- fixef(model) # fixed effect
  B <- as.matrix(bdiag(VarCorr(model))) # random effects matrice var-covar

  # creation matrix estimation random effects

  bi.mat <- matrix(NA, nrow = length(unique(data[,var.group])), ncol = ncol(B),
                   dimnames = list(unique(data[,var.group]), colnames(B)))

  # boucle pour chaque sujet

  for (ind_subject in unique(data[,var.group])){

    # historique des mesures repetees le sujet selectionne
    data.subject <- na.omit(data[which(data[,var.group]==ind_subject), all.vars(terms(model.formula))])

    if (nrow(data.subject)>0) {

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

    }else{

      bi.mat[which(rownames(bi.mat)==ind_subject),] <- NA

    }


  }

  return(list(b_i = bi.mat,
              beta = beta,
              call = model.formula,
              group = var.group))
}
