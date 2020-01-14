#' Function to predict the random effects on data using hlme output
#'
#' @param model A \code{hlme} object
#' @param data A dataframe containing longitudinal data
#'
#' @return A list containing the random effect for each subject, the fixed coefficients
#' from \code{hlme} and the formula used to compute the \code{hlme} mixed model
#'
#' @importFrom lme4 ranef fixef VarCorr
#'
#' @examples
#'
predRE <- function(model, data){

  model.formula <- formula(model)

  subject <- names(ranef(model))

  row_subject <- rownames(data)

  beta <- fixef(model) # fixed effect
  B <- as.matrix(bdiag(VarCorr(model))) # random effects matrice var-covar
  se <- sigma(model)^2 # residual variance error

  # random design matrix

  Z.formula <- as.formula(paste("~", as.character(findbars(model.formula)[[1]])[2]))
  Z <- model.matrix(Z.formula, data)

  row_subject <- intersect(row_subject, rownames(Z))

  # fixed design matrix

  X.formula <- as.formula(nobars(model.formula)[-2])
  X <- model.matrix(X.formula, data)

  row_subject <- intersect(row_subject, rownames(X))

  # outcome
  Y <- na.omit(data[,as.character(model.formula)[2], drop = FALSE])

  row_subject <- intersect(row_subject, rownames(Y))

  data <- data[row_subject,]

  predRE <- matrix(NA, nrow = length(unique(data[,subject])), ncol = ncol(B),
                   dimnames = list(unique(data[,subject]), colnames(Z)))
  predRE_row <- 1

  for (ind_subject in unique(data[,subject])){

    ind <- rownames(data[which(data[, subject] == ind_subject),])

    Lsubject <- nrow(data[ind,])

    Z_i <- matrix(Z[ind,], nrow = Lsubject)
    V_i <- Z_i%*%B%*%t(Z_i) + se*diag(Lsubject)
    Y_i <- Y[ind,]
    X_i <- X[ind,]
    b_i <- B%*%t(Z_i)%*%solve(V_i)%*%(Y_i-X_i%*%beta)

    predRE[predRE_row,] <- b_i

    predRE_row <- predRE_row + 1
  }

  return(list(b_i = predRE,
              beta = beta,
              call = model.formula,
              group = subject))

}
