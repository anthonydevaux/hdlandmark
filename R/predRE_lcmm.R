#' Function to predict the random effects on data using lcmm output
#'
#' @param model A \code{lcmm} object
#' @param formul
#' @param data A dataframe containing longitudinal data
#'
#' @return A list containing the random effect for each subject, the fixed coefficients
#' from \code{lcmm} and the formula used to compute the \code{lcmm} mixed model
#'
#' @export
#' @import lcmm
#' @importFrom Matrix bdiag
#'
#' @examples
#'
predRE.lcmm <- function(model, formul, data){

  subject <- formul$subject

  row_subject <- rownames(data)

  beta <- c(0, lcmm::fixef(model)[[2]]) # 0 (intercept not estimated) + fixed effect

  # Variance-covariance matrix of the random-effects

  B <- matrix(0, ncol = sum(model$idea0), nrow = sum(model$idea0))
  colnames(B) <- model$Xnames[model$idea0==1]
  rownames(B) <- model$Xnames[model$idea0==1]
  B[upper.tri(B,diag=TRUE)] <- model$best[(model$N[1]+model$N[2]+1):(model$N[1]+model$N[2]+model$N[3])]
  B <- t(B)
  B[upper.tri(B,diag=TRUE)] <- model$best[(model$N[1]+model$N[2]+1):(model$N[1]+model$N[2]+model$N[3])]

  se <- 1 # residual variance error (not estimated)

  # random design matrix

  Z.formula <- formul$random
  Z <- model.matrix(Z.formula, data)

  row_subject <- intersect(row_subject, rownames(Z))

  # fixed design matrix

  X.formula <- formul$fixed[-2]
  X <- model.matrix(X.formula, data)

  row_subject <- intersect(row_subject, rownames(X))

  # outcome
  Yraw <- na.omit(data[,as.character(formul$fixed)[2], drop = FALSE])
  Y <- lcmm::predictlink(model, Yvalues = Yraw[,1], ndraws = 0)$pred[order(order(Yraw)), 2, drop = FALSE]
  rownames(Y) <- rownames(Yraw)

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
              call = as.formula(paste(formul$fixed[2], formul$fixed[1], formul$fixed[3],
                                      "+ (", formul$random[2], "|", formul$subject, ")")),
              group = subject))

}
