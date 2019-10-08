predY_newdata <- function(predRE, newdata){

  formul <- predRE$formul

  beta <- predRE$beta
  b <- t(predRE$b_i)

  subject <- predRE$formul$subject

  num_fixed_effect <- length(beta)
  num_random_effect <- nrow(b)

  X_formula <- as.formula(as.character(formul$fixed)[-2])
  X <- model.matrix(X_formula, newdata)

  Z_formula <- as.formula(formul$random)
  Z <- model.matrix(Z_formula, newdata)

  predY <- matrix(NA, nrow = length(unique(newdata[,subject])), ncol = 1,
                  dimnames = list(unique(newdata[,subject]), "predY"))

  predY_row <- 1

  for (ind_subject in unique(newdata[,subject])){

    ind <- which(newdata[,subject]==ind_subject)

    b_i <- b[,ind]
    X_i <- X[ind,]
    Z_i <- Z[ind,]

    predY_i <- X_i%*%beta + Z_i%*%b_i

    predY[predY_row,] <- predY_i

    predY_row <- predY_row + 1

  }

  return(predY)
}
