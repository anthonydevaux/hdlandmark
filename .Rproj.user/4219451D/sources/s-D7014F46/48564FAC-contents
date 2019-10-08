derivY_newdata <- function(predRE, newdata, derivForm){

  subject <- predRE$formul$subject

  beta_deriv <- predRE$beta[derivForm$indFixed, , drop = FALSE]
  b_deriv <- t(predRE$b_i[,derivForm$indRandom])

  X_deriv_formula <- derivForm$fixed
  X_deriv <- model.matrix(X_deriv_formula, newdata)

  Z_deriv_formula <- derivForm$random
  Z_deriv <- model.matrix(Z_deriv_formula, newdata)

  Y_deriv <- matrix(NA, nrow = length(unique(newdata[,subject])), ncol = 1,
                    dimnames = list(unique(newdata[,subject]), "Y_deriv"))

  Y_deriv_row <- 1

  for (ind_subject in unique(newdata[,subject])){

    ind <- which(newdata[,subject]==ind_subject)

    b_i_deriv <- b_deriv[,ind]
    X_i_deriv <- X_deriv[ind,]
    Z_i_deriv <- Z_deriv[ind,]

    Y_i_deriv <- X_i_deriv%*%beta_deriv + Z_i_deriv%*%b_i_deriv

    Y_deriv[Y_deriv_row,] <- Y_i_deriv

    Y_deriv_row <- Y_deriv_row + 1

  }

  return(Y_deriv)
}
