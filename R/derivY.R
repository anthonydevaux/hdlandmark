#' Function to predict the slope of longitudinal outcome
#'
#' @param predRE A list object from \code{predRE} function
#' @param data A dataframe where each row containing some predictive variables for a specific subject
#' @param derivForm A list containing the derivation form
#'
#' @return A matrix containing the prediction slope of the longitudinal outcome for each subject
#'
#' @examples
#'
derivY <- function(predRE, data, derivForm){

  subject <- predRE$formul$subject

  na_subject <- NULL

  beta_deriv <- predRE$beta[derivForm$indFixed, , drop = FALSE]
  b_deriv <- t(predRE$b_i[,derivForm$indRandom])

  X_deriv_formula <- derivForm$fixed
  X_deriv <- model.matrix(X_deriv_formula, data)

  na_subject <- c(na_subject, setdiff(rownames(data), rownames(X_deriv)))

  Z_deriv_formula <- derivForm$random
  Z_deriv <- model.matrix(Z_deriv_formula, data)

  na_subject <- c(na_subject, setdiff(rownames(data), rownames(Z_deriv)))

  data <- data[which(!rownames(data)%in%na_subject),]

  data <- data[which(data[,subject]%in%intersect(data[,subject], colnames(b_deriv))),]

  Y_deriv <- matrix(NA, nrow = length(unique(data[,subject])), ncol = 1,
                    dimnames = list(unique(data[,subject]), "Y_deriv"))

  Y_deriv_row <- 1

  for (ind_subject in unique(data[,subject])){

    ind <- rownames(data[which(data[, subject] == ind_subject),])

    b_i_deriv <- b_deriv[,ind]
    X_i_deriv <- X_deriv[ind,]
    Z_i_deriv <- Z_deriv[ind,]

    Y_i_deriv <- X_i_deriv%*%beta_deriv + Z_i_deriv%*%b_i_deriv

    Y_deriv[Y_deriv_row,] <- Y_i_deriv

    Y_deriv_row <- Y_deriv_row + 1

  }

  return(Y_deriv)
}
