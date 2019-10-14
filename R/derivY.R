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

  id_subject <- unique(data[,subject])

  beta_deriv <- predRE$beta[derivForm$indFixed, , drop = FALSE]
  b_deriv <- predRE$b_i[,derivForm$indRandom, drop = FALSE]

  X_deriv_formula <- derivForm$fixed
  X_deriv <- model.matrix(X_deriv_formula, data)

  id_subject <- intersect(id_subject, data[rownames(X_deriv),subject])

  Z_deriv_formula <- derivForm$random
  Z_deriv <- model.matrix(Z_deriv_formula, data)

  id_subject <- intersect(id_subject, data[rownames(Z_deriv),subject])

  id_subject <- intersect(id_subject, rownames(b_deriv))

  data <- data[which(data[,subject]%in%id_subject),]

  Y_deriv <- matrix(NA, nrow = length(unique(data[,subject])), ncol = 1,
                    dimnames = list(unique(data[,subject]), "Y_deriv"))

  Y_deriv_row <- 1

  for (ind_subject in unique(data[,subject])){

    ind <- rownames(data[which(data[, subject] == ind_subject),])

    b_i_deriv <- b_deriv[which(rownames(b_deriv)==ind_subject), ]
    X_i_deriv <- X_deriv[ind,]
    Z_i_deriv <- Z_deriv[ind,]

    Y_i_deriv <- X_i_deriv%*%beta_deriv + Z_i_deriv%*%b_i_deriv

    Y_deriv[Y_deriv_row,] <- Y_i_deriv

    Y_deriv_row <- Y_deriv_row + 1

  }

  return(Y_deriv)
}
