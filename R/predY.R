#' Function to predict the value of longitudinal outcome for a specific time
#'
#' @param predRE A list object from \code{predRE} function
#' @param data A dataframe where each row containing some predictive variables for a specific subject
#'
#' @return A matrix containing the prediction value of the longitudinal outcome for each subject
#'
#' @examples
#'
predY <- function(predRE, data){

  formul <- predRE$formul

  beta <- predRE$beta
  b <- t(predRE$b_i)

  subject <- predRE$formul$subject

  na_subject <- NULL

  num_fixed_effect <- length(beta)
  num_random_effect <- nrow(b)

  X_formula <- as.formula(as.character(formul$fixed)[-2])
  X <- model.matrix(X_formula, data)

  na_subject <- c(na_subject, setdiff(rownames(data), rownames(X)))

  Z_formula <- as.formula(formul$random)
  Z <- model.matrix(Z_formula, data)

  na_subject <- c(na_subject, setdiff(rownames(data), rownames(Z)))

  data <- data[which(!rownames(data)%in%na_subject),]

  data <- data[which(data[,subject]%in%intersect(data[,subject], colnames(b))),]

  Y_pred <- matrix(NA, nrow = length(unique(data[,subject])), ncol = 1,
                  dimnames = list(unique(data[,subject]), "Y_pred"))

  Y_pred_row <- 1

  for (ind_subject in unique(data[,subject])){

    ind <- which(data[,subject]==ind_subject)

    b_i <- b[,ind]
    X_i <- X[ind,]
    Z_i <- Z[ind,]

    Y_i_pred <- X_i%*%beta + Z_i%*%b_i

    Y_pred[Y_pred_row,] <- Y_i_pred

    Y_pred_row <- Y_pred_row + 1

  }

  return(Y_pred)
}
