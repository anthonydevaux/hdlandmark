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
  b <- predRE$b_i

  subject <- predRE$formul$subject

  id_subject <- unique(data[,subject])

  num_fixed_effect <- length(beta)
  num_random_effect <- nrow(b)

  X_formula <- as.formula(as.character(formul$fixed)[-2])
  X <- model.matrix(X_formula, data)

  id_subject <- intersect(id_subject, data[rownames(X),subject])

  Z_formula <- as.formula(formul$random)
  Z <- model.matrix(Z_formula, data)

  id_subject <- intersect(id_subject, data[rownames(Z),subject])

  id_subject <- intersect(id_subject, rownames(b))

  data <- data[which(data[,subject]%in%id_subject),]

  Y_pred <- matrix(NA, nrow = length(unique(data[,subject])), ncol = 1,
                  dimnames = list(unique(data[,subject]), "Y_pred"))

  Y_pred_row <- 1

  for (ind_subject in unique(data[,subject])){

    ind <- rownames(data[which(data[, subject] == ind_subject),])

    b_i <- b[which(rownames(b)==ind_subject), ]
    X_i <- X[ind,]
    Z_i <- Z[ind,]

    Y_i_pred <- X_i%*%beta + Z_i%*%b_i

    Y_pred[Y_pred_row,] <- Y_i_pred

    Y_pred_row <- Y_pred_row + 1

  }

  return(Y_pred)
}
