#' Function to predict the random effects on data using hlme output
#'
#' @param model A \code{hlme} object
#' @param data A dataframe containing longitudinal data
#'
#' @return
#'
#' @examples
#'
predRE <- function(model, data){

  formul <- model$call

  subject <- formul$subject

  na_subject <- NULL

  num_fixed_effect <- length(model$Xnames)
  num_random_effect <- ncol(model$predRE) - 1

  num_param_fixed <- length(model$Xnames)
  num_param_effect <- (num_random_effect * (num_random_effect+1))/2

  # fixed regression coefficients
  beta <- matrix(model$best[1:num_param_fixed], ncol = 1)
  rownames(beta) <- names(model$best[1:num_param_fixed])
  colnames(beta) <- "beta"

  vcov_re <- model$best[(num_param_fixed+1):(num_param_fixed+num_param_effect)]

  # variance covariance random effects
  B <- matrix(NA, nrow = num_random_effect, ncol = num_random_effect, byrow = TRUE)

  ind_var_re <- 1

  for (ind_row in 1:ncol(B)){
    for (ind_col in 1:ind_row){
      B[ind_row, ind_col] <- vcov_re[ind_var_re]
      ind_var_re <- ind_var_re + 1
    }
  }

  B[upper.tri(B, diag = F)] <- B[lower.tri(B, diag = F)]

  # stardard error
  se <- model$best[length(model$best)]

  # random design matrix

  Z_formula <- as.formula(formul$random)
  Z <- model.matrix(Z_formula, data)

  na_subject <- c(na_subject, setdiff(rownames(data), rownames(Z)))

  # à automatiser en fonction des variables du model
  # fixed design matrix

  X_formula <- as.formula(as.character(formul$fixed)[-2])
  X <- model.matrix(X_formula, data)

  na_subject <- c(na_subject, setdiff(rownames(data), rownames(X)))

  # outcome
  Y_formula <- as.formula(paste(as.character(formul$fixed)[1], "-1+", as.character(formul$fixed)[2]))
  Y <- model.matrix(Y_formula, data)

  na_subject <- c(na_subject, setdiff(rownames(data), rownames(Y)))

  data <- data[which(!rownames(data)%in%na_subject),]

  predRE <- matrix(NA, nrow = length(unique(data[,subject])), ncol = num_random_effect,
                   dimnames = list(unique(data[,subject]), colnames(Z)))
  predRE_row <- 1

  for (ind_subject in unique(data[,subject])){

    ind <- which(data[,subject]==ind_subject)

    tot_subject <- nrow(data[ind,])

    Z_i <- matrix(Z[ind,], nrow = tot_subject)
    V_i <- Z_i%*%B%*%t(Z_i) + se^2*diag(tot_subject)
    Y_i <- Y[ind,]
    X_i <- X[ind,]
    b_i <- B%*%t(Z_i)%*%solve(V_i)%*%(Y_i-X_i%*%beta)

    predRE[predRE_row,] <- b_i

    predRE_row <- predRE_row + 1
  }

  return(list(b_i = predRE,
              beta = beta,
              formul = model$call))

}
