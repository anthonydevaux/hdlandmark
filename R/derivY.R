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

  if (is.null(derivForm$fixed) | is.null(derivForm$indFixed)){

    stop("Fixed or indFixed argument missing in derivForm !", "\n")

  }

  subject <- predRE$group

  id_subject <- unique(data[,subject])

  beta_deriv <- predRE$beta[derivForm$indFixed]

  X_deriv_formula <- derivForm$fixed
  X_deriv <- model.matrix(X_deriv_formula, data)

  id_subject <- intersect(id_subject, data[rownames(X_deriv),subject])

  if (!is.null(derivForm$random) & !is.null(derivForm$indRandom)){

    b_deriv <- predRE$b_i[,derivForm$indRandom, drop = FALSE]

    Z_deriv_formula <- derivForm$random
    Z_deriv <- model.matrix(Z_deriv_formula, data)

    id_subject <- intersect(id_subject, data[rownames(Z_deriv),subject])

    id_subject <- intersect(id_subject, rownames(b_deriv))

    add_Z <- TRUE

  }else{

    b_deriv <- NULL
    Z_deriv <- NULL

    add_Z <- FALSE
  }

  data <- data[which(data[,subject]%in%id_subject),]

  Y_deriv <- matrix(NA, nrow = length(unique(data[,subject])), ncol = 1,
                    dimnames = list(unique(data[,subject]), "Y_deriv"))

  Y_deriv_row <- 1

  for (ind_subject in unique(data[,subject])){

    ind <- rownames(data[which(data[, subject] == ind_subject),])

    b_i_deriv <- b_deriv[which(rownames(b_deriv)==ind_subject), ]
    X_i_deriv <- X_deriv[ind,]
    Z_i_deriv <- Z_deriv[ind,]

    Y_i_deriv <- X_i_deriv%*%beta_deriv

    if (add_Z){ # ajout random effects ?

      Y_i_deriv <- Y_i_deriv + Z_i_deriv%*%b_i_deriv

    }

    Y_deriv[Y_deriv_row,] <- Y_i_deriv

    Y_deriv_row <- Y_deriv_row + 1

  }

  return(Y_deriv)
}
