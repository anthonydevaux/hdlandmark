#' Function to predict the value of longitudinal outcome for a specific time
#'
#' @param predRE A list object from \code{predRE} function
#' @param data A dataframe where each row containing some predictive variables for a specific subject
#'
#' @return A matrix containing the prediction value of the longitudinal outcome for each subject
#' @export
#'
#' @importFrom lme4 nobars findbars
#'
#' @examples
#'
predY <- function(predRE, data, time, tLM){

  formul <- predRE$call

  beta <- predRE$beta
  b <- predRE$b_i

  subject <- predRE$group

  id_subject <- unique(data[,subject])

  data[,time] <- tLM

  X_formula <- as.formula(nobars(formul)[-2])
  X <- model.matrix(X_formula, data)

  id_subject <- intersect(id_subject, data[rownames(X),subject])

  #Z_formula <- as.formula(paste("~", as.character(findbars(formul)[[1]])[2]))
  Z_formula <- lapply(findbars(formul),
                      FUN = function(x){
                        return(as.formula(paste("~", as.character(x)[2])))
                      })

  if (length(Z_formula)>1){

    Z <- sapply(Z_formula, FUN = function(x) model.matrix(x, data))
    rownames(Z) <- rownames(X)

  }else{

    Z <- model.matrix(Z_formula[[1]], data)

  }

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

    if (!is.null(predRE$sigmae)){

      Y_i_pred <- Y_i_pred + rnorm(1,0,predRE$sigmae)

    }

    Y_pred[Y_pred_row,] <- Y_i_pred

    Y_pred_row <- Y_pred_row + 1

  }

  return(Y_pred)
}
