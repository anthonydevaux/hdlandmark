cumulY <- function(predRE, newdata, var_integrate, tHW, tLM){

  subject <- predRE$formul$subject

  Y_cumul <- matrix(NA, nrow = length(unique(newdata[,subject])), ncol = 1,
                    dimnames = list(unique(newdata[,subject]), "Y_cumul"))

  Y_cumul_row <- 1

  for (ind_subject in unique(newdata[,subject])){

    newdata_id <- newdata[newdata$id==num_id,]

    f_outcome <- function(t){
      predRE$beta[1] + predRE$beta[2]*t + predRE$beta[3]*newdata_id[,"sex"] + predRE$beta[4]*newdata_id[,"sex"]*t +
        predRE$b_i[rownames(predRE$b_i)==num_id,1] + predRE$b_i[rownames(predRE$b_i)==num_id,2]*t
    }

    Y_cumul[Y_cumul_row,] <- integrate(f_outcome, lower = tLM - tHW, upper = tLM)$value

    Y_cumul_row <- Y_cumul_row + 1

  }

  return(Y_cumul)

}
