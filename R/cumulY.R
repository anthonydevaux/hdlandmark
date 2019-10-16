cumulY <- function(predRE, data, time, tLM){

  browser()

  subject <- predRE$formul$subject

  Y_cumul <- matrix(NA, nrow = length(unique(data[,subject])), ncol = 1,
                    dimnames = list(unique(data[,subject]), "Y_cumul"))

  Y_cumul_row <- 1

  for (ind_subject in unique(data[,subject])){

    newdata_id <- data[data$id==ind_subject,]

    predRE$b_i <- predRE$b_i[which(rownames(predRE$b_i)==ind_subject),, drop = FALSE]

    F2 <- function(predRE, data, time){
      f2 <- function(t){
        predY(predRE, data, time, tLM = t)
      }
      return(f2)
    }
    int2 <- function(predRE, data, time){
      integrate(F2(predRE, data, time), lower = 0, upper = tLM)$value
    }

    Y_cumul[Y_cumul_row,] <- integrate(F2(predRE, data, time), lower = 0, upper = 4)$value

    Y_cumul_row <- Y_cumul_row + 1

  }

  return(Y_cumul)

}


# test <- function(tLM){
#   sapply(tLM, function(predRE, data, time){
#     predY(predRE, data, time)
#   })
# }
