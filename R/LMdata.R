#' Title
#'
#' @param data
#' @param time
#' @param time_event
#' @param event
#' @param tHor
#'
#' @return
#' @export
#'
#' @examples
#'
LMdata <- function(data, time, time_event, event, tHor){

  data[,time] <- data[,time] - 4
  data[,time_event] <- data[,time_event] - 4
  data[which(data[,time_event] > tHor), event] <- 0
  data[,time_event] <- pmin(data[,time_event], tHor)

  return(data)

}
