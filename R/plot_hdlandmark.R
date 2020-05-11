#' Title
#'
#' @param hdlandmark.output
#'
#' @return
#' @export
#'
#' @import ggplot2 viridis
#'
#' @examples
plot.hdlandmark <- function(hdlandmark.output){

  tLMs <- hdlandmark.output$tLMs
  tHors <- hdlandmark.output$tHors

  tab.plot <- data.frame()

  if (hdlandmark.output$kfolds > 1){

    for (tLM in tLMs){

      for (tHor in tHors){

        AUC.mean <- colMeans(hdlandmark.output$models[[as.character(tLM)]]$AUC[[as.character(tHor)]])
        BS.mean <- colMeans(hdlandmark.output$models[[as.character(tLM)]]$BS[[as.character(tHor)]])

        AUC.LB <- AUC.mean+qnorm(0.025)*apply(hdlandmark.output$models[[as.character(tLM)]]$AUC[[as.character(tHor)]], 2, sd)
        AUC.LB <- ifelse(AUC.LB<0,0,AUC.LB)

        AUC.UB <- AUC.mean+qnorm(0.975)*apply(hdlandmark.output$models[[as.character(tLM)]]$AUC[[as.character(tHor)]], 2, sd)
        AUC.UB <- ifelse(AUC.UB>1,1,AUC.UB)

        BS.LB <- BS.mean+qnorm(0.025)*apply(hdlandmark.output$models[[as.character(tLM)]]$BS[[as.character(tHor)]], 2, sd)
        BS.LB <- ifelse(BS.LB<0,0,BS.LB)

        BS.UB <- BS.mean+qnorm(0.975)*apply(hdlandmark.output$models[[as.character(tLM)]]$BS[[as.character(tHor)]], 2, sd)
        BS.UB <- ifelse(BS.UB>1,1,BS.UB)

        df <- data.frame(tLM = tLM, tHor = tHor, method = hdlandmark.output$models.name,
                         AUC.mean = AUC.mean, AUC.LB = AUC.LB, AUC.UB = AUC.UB,
                         BS.mean = BS.mean, BS.LB = BS.LB, BS.UB = BS.UB)

        tab.plot <- rbind(tab.plot, df)

      }

    }

    pd <- position_dodge(0.2)

    g.AUC <- ggplot(tab.plot, aes(x = tHor, y = AUC.mean, group = method, color = method)) +
      geom_line(size = 1) +
      xlab("Horizon time") +
      ylab("AUC") +
      scale_color_manual("Methods", values = viridis(length(unique(tab.plot$method)))) +
      facet_grid(tLM ~ .) +
      theme_bw()

    g.BS <- ggplot(tab.plot, aes(x = tHor, y = BS.mean, group = method, color = method)) +
      geom_line(size = 1) +
      xlab("Horizon time") +
      ylab("BS") +
      scale_color_manual("Methods", values = viridis(length(unique(tab.plot$method)))) +
      facet_grid(tLM ~ .) +
      theme_bw()

    # g.AUC <- ggplot(tab.plot, aes(x = tHor, y = AUC.mean, group = method, color = method)) +
    #   geom_point(position = pd, size = 3) +
    #   geom_errorbar(aes(ymin = AUC.LB, ymax = AUC.UB), width = 0.1, position = pd, size = 1) +
    #   xlab("Horizon time") +
    #   ylab("AUC (95% CI)") +
    #   scale_color_manual("Methods", values = viridis(length(unique(tab.plot$method)))) +
    #   facet_grid(tLM ~ .) +
    #   theme_bw()
    #
    # g.BS <- ggplot(tab.plot, aes(x = tHor, y = BS.mean, group = method, color = method)) +
    #   geom_point(position = pd, size = 3) +
    #   geom_errorbar(aes(ymin = BS.LB, ymax = BS.UB), width = 0.1, position = pd, size = 1) +
    #   xlab("Horizon time") +
    #   ylab("BS (95% CI)") +
    #   scale_color_manual("Methods", values = viridis(length(unique(tab.plot$method)))) +
    #   facet_grid(tLM ~ .) +
    #   theme_bw()

  }else{



  }

  print(g.AUC)
  print(g.BS)

}
