violinBox <- function(datViol, xVar, yVar, fillVar = NULL, violWidth = 1, # axes
    xLab = xVar, yLab = yVar, fillLab = fillVar, # axis labels
    legend_position = "top", legend_direction = "horizontal", # legend
    scale = c("viridis", "gradient"), mid = 0, # color scale
    logX = F, logY = F, # log axes
    plotKey = paste0(yVar, "_", xVar), # file name
    width = 5.5, height = 5, # file dimensions
    boxPlot = F # additional box plot
    ) {
    p <- ggplot(datViol, aes_string(x = xVar, y = yVar, fill = fillVar)) +
        labs(x = xLab, y = yLab, fill = fillLab)
    if (!is.null(fill)) {
        if (!is.factor(datViol[fill])) {
            datViol[fill] <- factor(datViol[fill], levels = unique(datViol[fill]) %>% sort)
        }
        if (scale[1] == "viridis") {
            p <- p + scale_fill_viridis_d()
        }
        else {
            p <- p + scale_fill_gradient2(midpoint = mid, low = "red", high = "blue")
        }   
    }
    p <- p + theme_Publication()
    if (logX) {
        p <- p + scale_x_log10()
    }
    if (logY) {
        p <- p + scale_y_log10()
    }
    
    pViol <- p + geom_violin(width = violWidth)
    pBox <- p + geom_boxplot()


    ggplot(datViol, aes(x = nLevels, y = Avg0, fill = PhenThresh)) +
        geom_violin(width = 1) +
        scale_y_log10() +
        labs(x = "Number of levels", y = "Steady state frequency", fill = "Phenotype") +
        theme_Publication() +
        theme(legend.position = c(0.8, 0.8),
            legend.direction = "vertical")
    ggsave(paste0(net,"_", plotKey, "_FreqDistThresh.png"), width = 5.5, height = 5)
}