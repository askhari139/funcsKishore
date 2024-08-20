formatRACIPE <- function(topoFile, force = F, gkNorm = T) {
    # setwd("RACIPE")
    #### Setup ------
    net <- topoFile %>% str_remove(".topo")
    paramNames <- read.delim(paste0(net, ".prs"), sep = "") %>%
        select(Parameter) %>%
        unlist
    nodes <- paramNames[str_detect(paramNames, "Prod_of")] %>%
        str_remove("Prod_of_")
    if (file.exists(paste0(net, "_formatRacipe.csv")) && !force) {
        print("Previous Formatting exists. Use force=T if recalculation is needed")
        dat <- read_csv(paste0(net, "_formatRacipe.csv"), show_col_types = F)
        normFactor <- read_csv(paste0(net, "_gkNormFactor.csv"), show_col_types = F)
    }
    ### Formatting ------
    else{
        topoDf <- read.delim(topoFile, sep = "", stringsAsFactors = F)
        deletables <- paste0(net, c("_solution_", "_T_test_", ".cfg"))
        fileList <- sapply(deletables, function(x) {
            fAll <- list.files(pattern = x)
            sapply(fAll, file.remove)
        })
        if(!file.exists(paste0(net, "_solution.dat"))) {
            print("Solution file does not exits. Simulate!")
            # setwd("..")
            return(NA)
        }
        dat <- read_delim(paste0(net, "_solution.dat"), 
            delim = "\t", col_names = F, col_types = "d")
        colnames(dat) <- c("paramID", "nStates", "Basin", nodes)
        parameters <- read_delim(paste0(net, "_parameters.dat"),
            delim = "\t", col_names = F, col_types = "d") %>%
            set_names(c("paramID", "nStates", paramNames))
        normFactor <- (parameters %>%
            select(all_of(paste0("Prod_of_", nodes))) %>% set_names(nodes))/(
                parameters %>% select(all_of(paste0("Deg_of_", nodes))) %>% set_names(nodes)) 
        normFactor <- normFactor %>%
            mutate(across(everything(), log2)) %>%
            mutate(paramID = parameters$paramID) %>%
            merge(dat %>% select(paramID), by = "paramID")
        
        write_csv(normFactor, paste0(net, "_gkNormFactor.csv"))
        write_csv(dat, paste0(net, "_formatRacipe.csv"))
    }
    if (gkNorm) {
        sapply(nodes, function(x) {
            dat[[x]] <<- dat[[x]] - normFactor[[x]]
        })
    }
    return(dat)
}

zNormDf <- function(df, nodes, nLevels = NULL) {
    df %>% select(all_of(nodes)) %>% mutate(across(everything(), .fns = function(x){
        ifelse(x-mean(x) >0, 1, 0)
        }))
}

minMaxDf <- function(df, nodes, nLevels=100) {
    df %>% select(all_of(nodes)) %>% mutate(across(everything(), .fns = function(x){
        x <- (x-min(x))/(max(x) - min(x))
        if (is.na(nLevels) || !is.numeric(nLevels))
            nLevels <- 100
        
        y <- round(x*nLevels)/nLevels
        return(y)
        }))

}

discretize <- function(topoFile, method = c("zScore", "minMax"), 
    nLevels = NULL, gkNorm = F, getFrequency = F) 
    {
        net <- topoFile %>% str_remove(".topo")
        formatData <- formatRACIPE(topoFile, gkNorm = gkNorm)
        gkNormFactor <- read_csv(paste0(net, "_gkNormFactor.csv"), show_col_types = F)
        nodes <- colnames(gkNormFactor)
        method <- method[1]
        if (method == "zScore")
            f <- zNormDf
        else if (method == "minMax") {
            if (is.null(nLevels))
                nLevels <- 100
            f <- minMaxDf
        }
        else {
            stop("Wrong discretization method")
        }
        df <- f(formatData, nodes, nLevels)
        if (getFrequency) {
            df$Basin <- formatData$Basin
            maxBasin <- max(max(df$Basin), 100)
            numParas <- max(formatData$paramID)
            df <- df %>% 
                mutate(Frequency = Basin/(numParas*maxBasin)) %>%
                group_by(across(all_of(nodes))) %>%
                summarise(Frequency = sum(Frequency), .groups = "drop")
        }
        return(df)
    }

teamScores <- function(df, teams) {
    nodes <- unlist(teams)
    if (is.null(names(teams)))
        names(teams) <- paste0("T_", 1:length(teams))
    scoreData <- sapply(teams, function(t) {
        df %>% select(all_of(t)) %>% rowMeans
    }) %>% data.frame %>% set_names(names(teams))
    return(cbind.data.frame(df, scoreData))
}

# Get the state table
# Organize the nodes by teams
# cluster the rows
# gather the data into long format
# assign team names for all nodes
# plot the heatmap faceted by the teams

StateFrequencyMap <- function(freqData, teams, key, figureDir = ".", 
    freqWidth = F, ord = NULL, orderBy = NULL) {
    # browser()
    if (!dir.exists(figureDir))
        dir.create(figureDir, recursive = T)
    nodeOrder <- teams %>% unlist %>% unname
    teamLabels <- rep(names(teams), sapply(teams, length))
    names(teamLabels) <- nodeOrder
    freqOrdering <- T
    if (is.null(freqData$Frequency)) {
        freqData$Frequency <- 1
        freqOrdering <- F
    }
    if (!is.null(ord)) {
        freqData$num <- ord
    }
    freqData <- freqData %>%
        group_by(across(all_of(nodeOrder))) %>%
        summarise(Frequency = sum(Frequency), 
            num = min(num),
            .groups = "drop")
        # select(all_of(nodeOrder), Frequency)
    if (is.null(ord)) {
        if (is.null(orderBy)) {
            orderBy <- nodeOrder
        }
        hc <- hclust(freqData %>% select(all_of(orderBy)) %>% dist)
        ord <- hc$order
        freqData$num <- ord
    }
    
    if (freqOrdering) {
        freqData <- freqData %>% arrange(desc(Frequency), num) %>%
            mutate(num = 1:nrow(.))
    }
    freqData <- freqData %>%
        pivot_longer(cols = all_of(nodeOrder), names_to = "Nodes", values_to = "Level") %>%
        mutate(Teams = factor(teamLabels[Nodes], levels = names(teams))) %>%
        mutate(Nodes = factor(Nodes, levels = nodeOrder))
    # browser()
    nLevels <- length(unique(freqData$Level))
    if (nLevels == 2) {
        levels <- unique(freqData$Level) %>% sort
        freqData$Level <- factor(freqData$Level, levels = levels)
        scaleFunc <- scale_fill_manual(values = c("red", "blue"))
    }
    else {
        scaleFunc <- scale_fill_gradient(low = "red", high = "blue")
    }
    if (freqWidth) {
        freqData$Frequency <- log10(freqData$Frequency)
        freqData$Frequency <- freqData$Frequency - min(freqData$Frequency)
        p <- ggplot(freqData, aes(x = Nodes, y = num, fill = Level)) +
        geom_tile(aes(height = 0.1 + Frequency / max(Frequency)))
    }
    else {
        p <- ggplot(freqData, aes(x = Nodes, y = num, fill = Level)) +
        geom_tile()
    }
    w <- min(2 + 3 * length(teams), 15)
    textSize <- min(0.9, 5 * (w - 2.5)/length(nodeOrder))
    
    p <- p +
        theme_Publication() +
        scaleFunc +
        scale_y_discrete(expand=c(0,0)) +
        theme(axis.text.y = element_blank(), 
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = rel(textSize)),
            legend.key.width = unit(0.8, "cm")) +
        labs(x = "Node", y = "", fill = "Expression Level  ") +
        facet_grid(cols = vars(Teams), scales = "free", space = "free")
    ggsave(paste0(figureDir, "/", key, ifelse(freqWidth, "_fw", ""), "_stateFrequencyMap.png"), 
            height = max(min(4.5 + nrow(freqData) * 0.003, 10), 20), width = min(2 + 3 * length(teams), 10))
}


scorePlots <- function(df, teams, key, figureDir = ".", freqWidth = F, hLim = 15) {
    # browser()
    if(!dir.exists(figureDir))
        dir.create(figureDir, recursive = T)
    if (is.null(df$Frequency)) {
        df$Frequency <- 1
        freqWidth <- F
    }
    nm <- names(teams)
    df <- df %>% 
        mutate(across(all_of(nm), .fns = function(x){round(x, 2)})) %>%
        group_by(across(all_of(nm))) %>%
        summarise(Frequency = sum(Frequency), .groups = "drop")
    dfOrig <- df
    # browser()
    ## heatmap
    # df <- df %>% 
    #     arrange(across(all_of(names(teams))), desc(Frequency)) %>%
    #     mutate(num = 1:nrow(.))
    nS <- nrow(df)
    
    h <- max(4.5 + nS * 0.003, hLim)
    if (h > hLim) {
        nS <- (hLim-4.5)/0.003
        nSnew <- nS/2
        df <- rbind.data.frame(tail(df, nSnew), head(df, nSnew))# %>%
            # arrange(across(all_of(names(teams))), desc(Frequency)) %>%
            # mutate(num = 1:nrow(.))
        h <- 4.5 + nS*0.003
    }
    ord <- hclust(df %>% select(-Frequency) %>% dist)
    ord <- ord$order
    df$num <- ord
    df <- df %>%
        arrange(desc(Frequency), num) %>%
        mutate(num = 1:nrow(.))
    df <- df %>%
        gather(key = "Team", value = "Score", -num, -Frequency)
    
    
    
    if (freqWidth) {
        p <- ggplot(df, aes(x = Team, y = reorder(num, Frequency), fill = Score))
        p <- p + geom_tile(aes(height = 0.1 + Frequency / max(Frequency)))
    }
        
    else {
        p <- ggplot(df, aes(x = Team, y = num, fill = Score))
       p <- p + geom_tile()
    }
    p <- p + 
        scale_fill_viridis_c() +
        scale_y_discrete(expand=c(0,0)) +
        theme_Heatmap()+
        theme(axis.text.y = element_blank(), 
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        labs(x = "Team", y = "", fill = "Score")
    # browser()
    ggsave(paste0(figureDir, "/", key, ifelse(freqWidth, "_fW", ""), "_scoreHeatmap.png"), 
    height = h, width = 2+ 0.5*length(teams))

    ## boxplot
    ggplot(df, aes(x = Team, y = Score)) +
        geom_boxplot() + 
        theme_Publication()
    ggsave(paste0(figureDir, "/", key, "_scoreBoxplot.png"), height = 5.2, width = 5.5)
    # browser()
    ## scatterplots
    if (length(teams) > 2) {
        png(paste0(figureDir, "/", key, "_pairWisePlots.png"), 
            width = 3.5*length(teams), height = max(2.5*length(teams), 10), 
            units = "in", res = 150)
        p <- ggpairs(dfOrig %>% 
            select(all_of(names(teams)), Frequency) %>%
            mutate(Frequency = as.character(round(log10(Frequency)))), 
            columns = 1:length(teams), 
            mapping = ggplot2::aes(color = Frequency, alpha = 0), 
            upper = "blank", legend = length(teams) + 1,
            diag = list(continuous = wrap("barDiag", bins = 15, 
                position = position_dodge2())))
        try(print(p))
        dev.off()
        # density
        png(paste0(figureDir, "/", key, "_pairWisePlotsDens.png"), 
            width = 3.5*length(teams), height = max(2.5*length(teams), 10), 
            units = "in", res = 150)
        p <- ggpairs(dfOrig %>% 
            select(all_of(names(teams)), Frequency) %>%
            mutate(Frequency = round(Frequency*10000)) %>%
            uncount(Frequency), 
            columns = 1:length(teams), 
            # mapping = ggplot2::aes(color = Frequency, alpha = 0), 
            # upper = "blank", 
            legend = length(teams) + 1,
            diag = list(continuous = wrap("barDiag", bins = 15)), 
            lower = list(continuous = "density"))
        try(print(p))
        dev.off()
    }
    

}

heatMapsRACIPE <- function(topoFile, nodeOrder = NULL, method = c("zScore", "minMax"), 
    nLevels = NULL, gkNorm = F, getFrequency = T, discretize = F, plotFrequency = F, 
    teams = NULL, teamKey = NULL, phenotypeFunction = NULL)
    {
        browser()
        net <- topoFile %>% str_remove(".topo")
        # formatData <- formatRACIPE(topoFile, gkNorm = gkNorm)
        gkNormFactor <- read_csv(paste0(net, "_gkNormFactor.csv"), show_col_types = F)
        key <- paste0(ifelse(gkNorm, "_gk", ""),
            ifelse(discretize, "_discret", ""), 
            paste0("_", method[1]), 
            teamKey)

        nodes <- colnames(gkNormFactor)
        if (is.null(nodeOrder)) {
            if (!file.exists(paste0(net, ".teams"))) {
                getGsVec()
            }
            nodeOrder <- readLines(paste0(net, ".teams")) %>% 
                str_split(",") %>% unlist
        }
        if (is.null(teams)) {
            teams <- readLines(paste0(net, ".teams")) %>% str_split(",")
        }
        if(is.null(names(teams))) {
            names(teams) <- paste0("T_", 1:length(teams))
        }
        method <- method[1]
        if (discretize) {
            df <- discretize(topoFile, method, nLevels, gkNorm, getFrequency) %>%
                select(all_of(nodes), Frequency)
        }
        else {
            formatData <- formatRACIPE(topoFile, gkNorm = gkNorm)
            df <- formatData %>% select(all_of(nodes)) %>%
                mutate(across(everything(), .fns = function(x){
                    if (method == "zScore")
                        (x - mean(x))/sd(x)
                    else if (method == "minMax") {
                        x <- x-min(x)/(max(x) - min(x))
                    }
                    else {
                        (x - mean(x))/sd(x)
                    }
                }))
        }
        df <- df %>%
            teamScores(teams)
        
        if (discretize) {
            df <- df %>% 
                arrange(across(all_of(names(teams))), desc(Frequency)) %>% 
                mutate(Phenotype = ifelse(is.null(phenotypeFunction), "", 
                phenotypFunction(df %>% select(all_of(names(teams))))))
            scoreDf <- df %>% select(all_of(names(teams)), Frequency)
            scorePlots(scoreDf, teams, key, 
                figureDir = paste0("ScoreFigures/", net, "/RACIPE"))
            scorePlots(scoreDf, teams, key, 
                figureDir = paste0("ScoreFigures/", net, "/RACIPE"), 
                freqWidth = T)
            StateFrequencyMap(df, teams, key,
                figureDir = paste0("ScoreFigures/", net, "/RACIPE"), 
                freqWidth = T)
            StateFrequencyMap(df, teams, key,
                figureDir = paste0("ScoreFigures/", net, "/RACIPE"))
            # HeatMapDiscrete(df)
        }
        else {
            df <- df %>% 
                arrange(across(all_of(names(teams)))) %>% 
                mutate(Phenotype = ifelse(is.null(phenotypeFunction), "", 
                phenotypFunction(df %>% select(all_of(names(teams))))))
            # HeatMapDiscrete(df, teams)
            scoreDf <- df %>% select(all_of(names(teams)))
            scorePlots(scoreDf, teams, key, 
                figureDir = paste0("ScoreFigures/", net, "/RACIPE"))
            StateFrequencyMap(df, teams, key,
                figureDir = paste0("ScoreFigures/", net, "/RACIPE"))
        }
    }
