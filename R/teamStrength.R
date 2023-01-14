getGsVec <- function(method = c("Cluster", "Assigned", "Brute"), nTeams = 2, lmax = 10) {
    method <- match.arg(method)
    topoFiles <- list.files(".", ".topo")
    teamsKey <- c("BruteTeams", "teams", "AssignedTeams")
    names(teamsKey) <- c("Brute", "Cluster", "Assigned")
    GVals <- sapply(topoFiles, function(topoFile) {
        print(topoFile)
        net <- topoFile %>% str_remove(".topo")
        teamsFile <- paste0(net, teamsKey[method])
        if (file.exists(teamsFile)) {
            group <- readLines(teamsFile) %>% str_split(",")
        }
        else if (method == "Brute")
            group <- findBruteTeams(topoFile, lmax)[[1]]
        else if (method == "Cluster")
            group <- findClusterTeams(topoFile, nTeams, lmax)
        else {
            message(paste("Teams are not assigned for", net))
            return()
        }
        
        ls <- TopoToIntMat(topoFile)
        intmat <- ls[[1]]
        nodes <- ls[[2]]
        inflMat <- InfluenceMatrix(topoFile, lmax)
        gL <- lapply(group, function(g1) {
            sapply(group, function(g2) {
                (inflMat[g1, g2] %>% sum)/(length(g1)*length(g2))
            })
        }) %>% unlist
        gs <- gL %>% abs %>% mean
        gW <- sapply(group, function(g) {
            (inflMat[g, g] %>% sum)/(length(g)*length(g))
        }) %>% abs %>% mean

        gL <- lapply(group, function(g1) {
            sapply(group, function(g2) {
                topoDf <- read_delim(topoFile, delim = " ", col_types = cols())
                if (all(g1 == g2)) {
                    topoDf <- topoDf %>% filter(Source %in% g1, Target %in% g1)
                }
                else {
                    topoDf1 <- topoDf %>% filter(Source %in% g1, Target %in% g2)
                    topoDf2 <- topoDf %>% filter(Source %in% g2, Target %in% g1)
                    topoDf <- rbind.data.frame(topoDf1, topoDf2)
                }
                write_delim(topoDf, paste0(net, "_gL.topo"), quote = "none")
                InfluenceMatrix(paste0(net, "_gL.topo"), lmax, write = F) %>% sum
            })
        }) %>% unlist

        if (is.null(nTeams) && method == "Cluster")
            return(c(gs, gW))
        else
            return(c(gL, gs, gW))
    }) %>% t
    
    df <- data.frame(Network = str_remove(topoFiles, ".topo"))
    gL <- NULL
    if (is.null(nTeams)) {
        nTeamz <- sqrt(ncol(GVals)-2)
        gL <- lapply(1:nTeamz, function(i) {
            sapply(1:nTeamz, function(j) {
                paste0("G", i, j)
            })
        })
    }
    
    df <- cbind.data.frame(df, GVals) %>% set_names(c("Network", gL, "Gs", "GW"))
    dir.create("CompiledData")
    write_csv(df, paste0("CompiledData/",method,"Teams.csv"), quote = "none")
}

plotTeams <- function(topoFile) {
    nOrder <- findTeams(topoFile) %>% unlist
    net <- topoFile %>% str_remove(".topo")
    ls <- TopoToIntMat(topoFile)
    intmat <- ls[[1]]
    nodes <- ls[[2]]
    inflMat <- InfluenceMatrix(topoFile)
    df2 <- data.frame(df) %>%
        mutate(nodes1 = nodes) %>%
        gather(key = "Nodes", value = "Influence", -nodes1) %>%
        mutate(nodes1 = factor(nodes1, levels = nOrder), Nodes = factor(Nodes, levels = nOrder))
    textSize <- ifelse(length(nOrder) > 28, 0.8, 1)
    p <- ggplot(df2, aes(x = Nodes, y = nodes1, fill = Influence)) + geom_tile() +
        theme_Publication() + scale_fill_gradient2(low = "blue", high = "red", limits = c(-1,1)) +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 90),
              axis.text = element_text(size = rel(textSize)),
              legend.position = "right",
              legend.direction = "vertical", legend.key.height = unit(0.5, "in"))

    DirectoryNav("MatrixPlots")
    ggsave(str_replace(topoFile, ".topo", "_group.png"), width = 7, height = 6)
    # logDf <- logDf %>% mutate(InfluencePlot = ifelse(Files == topoFile, "Yes", InfluencePlot))
    setwd("..")
}