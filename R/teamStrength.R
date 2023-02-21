getGsVec <- function(method = c("Cluster", "Assigned", "Brute"), nTeams = 2, lmax = 10, topoFiles = NULL) {
    method <- match.arg(method)
    if (is.null(topoFiles))
        topoFiles <- list.files(".", ".topo$")
    teamsKey <- c(".BruteTeams", ".teams", ".AssignedTeams")
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
        topoDfO <- read_delim(topoFile, delim = " ", col_types = cols())
        gL <- lapply(group, function(g1) {
            sapply(group, function(g2) {
                
                if (length(g1) == length(g2) && all(g1 == g2)) {
                    topoDf <- topoDfO %>% filter(Source %in% g1, Target %in% g1)
                }
                else {
                    topoDf1 <- topoDfO %>% filter(Source %in% g1, Target %in% g2)
                    topoDf2 <- topoDfO %>% filter(Source %in% g2, Target %in% g1)
                    topoDf <- rbind.data.frame(topoDf1, topoDf2)
                }
                if (nrow(topoDf) < 1)
                    return(NA)
                write_delim(topoDf, paste0(net, "_gL.topo"), quote = "none")
                df <- InfluenceMatrix(paste0(net, "_gL.topo"), lmax, write = F)
                # file.remove(paste0(net, "_gL.topo"))
                if (is.null(df)) return(NA)
                (df %>% sum)/(length(g1)*length(g2))
            })
        }) %>% unlist
        glFiles <- list.files(".", "_gL.topo")
        sapply(glFiles, file.remove)
        if (is.null(nTeams) && method == "Cluster")
            return(c(gs, gW))
        else
            return(c(gL, gs, gW))
    }) %>% t
    
    df <- data.frame(Network = str_remove(topoFiles, ".topo"))
    gL <- NULL
    if (!is.null(nTeams)) {
        nTeamz <- sqrt(ncol(GVals)-2)
        gL <- lapply(1:nTeamz, function(i) {
            sapply(1:nTeamz, function(j) {
                paste0("G", i, j)
            })
        }) %>% unlist
    }
    
    df <- cbind.data.frame(df, GVals) %>% set_names(c("Network", gL, "Gs", "GW"))
    dir.create("CompiledData")
    write_csv(df, paste0("CompiledData/",method,"Teams.csv"), quote = "none")
}

plotTeams <- function(topoFile, method = c("Cluster", "Assigned", "Brute"), 
    lmax = 10, force = F, write = F) {
    method <- match.arg(method)
    teamsKey <- c(".BruteTeams", ".teams", ".AssignedTeams")
    names(teamsKey) <- c("Brute", "Cluster", "Assigned")
#     nOrder <- findTeams(topoFile) %>% unlist
    net <- topoFile %>% str_remove(".topo")
    teamsFile <- paste0(net, teamsKey[method])
    if (file.exists(teamsFile)) {
        nOrder <- readLines(teamsFile) %>% str_split(",") %>% unlist
    }
    else {
        print("GTFO")
        return()
    }
    ls <- TopoToIntMat(topoFile)
    intmat <- ls[[1]]
    nodes <- ls[[2]]
    inflMat <- InfluenceMatrix(topoFile, lmax = lmax, force = force, write = write)
    df2 <- data.frame(inflMat) %>%
        mutate(nodes1 = rownames(.)) %>%
        gather(key = "Nodes", value = "Influence", -nodes1) %>%
        mutate(nodes1 = factor(nodes1, levels = nOrder), Nodes = factor(Nodes, levels = nOrder))
    textSize <- ifelse(length(nOrder) > 28, 0.2, 1)
    p <- ggplot(df2, aes(x = Nodes, y = nodes1, fill = Influence)) + geom_tile() +
        theme_Publication() + scale_fill_gradient2(low = "blue", high = "red", limits = c(-1,1)) +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              axis.text = element_text(size = rel(textSize)),
              legend.position = "right",
              legend.direction = "vertical", legend.key.height = unit(0.5, "in"))

    DirectoryNav("MatrixPlots")
    ggsave(paste0(net, "_", lmax, ".png"), width = 7.5, height = 6)
    # logDf <- logDf %>% mutate(InfluencePlot = ifelse(Files == topoFile, "Yes", InfluencePlot))
    setwd("..")
}

plotCustomColors <- function(topoFiles = NULL, negCol = "#696969", 
    T1 = "#ff0000", T2 = "#0000ff", zeroColr = "#ffffff",
    alpha = T) {
    # DirectoryNav("influenceTest")
    # file.copy(paste0("../", topoFile), ".")
    if (is.null(topoFiles)) 
        topoFiles <- list.files(".", ".topo$")
    if (length(topoFiles) == 0 || !all(topoFiles %>% sapply(file.exists)))
    {
        message("Invalid input topofiles. Do topofiles exist in your folder?")
        return()

    }
    getGsVec(method = "Cluster", topoFiles = topoFiles)
    sapply(topoFiles, function(topoFile) {
        inflMat <- read_csv(paste0("Influence/", 
            str_replace(topoFile, ".topo", "_reducedInfl.csv")))
        teams <- readLines(str_replace(topoFile, ".topo", ".teams")) %>%
            str_split(",")
        teamsOrder <- teams %>% unlist
        df <- inflMat %>% 
            gather(key = "Target", value = "Influence", -Source) %>%
            mutate(colorVal = negCol, alphaVal = 1)
        if (alpha) {
            df <- df %>%
                mutate(alphaVal = abs(Influence),
                    colorVal = ifelse(Source %in% teams[[1]] & 
                        Target %in% teams[[1]], T1, colorVal)) %>%
                mutate(colorVal = ifelse(Source %in% teams[[2]] & 
                    Target %in% teams[[2]], T2, colorVal))
        }
        df <- df  %>%
            mutate(Target = factor(Target, levels = teamsOrder),
                Source = factor(Source, levels = teamsOrder))
        ggplot(df, aes(x = Target, y = Source, fill = colorVal, 
            alpha = alphaVal)) +
            geom_tile() + 
            theme_Publication() + 
            scale_fill_identity() + 
            scale_alpha_identity() +
            theme(axis.text.x = element_text(angle = 90, 
                hjust = 1, vjust = 0.5))

        DirectoryNav("MatrixPlots")
        ggsave(str_replace(topoFile, 
            ".topo", "_customInfluence.png"), width = 6.5, height = 6)
        setwd("..")
    })
    
        
}
