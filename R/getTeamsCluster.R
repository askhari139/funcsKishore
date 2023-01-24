teamStrengthF <- function(group, inflMat) {
    gL <- lapply(group, function(g1) {
            sapply(group, function(g2) {
                (inflMat[g1, g2] %>% sum)/(length(g1)*length(g2))
            })
        }) %>% unlist
        gL %>% abs %>% mean
}

findClusterTeams <- function(topoFile, nTeams = 2, lmax = 10, nTeamSample = 2:10) {
    # browser()
    net <- topoFile %>% str_remove(".topo")
    ls <- TopoToIntMat(topoFile)
    intmat <- ls[[1]]
    nodes <- ls[[2]]
    inflMat <- InfluenceMatrix(topoFile, lmax)
    nodes <- rownames(inflMat)
    df <- inflMat
    df1 <- apply(df, 2, function(x) {
        ifelse(x > 0, 1, -1)
    })
    df1 <- cbind(df1, t(df1))
    hc <- hclust(dist(df1))
    if (is.null(nTeams)) {
        gs <- sapply(nTeamSample, function(nT) {
            clust <- cutree(hc, nT)
            group <- lapply(1:nTeams, function(x) {
                nodes[clust==x] %>% sort
            }) 
            teamStrengthF(group, df)
        })
        nTeams <- nTeamSample[which.max(gs)]
    }
    clust <- cutree(hc, nTeams)
    group <- lapply(1:nTeams, function(x) {
        nodes[clust==x] %>% sort
    })
    l <- group %>% sapply(paste0, collapse = ",") %>% sort
    writeLines(l, paste0(net, ".teams"))
    return(group)
}
