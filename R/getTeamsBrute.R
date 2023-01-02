# Iterates over all nodes in the influence matrix to identify a team of nodes,
# using the definition that all the nodes in a team must have positive influence
# on each other

getTeam <- function(inflMat) {
    nodes <- colnames(inflMat)
    groupNodes <- c()
    sapply(nodes, function(node) {
        nodeList <- c(groupNodes, node)
        flag <- 0
        if (all(inflMat[nodeList, nodeList] > 0))
            flag <- 1
        if (flag)
            groupNodes <<- nodeList
    })
    groupNodes <- groupNodes %>% sort
    return(list(groupNodes, nodes[!(nodes %in% groupNodes)]))
}

# Finds teams by iteratively applying getTeams on the influence matrix and
# reducing the matrix with each iteration

findBruteTeams <- function(topoFile, twoTeams = F, lmax = 10) {
    net <- topoFile %>% str_remove(".topo")
    ls <- TopoToIntMat(topoFile)
    intmat <- ls[[1]]
    nodes <- ls[[2]]
    inflMat <- InfluenceMatrix(topoFile, lmax)
    nodes <- rownames(inflMat) %>% sort
    df <- inflMat
    df1 <- apply(df, 2, function(x) {
        ifelse(x > 0, 1, -1)
    })
    groupsData <- list()
    nGroups <- 0
    group <- getTeam(inflMat)
    if (twoTeams) {
        groupsData[[1]] <- group[[1]]
        inflMat <- inflMat[group[[2]], group[[2]]]
        group <- getTeam(inflMat)
        groupsData[[2]] <- group[[1]]
        stray <- group[[2]]
        l <- groupsData %>% sapply(paste0, collapse = ",") %>% sort
        writeLines(l, paste0(net, ".BruteTeams"))
        return(list(groupsData, stray))

    }
    while(length(group[[1]])>1) {
        nGroups <- nGroups + 1
        groupsData[[ngroups]] <- group[[1]]
        inflMat <- inflMat[group[[2]], group[[2]]]
        group <- getTeam(inflMat)
    }
    stray <- unlist(group)
    l <- groupsData %>% sapply(paste0, collapse = ",") %>% sort
    writeLines(l, paste0(net, ".BruteTeams"))
    return(list(groupsData, stray))
}


