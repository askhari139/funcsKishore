library(tidyverse)
library(compiler)


randMat <- function(mat, num, sparsity = 0.3)
{
    n <- length(mat)
    k <- round(n*sparsity)
    id <- sample(1:n, k)
    ed <- sign(num)
    if (floor(num) != num) {
        edOp <- ed*-1
        pure <- round(length(id)*abs(num))
        idPure <- sample(id, pure)
        idImpure <- id[!(id %in% idPure)]
        mat[idPure] <- ed
        mat[idImpure] <- edOp
    }
    else {
        mat[id] <- ed
    }
    
    mat
}
randMat <- cmpfun(randMat)

breaker <- function(vec, parts) {
    if(sum(parts) != length(vec)) {
        parts[length(parts)] <- length(vec) - sum(parts[-length(parts)])
    }
    i <- 1
    groups <- lapply(parts, function(x) {
        i0 <- i
        i <<- i + x
        vec[i0:(i-1)]
    })
    groups
}

setGroupMatrix <- function(groups, 
                           sparsity,
                           sparsityMatrix = NULL, 
                           orderMatrix = NULL) {
    nGroups <- length(groups)
    nodes <- unlist(groups) %>% sort
    nNodes <- length(nodes)
    if (is.null(sparsityMatrix)) {
        sparsityMatrix <- rep(sparsity, nGroups*nGroups) %>% 
            matrix(ncol = nGroups)
    }
    if (is.null(orderMatrix)) {
        orderMatrix <- sapply(1:nGroups, function(i) {
            y <- rep(-1, nGroups)
            y[i] <- 1
            y
        })
    }
    intMat <- rep(0, nNodes*nNodes) %>% matrix(ncol = nNodes)
    rownames(intMat) <- colnames(intMat) <- nodes
    intMat <- lapply(1:nGroups, function(i) {
        lapply(1:nGroups, function(j) {
            g1 <- groups[[i]]
            g2 <- groups[[j]]
            randMat(intMat[g1, g2], orderMatrix[i,j],
                    sparsityMatrix[i,j])
        }) %>% reduce(cbind)
    }) %>% reduce(rbind)
    intMat
}

netGen <- function(numNodes, groupSizes, orderMatrix = NULL, sparsityMatrix = NULL,
                   sparsity = 0.3, id = sample(1:100, 1))
{
    nodes <- paste0("A", 1:numNodes)
    groups <- breaker(nodes, groupSizes)
    intMat <- setGroupMatrix(groups, sparsity, sparsityMatrix, orderMatrix)
    rownames(intMat) <- colnames(intMat) <- nodes
    df <- data.frame(intMat) %>% mutate(Source = rownames(intMat)) %>%
        gather(key = "Target", value = "Type", -Source) %>%
        filter(Type != 0) %>% mutate(Type = ifelse(Type == 1, 1, 2))
    netName <- paste0("ArtNet_", id, "_", length(groupSizes), ".topo")
    write_delim(df, netName, " ")
    groupLines <- groups %>% sapply(function(x) {
        x %>% paste0(collapse = ",")
    })
    writeLines(groupLines, 
               paste0("ArtNet_", id, "_", length(groupSizes), ".AssignedTeams"))
}
netGen <- cmpfun(netGen)