
randMat <- function(mat, num, impurity = 0.0, density = 0.3)
{
    n <- length(mat)
    k <- round(n*density)
    id <- sample(1:n, k)
    ed <- sign(num)
    if (impurity != 0) {
        idImpure <- sample(id, round(k*impurity))
    }
    edOp <- ed*-1
    mat[id] <- ed
    if (impurity != 0) {
        mat[idImpure] <- edOp
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

# checkInteger <- function(x) {
#     if (floor(x) == x)
# }

setGroupMatrix <- function(groups,
                           density,
                           impurity,
                           densityMatrix = NULL,
                           impurityMatrix = NULL,
                           orderMatrix = NULL) {
    nGroups <- length(groups)
    nodes <- unlist(groups) %>% sort
    nNodes <- length(nodes)
    if (is.null(densityMatrix)) {
        densityMatrix <- rep(density, nGroups*nGroups) %>% 
            matrix(ncol = nGroups)
    }
    if (is.null(impurityMatrix)) {
        impurityMatrix <- rep(impurity, nGroups*nGroups) %>% 
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
            m <- randMat(intMat[g1, g2], orderMatrix[i,j],impurityMatrix[i,j], 
                    densityMatrix[i,j])
        }) %>% reduce(cbind)
    }) %>% reduce(rbind)
    intMat
}

netGen <- function(groupSizes, numNodes = 0, orderMatrix = NULL, 
                   densityMatrix = NULL,density = 0.3, 
                   impurityMatrix = NULL, impurity = 0.0, 
                   id = sample(1:100, 1), key = "ArtNet")
{
    intCheck <- all(floor(groupSizes) == groupSizes)
    if ((numNodes == 0) & !intCheck) {
        stop("Either provide an integer vector of group sizes or a number of nodes")
    }
    else if (intCheck & (numNodes == 0)) {
        numNodes <- sum(groupSizes)
    }
    else {
        if ((sum(groupSizes) != 1) & (sum(groupSizes) != numNodes)) {
            groupSizes <- groupSizes/sum(groupSizes)
        }
        groupSizes <- round(numNodes*groupSizes)
        groupSizes[length(groupSizes)] <- numNodes - sum(groupSizes[-length(groupSizes)])
    }
    names(groupSizes) <- LETTERS[1:length(groupSizes)]
    groups <- lapply(LETTERS[1:length(groupSizes)], function(x) {
        paste0(x, 1:groupSizes[x])
    })
    nodes <- unlist(groups)
    intMat <- setGroupMatrix(groups, density, impurity, densityMatrix, 
        impurityMatrix, orderMatrix)
    rownames(intMat) <- colnames(intMat) <- nodes
    df <- data.frame(intMat) %>% mutate(Source = rownames(intMat)) %>%
        gather(key = "Target", value = "Type", -Source) %>%
        filter(Type != 0) %>% mutate(Type = ifelse(Type == 1, 1, 2))
    netName <- paste0(key, "_", id, "_", length(groupSizes), ".topo")
    write_delim(df, netName, " ")
    groupLines <- groups %>% sapply(function(x) {
        x %>% paste0(collapse = ",")
    })
    writeLines(groupLines, 
               paste0(key, "_", id, "_", length(groupSizes), ".AssignedTeams"))
}
netGen <- cmpfun(netGen)
