### my algo
## Idea is as follows:
# count the number of possible unique combinations / input the number of combinations required (m)
# count the number of 1's and 2's and identify the lesser count (n)
# define unique perm counter
# while the counter is less than m, sample n positions, put the corresponding digits in it,
#  convert the positions to binary, store it if the combination is unique.

uniquePerms <- function(x, max = 0){
    N <- length(x)
    n1 <- sum(x == 1)
    n2 <- sum(x == 2)
    nVec <- c(n1,1)
    if (n1>n2) nVec <- c(n2,2)
    notOne <- ifelse(nVec[2] == 1, 2, 1)
    if (max == 0)
        max <- factorial(N)/(factorial(n1)*factorial(n2))
    if (max > 2000)
        max <- 2000
    perms <- rep(0,max+1)
    perms[1] <- sum(2^(which(x == nVec[2])))
    lists <- matrix(rep(notOne, max*N), ncol = N)
    count <- 1
    while(count<max){
        nSamp <- sample(1:N, nVec[1])
        perm <- sum(2^nSamp)
        if (!any(perms == perm))
        {
            perms[count+1] <- perm
            lists[count,nSamp] <- nVec[2]
            count <- count + 1
            print(perm)
        }
    }
    return(lists)
}

randNetGen <- function(topo_files, maxNets = 500) {
    sapply(topo_files, function(file) {
        name <- str_remove(file, ".topo")
        dir.create(name)
        df <- read.delim(file, sep = " ")
        setwd(name)
        wt <- df$Type
        onetwo <- df[, 1:2]
        rand_orders <- uniquePerms(wt, max = maxNets)
        dummy <- sapply(1:nrow(rand_orders), function(x){
            y <- rand_orders[x,]
            df1 <- cbind.data.frame(onetwo, y)
            colnames(df1) <- c("Source", "Target", "Type")
            write_delim(df1, paste0(name, "_", x, ".topo"), delim = " ")
        })
        setwd("..")
    })
}
# topo_files <- list.files("topo_files", pattern = ".topo")
#
# setwd("topo_files")
