netx <- import("networkx")

cleanTopo <- function(topoFile, delim = "") { # remove un-necessary hiphens, use single space delimiter
    df <- read.delim(topoFile, sep = delim) %>%
        mutate(Source = Source %>% str_replace_all(regex("\\W+"), ""),
               Target = Target %>% str_replace_all(regex("\\W+"), ""))
    write_delim(df, topoFile, delim = " ", quote = "none")
    
}

getPeripheral <- function(topoFile) {
    ls <- TopoToIntMat(topoFile)
    intmat <- ls[[1]]
    nodes <- ls[[2]]
    signal <- which(apply(intmat, 2, function(x) {
        all(x == 0)
    }))
    output <- which(apply(intmat, 1, function(x) {
        all(x == 0)
    }))
    secondary_signal <- which(apply(intmat, 2, function(x) {
        if (length(signal) != 0) {
            all(x[-signal] == 0)
        } else {
            F
        }
    }))
    secondary_output <- which(apply(intmat, 1, function(x) {
        if (length(output) != 0) {
            all(x[-output] == 0)
        } else {
            F
        }
    }))
    nonEssentials <-
        c(signal, output, secondary_signal, secondary_output) %>% unique
    nodes[-nonEssentials]
}

SecondarySignals <- function(topoDf, sig) {
    targets <- topoDf %>% filter(Source == sig) %>%
        select(Target) %>% unlist
    secSigs <- c()

}


mergeNets <- function(nets, nam) {
    setwd(topoFileFolder)
    df <- lapply(nets, function(net) {
        read.delim(paste0(net, ".topo"), sep = "") %>%
            mutate(Source = tolower(Source),
                   Target = tolower(Target))
    }) %>% reduce(rbind.data.frame) %>% distinct
    write_delim(df,
                paste0(nam, ".topo"),
                delim = " ",
                quote = "none")
}

getImpurities <- function(topoFile){
    net <- topoFile %>% str_remove(".topo")
    ls <- TopoToIntMat(topoFile)
    intMat <- ls[[1]]
    NodesDef <- ls[[2]]
    rownames(intMat) <- colnames(intMat) <- NodesDef
    teams <- readLines(paste0(net, ".teams")) %>% str_split(",")
    nTeams <- length(teams)
    lengths <- sapply(teams, length)
    impurities <- sapply(1:nTeams, function(i) {
        sapply(1:nTeams, function(j) {
            nMax <- lengths[i]*lengths[j]
            subMat <- intMat[teams[[i]], teams[[j]]]

            if (i == j) {
                imp <- sum(subMat == -1)
            }
            else {
                imp <- sum(subMat == 1)
            }

        })
    })
}

getLoopData <- function(topoFile, size_limit = NULL, netx = netx) {
    topoDf <- read_delim(topoFile, delim = " ", show_col_types = F) %>% 
        mutate(Sign = ifelse(Type == 2, -1, Type), 
        Edges = paste0(Source, "_", Target))
    g <- netx$from_pandas_edgelist(topoDf, source = 'Source', target = 'Target', 
        edge_attr = "Sign", create_using = netx$DiGraph())
    # g <- netx$DiGraph()

    cycles <- g %>%
        netx$simple_cycles(length_bound = size_limit) %>%
        pyBuiltins$list()
    loopSigns <- sapply(cycles, function(cycle) {
        cycle <- c(cycle, cycle[1])
        edges <- data.frame(Source = cycle[1:(length(cycle)-1)], 
            Target = cycle[2:length(cycle)]) %>%
            mutate(Edges = paste0(Source, "_", Target))
        edgeSigns <- topoDf %>% 
            filter(Edges %in% edges$Edges) %>%
            pull(Sign)
        ifelse(prod(edgeSigns) == -1, "N", "P")
    })
    loopLengths <- sapply(cycles, length)
    cycleSeq <- sapply(cycles, paste0, collapse = "_")
    df <- data.frame(Cycles = cycleSeq, Nature = loopSigns, Edge_count = loopLengths)
    write_csv(df, str_replace(topoFile, ".topo", "_directedLoops.csv"), quote = "none")
}