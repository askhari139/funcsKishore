loopDatGen <- function(folder){
    setwd(folder)
    command <- paste0("python ../signConsistency.py undirectedLoops")
    system(command)
    command <- paste0("python ../Positive_loops.py directedLoops")
    system(command)
    setwd("..")
}

consistencyFormatting <- function(net, edges, edgeNames)
{
    cycles <- read_csv(paste0(net, "_undirectedLoops.csv"), lazy = F) %>% select(Cycles) %>% unlist
    cycleDat <- lapply(cycles, function(cyc){
        nodes <- str_split(cyc, ",")[[1]]
        if(length(nodes) == 1)
        {
            ed <- data.frame(Edge1 = paste0(nodes, ",", nodes))
        }
        else
        {
            nodes <- c(nodes, nodes[1])
            ed <- lapply(1:(length(nodes)-1), function(i){
                e <- c(nodes[i], nodes[i+1])
                e <- list(e, rev(e)) %>% sapply(paste0, collapse = ",")
                e <- e[which(e %in% edgeNames)]
            }) %>% expand.grid %>% set_names(paste0("Edge", 1:ncol(.)))
        }
        
        ed %>% mutate(Nature = ed %>% apply(1, function(e){edges[e] %>% prod})) %>%
            mutate(Nature = ifelse(Nature == 1, "P", "N")) %>% 
            unite("EdgeList",contains("Edge"), sep = ";") %>%
            mutate(Cycles = cyc, Edge_count = str_count(EdgeList, ";") + 1)
    }) %>% reduce(rbind.data.frame)
    cycleDat
}

consistency <- function(net, compute = F)
{
    loopFile <- paste0(net, "_Consistency.csv")
    if (!file.exists(loopFile) || compute)
    {
        topoFiles <- list.files(".", ".topo")
        nets <- topoFiles %>% str_remove(".topo")
        loopData <- sapply(nets, function(x){
            topoDf <- read.delim(paste0(x, ".topo"), sep = "") %>% 
                unite(col = "Edge",sep = ",", Source, Target) %>% 
                mutate(Type = ifelse(Type == 2, -1, 1))
            edges <- topoDf$Type
            names(edges) <- topoDf$Edge
            edgeNames <- names(edges)
            setwd("undirectedLoops")
            df <- consistencyFormatting(x, edges, edgeNames)
            setwd("..")
            # browser()
            # d <- read.csv(paste0("undirectedLoops/",x, "_undirectedLoops.csv"))
            posLoops <- c(sum(df$Nature == "P"), 
                          df %>% filter(Nature == "P") %>% mutate(Edge_count = 1/Edge_count) %>% 
                              select(Edge_count) %>% unlist %>% sum)
            negLoops <- c(sum(df$Nature == "N"), 
                          df %>% filter(Nature == "N") %>% mutate(Edge_count = 1/Edge_count) %>% 
                              select(Edge_count) %>% unlist %>% sum)
            fracNeg <- negLoops[1]/(negLoops[1] + posLoops[1])
            print(x)
            write.csv(df, paste0("undirectedLoops/",x, "_undirectedLoopsF.csv"), row.names = F)
            df$ID <- 1:nrow(df)
            edgeDat <- lapply(edgeNames, function(e){
                df %>% filter(str_detect(EdgeList, e)) %>% group_by(Nature) %>%
                    summarise(Count = n(), loops = paste0(ID, collapse = ",")) %>% 
                    mutate(Edge = e)
            }) %>% reduce(rbind.data.frame) %>% 
                mutate(Nature = factor(Nature, levels = c("P", "N"))) %>% 
                spread(key = Nature, value = Count)
            if (!("P" %in% colnames(edgeDat)))
                edgeDat$P <- 0
            if (!("N" %in% colnames(edgeDat)))
                edgeDat$N <- 0
            edgeDat[is.na(edgeDat)] <- 0
            negDat <- edgeDat %>% filter(P == 0) %>% arrange(-N)
            posDat <- edgeDat %>% filter(N == 0)
            mixDat <- edgeDat %>% filter(P != 0, N != 0)
            nLoops <- df$ID[df$Nature == "N"]
            coveredLoops <- c()
            swapEdges <- 0
            if(nrow(negDat) != 0)
            {
                dummy <- sapply(negDat$loops, function(l){
                    l <- l %>% str_split(",") %>% unlist %>% as.integer
                    i <- intersect(coveredLoops, l)
                    if (is.null(i) || length(i) < length(l))
                    {
                        coveredLoops <<- c(coveredLoops, l) %>% unique
                        swapEdges <<- swapEdges + 1
                    }
                })
                if (length(coveredLoops) < length(nLoops))
                {
                    remaining <- nLoops[!(nLoops %in% coveredLoops)]
                    dat <- edgeDat %>% 
                        filter(sapply(loops, function(l){
                            any(str_detect(l, as.character(remaining)))
                        }))
                    loops <- paste0(dat$loops, collapse = ",") %>% str_split(",") %>%
                        unlist %>% as.integer
                    pLoops <- loops[!(loops %in% remaining)]
                    dummy <- sapply(dat$loops, function(l){
                        l <- l %>% str_split(",") %>% unlist %>% as.integer
                        l <- l[!(l %in% pLoops)]
                        i <- intersect(coveredLoops, l)
                        if (is.null(i) || length(i) < length(l))
                        {
                            coveredLoops <<- c(coveredLoops, l) %>% unique
                            swapEdges <<- swapEdges + 1
                            return(T)
                        }
                        return(F)
                    })
                    ed <- dat$Edge[dummy]
                    dat <- edgeDat %>% 
                        filter(sapply(loops, function(l){
                            any(str_detect(l, as.character(pLoops)))
                        })) %>% filter(!(Edge %in% ed))
                    coveredLoops <- c()
                    dummy <- sapply(dat$loops, function(l){
                        l <- l %>% str_split(",") %>% unlist %>% as.integer
                        i <- intersect(coveredLoops, l)
                        if (is.null(i) || length(i) < length(l))
                        {
                            coveredLoops <<- c(coveredLoops, l) %>% unique
                            swapEdges <<- swapEdges + 1
                            return(T)
                        }
                        return(F)
                    })
                    
                }
            }
            c(posLoops, negLoops, swapEdges)
        }) %>% t %>% data.frame %>% 
            set_names(c("posLoopsUnd", "posLoopsUndWeighted","negLoopsUnd", "negLoopsUndWeighted",
                        "inconsistency")) %>%
            mutate(Net = nets)
        write.csv(loopData, paste0(net, "_Consistency.csv"), row.names = F)
    }
    read.csv(paste0(net, "_Consistency.csv")) %>% mutate(Network = Net) %>% 
        select(-Net)
}

