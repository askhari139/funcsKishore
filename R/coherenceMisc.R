retProb <- function(sInit, u)
{
    n <- length(sInit)
    f <- 1/n
    sapply(1:length(sInit), function(r) {
        s <- sInit
        s_dummy <- s%*%u %>% sign
        s[r] <- s_dummy[r]
        ifelse(all(s == sInit), f, 0)
    }) %>% sum
}
retProb <- cmpfun(retProb)

retProbPert <- function(sInit, u)
{
    n <- length(sInit)
    f <- 1/n
    sapply(1:length(sInit), function(r) {
        s <- sInit
        i <- sample(1:n, 3)
        s[i] <- -1*s[i]
        s_dummy <- s%*%u %>% sign
        s[r] <- s_dummy[r]
        ifelse(all(s == sInit), f, 0)
    }) %>% sum
}
retProbPert <- cmpfun(retProb)

Score <- function(s, id) {
    sum(replace(s, s == -1, 0)[id])/length(id)
}
Score <- cmpfun(Score)

potential <- function(topoFile, nStates = 100000) {
    net <- topoFile %>% str_remove(".topo")
    ls <- TopoToIntMat(topoFile)
    update_matrix <- ls[[1]]
    nodes <- ls[[2]]
    colnames(update_matrix) <- rownames(update_matrix) <- nodes
    nodeOrder <- readLines(paste0(net, "_nodes.txt"))
    update_matrix <- update_matrix[nodeOrder, nodeOrder]
    update_matrix <- 2*update_matrix + diag(length(nodes))
    groupLabels <- readLines(str_replace(topoFile, ".topo", ".teams")) %>% str_split(",")
    nNodes <- length(nodes)
    states <- sample(c(-1,1), nStates*nNodes, replace = T) %>%
        matrix(ncol = nStates, nrow = nNodes) %>% data.frame
    rownames(states) <- nodeOrder
    df <- data.frame(Escore = sapply(states, Score, id = which(nodeOrder %in% groupLabels[[1]])),
                     Mscore = sapply(states, Score, id = which(nodeOrder %in% groupLabels[[2]])),
                     prob = sapply(states, retProbPert, u = update_matrix)) %>%
        mutate(Potential = -log(prob - min(prob) + 1e-3)) %>%
        mutate(Potential = replace(Potential, Potential > median(Potential), 4))
        # mutate(Potential = 1-prob)
    write_csv(df %>% arrange(Escore, desc(Mscore)), paste0(net, "_landscapeData.csv"), quote = "none")
    ggplot(df %>% group_by(Escore,Mscore) %>% summarise(Potential=mean(Potential)),
           aes(x = Escore, y = Mscore, fill = Potential)) + geom_tile() +
        theme_Publication() + scale_fill_viridis_c() +
        labs(x = "Epithelial Score", y = "Mesenchymal Score", fill = "Pseudo\nPotential") +
        theme(legend.position = "right", legend.direction = "vertical", legend.key.height = unit(0.8, "cm"))
    ggsave(paste0(net, "_landScape.png"), width = 6, height = 5)
}

# potential("wild.topo")
# potential("EMT_RACIPE_102.topo")
