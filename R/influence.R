#' Title
#'
#' @param topoFile
#' @param plotOut
#' @param nOrder
#'
#' @return
#' @export
#'
#' @examples
TopoToIntMat <-
    function(topoFile,
             plotOut = F,
             nOrder = NULL) {
        df <- read.delim(topoFile, sep = " ", stringsAsFactors = F)
        df <- df %>%
            mutate(Type = ifelse(Type == 2,-1, 1))

        nodes <- unique(c(df$Source, df$Target)) %>%
            sort(decreasing = T)
        n_nodes <- length(nodes)
        intmat <- rep(0, n_nodes * n_nodes) %>%
            matrix(ncol = n_nodes)
        df1 <- df %>%
            mutate(Source = sapply(Source, function(x) {
                which(nodes == x)
            })) %>%
            mutate(Target = sapply(Target, function(x) {
                which(nodes == x)
            }))
        dummy <- apply(df1, 1, function(x) {
            i <- x[1]
            j <- x[2]
            k <- x[3]
            intmat[i, j] <<- k
        })
        if (plotOut) {
            if (is.null(nOrder)) {
                nOrder <- getEMSONodes(topoFile)
            }
            df <- data.frame(intmat) %>%
                set_names(nodes) %>%
                mutate(nodes1 = nodes) %>%
                gather(key = "Nodes", value = "Influence",-nodes1) %>%
                mutate(
                    nodes1 = factor(nodes1, levels = nOrder),
                    Nodes = factor(Nodes, levels = nOrder)
                )
            ggplot(df, aes(
                x = Nodes,
                y = nodes1,
                fill = Influence
            )) + geom_tile() +
                theme_Publication() + scale_fill_gradient2(low = "blue",
                                                           high = "red",
                                                           limits = c(-1, 1)) +
                labs(fill = "Edge") +
                theme(
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.x = element_text(angle = 90),
                    legend.position = "right",
                    legend.direction = "vertical",
                    legend.key.height = unit(0.5, "in")
                )
            DirectoryNav("MatrixPlots")
            ggsave(paste0(str_replace(
                topoFile, ".topo", "_interaction.png"
            )),
            width = 7,
            height = 6)
            setwd("..")
        }
        return(list(intmat, nodes))
    }

#' Title
#'
#' @param mat
#' @param power
#'
#' @return
#' @export
#'
#' @examples
ComputePowerMatrix <- function(mat, power) {
    res <- mat
    if (power == 1) {
        return(res)
    }
    for (i in 2:power) {
        res <- res %*% mat
    }
    return(res)
}
ComputePowerMatrix <- cmpfun(ComputePowerMatrix)

#' Title
#'
#' @param net
#' @param intmat
#' @param nodes
#' @param lmax
#' @param write
#'
#' @return
#' @export
#'
#' @examples
InfluenceMatrix <-
    function(topoFile,
             lmax = 10,
             write = T) {
        if (file.exists(paste0("Influence/", net, "_reducedInfl.csv"))) {
            influence_reduced <-
                read_csv(paste0("Influence/", net, "_reducedInfl.csv"),
                         col_types = cols())
            rn <- influence_reduced[, 1] %>% unlist
            influence_reduced <- influence_reduced[,-1] %>% as.matrix
            rownames(influence_reduced) <- rn
            return(influence_reduced)
        }

        net <- topoFile %>% str_remove(".topo")
        ls <- TopoToIntMat(topoFile)
        intmat <- ls[[1]]
        nodes <- ls[[2]]
        intmax <- intmat
        intmax[which(intmax == -1)] <- 1
        res <- 0
        for (l in 1:lmax) {
            intM <- ComputePowerMatrix(intmat, l)
            maxM <- ComputePowerMatrix(intmax, l)
            r1 <- intM / maxM
            r1[is.nan(r1)] <- intM[is.nan(r1)]
            res <- res + r1
        }
        res <- res / lmax

        nodes <- nodes

        influence_mat <- res
        colnames(influence_mat) <- rownames(influence_mat) <- nodes
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
        if (length(nonEssentials)) {
            influence_reduced <- influence_mat[-nonEssentials,-nonEssentials]
            nodes_reduced <- nodes[-nonEssentials] %>%
                c() %>%
                str_replace_all(regex("\\W+"), "")
        } else {
            influence_reduced <- influence_mat
            nodes_reduced <- nodes %>% str_replace_all(regex("\\W+"), "")
        }
        if (length(nodes_reduced) < 2) {
            return()
        }
        rownames(influence_reduced) <-
            colnames(influence_reduced) <- nodes_reduced
        if (write) {
            DirectoryNav("Influence")
            write.csv(influence_mat, paste0(net, "_fullInfl.csv"))
            write.csv(influence_reduced, paste0(net, "_reducedInfl.csv"))
            setwd("..")
        }
        influence_reduced
    }
InfluenceMatrix <- cmpfun(InfluenceMatrix)