formatRACIPE <- function(topoFile, force = F, gkNorm = T) {
    # setwd("RACIPE")
    #### Setup ------
    net <- topoFile %>% str_remove(".topo")
    paramNames <- read.delim(paste0(net, ".prs"), sep = "") %>%
        select(Parameter) %>%
        unlist
    nodes <- paramNames[str_detect(paramNames, "Prod_of")] %>%
        str_remove("Prod_of_")
    if (file.exists(paste0(net, "_format.csv")) && !force) {
        stop("Previous Formatting exists. Use force=T if recalculation is needed")
        dat <- read_csv(paste0(net, "_format.csv"), show_col_types = F)
        normFactor <- read_csv(paste0(net, "_gkNormFactor.csv"), show_col_types = F)
    }
    ### Formatting ------
    else{
        topoDf <- read.delim(topoFile, sep = "", stringsAsFactors = F)
        deletables <- paste0(net, c("_solution_", "_T_test_", ".cfg"))
        fileList <- sapply(deletables, function(x) {
            fAll <- list.files(pattern = x)
            sapply(fAll, file.remove)
        })
        if(!file.exists(paste0(net, "_solution.dat"))) {
            print("Solution file does not exits. Simulate!")
            # setwd("..")
            return(NA)
        }
        dat <- read_delim(paste0(net, "_solution.dat"), 
            delim = "\t", col_names = F, col_types = "d")
        colnames(dat) <- c("paramID", "nStates", "Basin", nodes)
        parameters <- read_delim(paste0(net, "_parameters.dat"),
            delim = "\t", col_names = F, col_types = "d") %>%
            set_names(c("paramID", "nStates", paramNames))
        normFactor <- (parameters %>%
            select(all_of(paste0("Prod_of_", nodes))) %>% set_names(nodes))/(
                parameters %>% select(all_of(paste0("Deg_of_", nodes))) %>% set_names(nodes)) 
        normFactor <- normFactor %>%
            mutate(across(everything(), log2)) %>%
            mutate(paramID = parameters$paramID) %>%
            merge(dat %>% select(paramID), by = "paramID")
        
        write_csv(normFactor, paste0(net, "_gkNormFactor.csv"))
        write_csv(dat, paste0(net, "_formatRacipe.csv"))
    }
    if (gkNorm) {
        sapply(nodes, function(x) {
            dat[[x]] <<- dat[[x]] - normFactor[[x]]
        })
    }
    return(dat)
}

zNormDf <- function(df, nodes) {
    df %>% select(all_of(nodes)) %>% mutate(across(everything(), .fns = function(x){
        ifelse(x-mean(x) >0, 1, 0)
        }))
}

minMaxDf <- function(df, nodes, nLevels=100) {
    df %>% select(all_of(nodes)) %>% mutate(across(everything(), .fns = function(x){
        x <- x-min(x)/(max(x) - min(x))
        if (is.na(nLevels) || !is.numeric(nLevels))
            nLevels <- 100
        
        y <- round(x*nLevels)/nLevels
        return(y)
        }))

}

discretize <- function(topoFile, method = c("zScore", "minMax"), 
    nLevels = NULL, gkNorm = F) {
        net <- topoFile %>% str_remove(".topo")
        formatData <- formatRACIPE(topoFile)
        gkNormFactor <- read_csv(paste0(net, "_gkNormFactor.csv"), show_col_types = F)
        method <- method[1]
        if (method == "zScore")
            f <- zNormDf
        else if (method == "minMax") {
           f <- minMaxDf
        }
        else {
            stop("Wrong discretization method")
        }
    }