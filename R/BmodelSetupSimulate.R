BmodelSetup <- function() {
    x <- system("julia -e \"using Boolean\"", ignore.stderr = TRUE)
    if (x == 1)
        system("julia -e \" using Pkg; Pkg.add(url = \"https://github.com/askhari139/Boolean.jl\")\"")
}

simulation <- function(setupBmodel = F, simPackage = "../Bmodel", numThreads = 3, nInit = 100000) {
    nThreads <- parallel::detectCores()
    if(numThreads > 0.8*nThreads) {
        numThreads <- round(0.8*nThreads)
    }
    if (setupBmodel) {
        BmodelSetup()
    }
    wd <- getwd()
    setwd(simPackage)
    simPackage <- getwd()
    setwd(wd)
    topoFiles <- list.files(".", ".topo$")
    size <- floor(length(topoFiles) / numThreads)
    topoList <- lapply(1:numThreads, function(x) {
        k <- (x - 1) * size
        id <- 1:size + k
        if (size + k > length(topoFiles))
            id <- id[1]:(length(topoFiles))
        topoFiles[id] %>% paste0(collapse = " ")
    })
    future::plan(future::multisession, workers = numThreads)
    simulater <- future.apply::future_lapply(topoList, function(x) {
        topoFiles <- paste0("[\"", paste0(x, collapse = "\",\""), "\"]")
        command <- paste0("julia -e 'using Boolean; topoFiles = ", topoFiles, "; for topoFile in topoFiles; y = @elapsed x = bmodel_reps(topoFile; nInit = ",nInit,"); end'")
        system(command)
    })
    future:::ClusterRegistry("stop")
}

simulateRACIPE <- function(racipePath, gkNorm = F, multiThread = T, numThreads = 8) {
    topoFiles <- list.files(".", ".topo")
    DirectoryNav("RACIPE")
    sapply(topoFiles, function(x) {

        net <- str_remove(x, ".topo")
        if(file.exists(paste0(net, "_C_solution.dat")))
            return()
        file.copy(paste0("../", x), paste0(net, "_A.topo"))
        file.copy(paste0("../", x), paste0(net, "_B.topo"))
        file.copy(paste0("../", x), paste0(net, "_C.topo"))
        tpFls <- str_replace(x, ".topo", paste0("_", c("A", "B", "C"), ".topo"))
        cmds <- sapply(tpFls, function(t) {
            paste0(racipePath, " ", t, " -num_paras 10000",
                ifelse(multiThread, paste0(" -threads ", numThreads), ""),
                " > ", str_replace(t, ".topo", ".log")) %>% system()
            prs <- str_replace(t, ".topo", ".prs") %>% read.delim
            nodes <- prs %>% filter(str_detect(Parameter, "Prod")) %>%
                select(Parameter) %>% unlist %>% str_remove("Prod_of_")
            solnNames <- c("parIndex", "nStates", "basin", nodes)
            solnDf <- read_delim(str_replace(t, ".topo", "_solution.dat"),
                col_types = "d", col_names = solnNames, delim = "\t")
            parNames <- prs %>% select(Parameter) %>%
                    unlist %>% c("parIndex", "nStates", .)
            parDf <- read_delim(str_replace(t, ".topo", "_parameters.dat"),
                col_types = "d", col_names = parNames, delim = "\t")
            write_delim(parDf, str_replace(t, ".topo", "_parameters.dat"),
                delim = "\t", quote = "none")
            if (gkNorm) {
                p <- sapply(nodes, function(node) {
                    parDf %>% select(all_of(paste0(c("Prod_of_", "Deg_of_"), node))) %>%
                        set_names(c("p", "d")) %>%
                        mutate(r = p/d) %>%
                        select(r) %>% unlist %>% log2
                }) %>% data.frame %>% set_names(paste0(c("Norm_"), nodes)) %>%
                    mutate(parIndex = parDf$parIndex)
                solnDf <- solnDf %>%
                    merge(p, by = "parIndex", all = T)
                sapply(nodes, function(node) {
                    solnDf[[node]] <<- solnDf[[node]] - solnDf[[paste0("Norm_", node)]]
                })
            }
            write_delim(solnDf, str_replace(t, ".topo", "_solution.dat"),
                delim= "\t", quote = "none")
        })
        list.files(".", "solution_\\d+.dat") %>% sapply(file.remove)
        list.files(".", ".cfg") %>% sapply(file.remove)
        list.files(".", "T_test") %>% sapply(file.remove)
    })

    setwd("..")
}

simulateNetworkBmodel <- function(topoFile, simpackage = "./Bmodel", shubham = "false",
        nInit = "100000", nIter = "1000", mode = "Async", stateRep = -1,
        randSim = "false", discrete = "false", nLevels = 2, vaibhav = "false") {
    wd <- getwd()
    if (!file.exists(topoFile)) {
        print("topoFile not found!")
        return()
    }
    setwd(simpackage)
    simpackage <- getwd()
    setwd(wd)
    script <- c()
    script <- c(script, paste0("include(\"", simpackage, "/bmodel.jl\")"))
    # script <- c(script,
    #     paste0('y1 = @elapsed x = bmodel_reps("',
    #         topoFile, "\"",
    #         "; nInit = ", nInit, ", nIter = ", nIter,
    #         ", mode = \"", mode, "\", stateRep = ", stateRep,
    #         ", randSim = ", randSim, ", shubham = true, discrete = ",
    #         discrete,", nLevels= ", nLevels,  ")"))
    script <- c(script,
        paste0("y1 = @elapsed x = bmodel_reps(\"",
            topoFile, "\"",
            "; nInit = ", nInit, ", nIter = ", nIter,
            ", mode = \"", mode, "\", stateRep = ", stateRep,
            ", randSim = ", randSim, ", shubham = ", shubham, ", discrete = ",
            discrete, ", nLevels = ", nLevels, ", vaibhav = ", vaibhav, ")"))
    script <- c(script,
        paste0("print(\"", topoFile, "\", \" - \", y1, \" seconds.\")"))
    writeLines(script, "dummy.jl")
    system("julia dummy.jl")
}
