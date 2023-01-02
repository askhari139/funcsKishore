BmodelSetup <- function(simPackage = "../Bmodel") {
    wd <- getwd()
    DirectoryNav(simPackage)
    simPackage <- getwd()
    download.file(
        "https://github.com/askhari139/Boolean.jl/archive/refs/heads/main.zip",
        destfile = "Bmodel.zip"
    )
    zipFile <- "Bmodel.zip"
    unzip(zipFile, junkpaths = T)
    file.remove(zipFile)
    script <- readLines("script.jl")
    script[1] <- paste0("include(\"", simPackage, "/bmodel.jl\")")
    writeLines(script, "script.jl")
    
    script <- readLines("scriptWindows.jl")
    script[1] <- paste0("include(\"", simPackage, "/bmodel.jl\")")
    writeLines(script, "scriptWindows.jl")
    JuliaCall::julia_source("dependencyInstaller.jl")
    setwd(wd)
}

simulation <- function(setupBmodel = F, simPackage = "../Bmodel", numThreads = 3) {
    nThreads <- parallel::detectCores()
    if(numThreads > 0.8*nThreads) {
        numThreads <- round(0.8*nThreads)
    }
    if (setupBmodel) {
        if (!file.exists(paste0(simPackage, "/bmodel.jl")))
            BmodelSetup(simPackage)
    }
    if (!dir.exists(simPackage)|| !file.exists(paste0(simPackage, "/bmodel.jl"))) {
        message("Bmodel package not found. Run with setupBmodel = T if you would like to set it up")
        return()
    }
    wd <- getwd()
    setwd(simPackage)
    simPackage <- getwd()
    setwd(wd)
    os <- .Platform$OS.type
    script <- "script.jl"
    if (os == "windows")
        script <- "scriptWindows.jl"
    file.copy(paste0(simPackage, "/", script), ".", overwrite = T)
    if (os != "windows") {
        Sys.setenv(JULIA_NUM_THREADS = as.character(numThreads))
        command <- "julia script.jl"
        system(command)
    }
    else {
        topoFiles <- list.files(".", ".topo$")
        size <- floor(length(topoFiles) / numThreads)
        topoList <- lapply(1:numThreads, function(x) {
            k <- (x - 1) * size
            id <- 1:size + k
            if (size + k > length(topoFiles))
                id <- id[1]:(length(topoFiles))
            topoFiles[id] %>% paste0(collapse = " ")
        })
        Sys.setenv(JULIA_NUM_THREADS = as.character(1))
        future::plan(future::multisession, workers = numThreads)
        simulater <- future.apply::future_lapply(topoList, function(x) {
            command <- paste0("julia ", script, " ", x)
            system(command)
        })
        future:::ClusterRegistry("stop")
    }
}
