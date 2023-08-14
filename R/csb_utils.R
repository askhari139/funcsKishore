idFileGen <- function(topo) {
    cleanTopo(topo)
    df <- read_delim(topo, delim = " ")
    nodes <- c(df$Source, df$Target) %>% unique %>% sort
    df <- data.frame(node = nodes, id = (1:length(nodes))-1)
    write_delim(df, str_replace(topo, ".topo", ".ids"), delim = " ", quote = "none")
}

getInit <- function(topo, paramAdditional = NULL, valueAdditional = NULL) {
    initFormat <- data.frame(parameter = c("input_folder_name", "otput_folder_name", 
        "input_filenames", "num_runs", "num_simulations", "maxtime", "constant_node_count",
        "time_step", "updation_rule"), 
        value = c(".", ".", str_remove(topo, ".topo"), 3, 10000, 1000, 0, 0.1, 2))
    if (!is.null(paramAdditional) && !is.null(valueAdditional) && length(paramAdditional) == length(valueAdditional)) {
        ids <- paramAdditional %in% initFormat$parameter
        if(any(ids))
        initFormat <- initFormat %>% 
            filter(!parameter %in% paramAdditional[ids]) %>%
            bind_rows(data.frame(parameter = paramAdditional[ids], value = valueAdditional[ids]))
    }
    write_delim(initFormat, "init.txt", delim = " ", quote = "none")
}

processCsbOutput <- function(topo, nodeOrder = NULL) {
    net <- topo %>% str_remove(".topo")
    ssFils <- list.files(".", paste0(net, "_ss_run"))
    if (length(ssFils) == 0) {
        print("No steady state files found!")
        return()
    }
    dfAll <- lapply(ssFils, function(x) {
        df <- read_delim(x, delim = " ", col_types = "d") %>%
            select(-ID)
        if (!is.null(nodeOrder)) {
            df <- df %>% select(all_of(nodeOrder))
        }
        df <- df %>%
            mutate(across(everything(), function(x) as.character(ifelse(x == -1, 0, x)))) %>%
            unite("State", everything(), sep = "") %>%
            group_by(State) %>%
            summarise(Count = n()) %>%
            mutate(Prob = Count/sum(Count)) %>%
            select(State, Prob)
    }) %>%
        reduce(full_join, by = "State") %>%
        mutate(across(is.numeric, function(x) ifelse(is.na(x), 0, x))) %>%
        mutate(Avg = rowMeans(select(., -State)),
                Std = select(., -State) %>% apply(1, sd)) %>%
        select(State, Avg, Std)
    write_delim(dfAll, paste0(net, "_ss.csv"), delim = ",", quote = "none")
}