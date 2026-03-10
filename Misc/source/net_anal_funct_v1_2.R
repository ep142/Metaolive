# Network analysis functions v1.1 
# 1/11/2024
# these functions are adapted from MAN_fucntions.R
# see GitHub repo MAN_in_cheese

# microbial association network inference ---------------------------------
# infers microbial association network with NetCoMi::netConstruct using 
# error catching (returning try-error object if it fails, the microNet object if
# successful); the method and the parameter list are provided with the function 
# (and are set as options in the main script)

infer_MAN <- function(myphyseq, inf_method, method_parameters){
  mnet_object <- try(NetCoMi::netConstruct(myphyseq,
                                           measure = inf_method,
                                           normMethod = method_parameters$normmethodPar,
                                           zeroMethod = method_parameters$zeromethodPar,
                                           sparsMethod = method_parameters$sparMethod,
                                           measurePar = method_parameters$measureParList,
                                           alpha = method_parameters$alpha,
                                           dissFunc = method_parameters$dissFuncPar,
                                           verbose = method_parameters$verbosePar,
                                           cores = ncores,
                                           seed = 123456
  )
  )
  return(mnet_object)
}


# calculate network stats -------------------------------------------------
# calculates net stats if possible, otherwise returns a try-error object. 
# I am just using defaults for most parameters, 
# in the future I might want to set more parameters
# note that centralities are calculated for the LCC by default which
# is the most reasonable choice; if you want to change that, set lcc to FALSE

calculate_net_stats <- function(microNet_obj, verbosePar = 1, lcc = TRUE){
  if(class(microNet_obj) == "microNet"){
    netstats <- try(
      netAnalyze(
        microNet_obj,
        centrLCC = lcc,
        avDissIgnoreInf = FALSE,
        sPathAlgo = "dijkstra",
        sPathNorm = TRUE,
        normNatConnect = TRUE,
        connectivity = TRUE,
        clustMethod = NULL,
        clustPar = NULL,
        clustPar2 = NULL,
        weightClustCoef = TRUE,
        hubPar = "eigenvector",
        hubQuant = 0.95,
        lnormFit = FALSE,
        weightDeg = FALSE,
        normDeg = TRUE,
        normBetw = TRUE,
        normClose = TRUE,
        normEigen = TRUE,
        verbose = verbosePar
      )
    )
    } else {
      netstats <- "no stats to return"
    }
  return(netstats)
}


# extract node properties -------------------------------------------------
# a function to extract, as a data frame, the node properties from a
# "microNetProps" object


# extracting the node stats 

extract_node_stats <- function(net_stat_list, nodestat){
  node_props <- data.frame(
    pos_degree = apply(net_stat_list$input$assoMat1, 2, function(x) sum(x>0)-1),
    neg_degree = apply(net_stat_list$input$assoMat1, 2, function(x) sum(x<0)),
    clust_memb = net_stat_list[["clustering"]][["clust1"]],
    degree = net_stat_list[["centralities"]][["degree1"]],
    between = net_stat_list[["centralities"]][["between1"]],
    close = net_stat_list[["centralities"]][["close1"]],
    eigenv = net_stat_list[["centralities"]][["eigenv1"]]
  )
  node_props <- node_props %>% 
    mutate(is_hub = rownames(node_props) %in% net_stat_list$hubs$hubs1) %>%
    rownames_to_column("label")
  node_props <- left_join(node_props, nodestat)
  return(node_props)
}


# converting microNet in tidygraph objects --------------------------------

# note that some of the conditions handled by if-else structures are redundant
# and unlikely to happen in practice. The function can probably be simplified

microNet_to_tidygraph <- function(micronet_obj,
                                  net_to_use = 1,
                                  add_names_to_edges = T,
                                  use_asso_matrix = T,
                                  fail_w_err = T
) {
  # assuming you have NetCoMi installed and loaded
  needed_packages <- c("dplyr" , "igraph", "tidygraph")
  if (!all(needed_packages %in% .packages())) {
    cat("installing/loading needed packages\n")
    .to_be_loaded <-
      needed_packages[!(needed_packages %in% .packages())]
    # check if all are in installed packages and install missing
    .inst2 <- .to_be_loaded %in% installed.packages()
    if (any(!.inst2))
      install.packages(.to_be_loaded[!.inst2])
    sapply(.to_be_loaded, require, character.only = TRUE)
  }
  # check the class of the object
  if ((class(micronet_obj) != "microNet")) {
    if (fail_w_err) {
      stop("this function only handles microNet objects\n")
    } else {
      warning("this function only handles microNet objects\n")
      tidygraph_from_micronet <- "not a microNet object"
      return(tidygraph_from_micronet)
    }
  } else {
    # get the adjacency matrix from microNet object
    if (net_to_use == 1) {
      adja_mat <- micronet_obj$adjaMat1
    } else {
      adja_mat <- micronet_obj$adjaMat2
    }
    # check if the adjacency matrix is NULL
    if (is.null(adja_mat) | !is.matrix(adja_mat)) {
      if (fail_w_err) {
        stop("cannot find the adjacency matrix\n")
      } else {
        warning("cannot find the adjacency matrix\n")
        tidygraph_from_micronet <- "not a microNet object"
        return(tidygraph_from_micronet)
      }
    } else {
      # create igraph object from the adjacency matrix
      igraph_from_micronet <- graph_from_adjacency_matrix(
        adja_mat,
        mode = "lower",
        weighted = T,
        diag = F,
        add.colnames = NULL
      )
      # create tidygraph object (from the adjacency matrix) note that this also has class igraph
      tidygraph_from_micronet <- as_tbl_graph(igraph_from_micronet)
      # check if there is at least one edge
      has_edges <- nrow((tidygraph_from_micronet %>% activate(edges) %>% as_tibble())) >= 1
      if(has_edges){
        # optionally add names to edges
        if (add_names_to_edges) {
          tidygraph_from_micronet <- tidygraph_from_micronet %>%
            activate(edges) %>%
            mutate(from_name = .N()$name[from],
                   to_name = .N()$name[to])
        }
        if (use_asso_matrix) {
          # create a network from the association matrix
          if (net_to_use == 1) {
            asso_mat <-
              micronet_obj$assoMat1
          } else {
            asso_mat <- micronet_obj$assoMat2
          }
          if (is.null(asso_mat) | !is.matrix(asso_mat)) {
            if (fail_w_err) {
              stop("cannot find the association matrix\n")
            } else {
              warning("cannot find the association matrix\n")
              return(tidygraph_from_micronet)
            }
          } else {
            asso_graph <- graph_from_adjacency_matrix(
              asso_mat,
              mode = "lower",
              weighted = T,
              diag = F,
              add.colnames = NULL
            )
            tidy_asso_graph <- as_tbl_graph(asso_graph)
            edge_tidy_asso_graph <- tidy_asso_graph %>%
              activate(edges) %>%
              as_tibble() %>%
              dplyr::rename(asso_est = weight) %>%
              mutate(asso_type = if_else(asso_est > 0, "copres", "mut_ex"))
            # join the association measure
            tidygraph_from_micronet <- tidygraph_from_micronet %>%
              activate(edges) %>%
              left_join(., edge_tidy_asso_graph)
            # create an extra column with node labels
            tidygraph_from_micronet <- tidygraph_from_micronet %>%
              activate(nodes) %>%
              mutate(nlabel = name)
          }
        }
        return(tidygraph_from_micronet)
      } else {
        warning("the network has 0 edges\n")
        tidygraph_from_micronet <- "no edges"
        return(tidygraph_from_micronet)
      }
    }
  }
}


# get node stats + calculate edge stats ----------------------------------------------------------
# A function for calculating edge stats (here it does only centrality_edge_betweennes)
# takes a tidygraph object and a further argument which is either a logical or a data frame
# with node stats to be merged with nodes. Returns a modified tidygraph object
# tot_degree is calculated from pos_degree+neg_degree
merge_stats <- function(tg, node_stats = F, ebetw = calc_e_betw){
  if(!is.tbl_graph(tg)) stop("This is not a tidygraph")
  if(ebetw){
    tg <- tg %>% activate(edges) %>% mutate(edge_betw = centrality_edge_betweenness())
  }
  # check if node_stats is a data.frame to decide to merge nodes
  if("data.frame" %in% class(node_stats)){
    # renaming a bariable befor joining
    node_stats <- node_stats %>% rename(name = label)
    tg <- tg %>% activate(nodes) %>%
      left_join(.,node_stats) %>%
      mutate(tot_degree = pos_degree+neg_degree)
  }
  return(tg)
}


# plot networks -----------------------------------------------------------
# c0l0r is either phylum or clust_memb or phylumC
# with phylum, phylum is used for colors but ,ggraphs for different data sets may have different colors
# with clust_memb cluster membership is used
# with phylumC, if available, a fixed scale is used (if available otehrwise phylum is used) 
# a sets alpha for the names of the nodes
# edge_cols changes the color scheme for edges, default is green red but can be cyan magenta
# should be named
# name_lenght sets the length of the node names after which truncation will occur
# this is patchy I should find a way to do it programmatically

# note for self: I could use aes_string
# additionally, I should be more flexible in phylum, may be also allow class

plot_ggraph <- function(tidy_graph, stname = "", method = "", 
                        c0l0r = "phylum", lp = "bottom", 
                        clp = "off", a =0.6, node_label = "genus",
                        name_length = 15,
                        edge_cols = c(copres = "green", mut_ex = "red")){
  # remove nodes with 0 degree (could become an option in the future)
  # make cluster membership a factor
  # safety measure
  if(!(c0l0r %in% c("phylum", "phylumC", "clust_memb"))) stop("c0l0r must be one of 'phylum', 'phylumC', 'clust_memb'")
  g2plot <- tidy_graph %>% 
    activate(nodes) %>%
    dplyr::filter(pos_degree>0 | neg_degree>0) %>%
    mutate(clust_memb = as_factor(clust_memb))
  if(c0l0r == "phylumC"){
    nnames <- colnames(g2plot %>% activate(nodes) %>% as_tibble())
    if(!("phylum_color" %in% nnames)) stop("If you use phylumC for node color you need a phylum_color variable")
  }
  g2plot_title <- paste(stname, method, sep = " ")
  if(c0l0r == "phylum" | c0l0r == "clust_memb"){
    # note that using check_overlap = T may remove the names of some nodes
    ggraph_plot <- ggraph(g2plot, layout = 'fr', weights = weight) + 
      geom_edge_link(mapping = aes(edge_colour = asso_type, edge_width = weight),
                     alpha = I(0.5), show.legend = F) + 
      geom_node_point(mapping = aes(colour = .data[[c0l0r]], size = tot_degree)) +
      geom_node_text(mapping = aes(label = str_trunc(.data[[node_label]], name_length,"center", ellipsis = ".")), check_overlap = F, alpha = a) +
      labs(title = g2plot_title, size = "degree") +
      scale_edge_color_manual(values = edge_cols) +
      scale_edge_width_continuous(range = c(1,4)) +
      coord_cartesian(clip = clp) +  
      theme_graph() +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = lp)
  } else {
    node_color_scale_df <- tidy_graph %>%
      activate(nodes) %>%
      as_tibble() %>% 
      select(phylumC, phylum_color) %>%
      mutate(phylumC = as.character(phylumC)) %>%
      distinct()
    node_color_vector <- node_color_scale_df$phylum_color
    names(node_color_vector) <- node_color_scale_df$phylumC
    g2plot <- g2plot %>% activate(nodes) %>% mutate(phylumC = fct_drop(phylumC))
    ggraph_plot <- ggraph(g2plot, layout = 'fr', weights = weight) + 
      geom_edge_link(mapping = aes(edge_colour = asso_type, edge_width = weight),
                     alpha = I(0.5), show.legend = F) + 
      geom_node_point(mapping = aes(colour = .data[[c0l0r]], size = tot_degree)) +
      geom_node_text(mapping = aes(label = str_trunc(.data[[node_label]], name_length,
                                                     "center", ellipsis = ".")), 
                     check_overlap = F) +
      labs(title = g2plot_title, size = "degree") +
      scale_edge_color_manual(values = edge_cols) +
      scale_edge_width_continuous(range = c(1,4)) +
      scale_color_manual(values = node_color_vector) +
      theme_graph(base_family="sans") +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = lp) 
  }
  return(ggraph_plot)
}


# oddsratio_assotype ------------------------------------------------------

# A function to calculate odds ratio for copresence and mutual exclusion
# associations as a function of taxonomic relationships of the two nodes
# connected by an association edge. 
# Known issues: epi2by2() function has changed with epiR 2.0.26; 
# The odds ratio function

odds_ratio <- function(inputdf, taxo_level = "family"){
  epiR_version <- packageVersion("epiR")
  if(epiR_version < "2.0.26"){
    stop("\nYou must update epiR to version 2.0.26 or higher!")
  }
  mini_df <- switch(taxo_level,
                    family = select(inputdf, asso_type, sf),
                    order = select(inputdf, asso_type, so),
                    class = select(inputdf, asso_type, sc)
  )
  colnames(mini_df)[2] <- "same_taxon"
  # copresence
  epi_list_cop <- epiR::epi.2by2(xtabs(~same_taxon + asso_type, data = mini_df))
  # extract in a data frame the Wald incidence risk ratio
  # extracts the appropriate chisq values
  chicop <- if(epi_list_cop$massoc.detail$chi2.correction){
    epi_list_cop$massoc.detail$chi2.strata.yates
  } else {
    epi_list_cop$massoc.detail$chi2.strata.uncor
  }
  OR <- cbind(epi_list_cop$massoc.detail$OR.strata.wald, chicop)
  names(OR) <- str_c("OR", names(OR), sep = "_")
  RR <- epi_list_cop$massoc.detail$RR.strata.wald
  names(RR) <- str_c("RR", names(RR), sep = "_")
  assort_results <- cbind(
    OR,
    RR
  )
  return(assort_results)
}

# df_return
# checks a list containing data frames or try-error objects, drops the try-error objects
# and binds the data frames by row
df_return <- function(input_list){
  classes <- map_chr(input_list, class)
  input_list<-input_list[which(classes == "data.frame")]
  df <- bind_rows(input_list, .id = "method")
  return(df)
}


# a report for the functions, a tribble
print_function_report <- function(){
  ext_functions <- tribble(
    ~funct_name,              ~description,
    "load metadata",          "loads and checks the metadata",
    "report_step_0",          "create a report for step 0 (original phyloseqs)",
    "prune_samples_by_size",  "removes samples with less than min_seqs, optionally prints ECDF",
    "rem_low_sample_obj",     "removes phyloseq objects with less than min_samples",
    "report_step_n",          "create a report for step n (from sample pruning onward)",
    "gset_rank_names",        "get and set rank names",
    "remove_Chl_Mit",         "removes chloroplasts and mitochondria",
    "tax_glom_name_change",   "performs taxa agglomeratio and changes name",
    "filter_by_prev_ab",      "performs prevalence and abundance filtering, with optional graphs and tables",
    "infer_MAN",              "infers microbial association networks using netConstruct",
    "calculate_net_stats",    "calculate network and node statistics",
    "extract_node_stats",     "extracts node properties from a microNetProps object",
    "microNet_to_tidygraph",  "convert a microNet object to a tidygraph object",
    "merge_stats",            "merge node stats in the tidygraph object, optionally calculates edge betweenness",
    "plot_ggraph",            "use ggraph to plot the tidygraph, and returns the object",
    "odds_ratio",             "calculates odds ratios for taxonomic assortativity",
    "df_return",              "takes a list of data frames and try-error objects, removes the latter and binds the data frames"
  )
  return(ext_functions)
}