# filter_my_taxa_1_2_1 17/6/25
# filter_my_taxa ----------------------------------------------------------

# a function for filtering taxa based on prevalence and abundance thresholds or
# permutation filtering (see https://doi.org/10.1093/biostatistics/kxy020)
# takes a phyloseq object 
# and other options and returns a list
# myphyseq is the phyloseq object to process
# name is the physeq object name (will be taken from the list during processing)
# prevfilter is either "prevab" (for a prevalence and abundance filter)
# or "PERFect" for a permutation filter
# PERFect_options are the options for the permutation filter, can be passed as list
# prevthreshold the prevalence threshold, as fraction of samples
# passboth a flag, T if both the prevalence and abundance filter should be passed
# abthreshold the abundance threshold
# filenm a prefix which will be used in filenames
# colvar_pattern is the level of taxonomic aggregation in the prevalence
# and abundance plot, can be "[P|p]hylum" (the default) or "[C|c]lass",
# gres the resolution of graphs
# gtype graph type ("jpeg", "tiff", depending on your system)
# saveplot a flag: should the prevalence and abundance plot be saved
# printplot a flag: should the prevalence and abundance plot be printed
# outfolder the output folder
# savepat a flag: should the prevalence and abundance table be saved

# the function works with defaults and returns a list

# install/load packages 
.scran_packages <- c("tidyverse", "tictoc")
.sbioc_packages <- c("BiocManager", "phyloseq")
.attached_pckgs <- .packages()

# install and load Bioconductor packages if necessary
if(!all(.sbioc_packages %in% .attached_pckgs)){
  .sinst <- .sbioc_packages %in% installed.packages()
  if(any(!.sinst)) {
    if(!.sinst[1]) {
      install.packages("BiocManager")
      .sinst <- .sbioc_packages %in% installed.packages()
    }
    if(any(!.sinst[2:length(.sinst)])) {
      BiocManager::install(.sbioc_packages[!.sinst], ask = F)
    }
  }
  sapply(.sbioc_packages, require, character.only = TRUE)
}

# install and load CRAN packages if necessary 
if(!all(.scran_packages %in% .attached_pckgs)){
  .sinst <- .scran_packages %in% installed.packages()
  if(any(!.sinst)) {
    install.packages(.scran_packages[!.sinst])
  }
  sapply(.scran_packages, require, character.only = TRUE)
}

# PERFect is not available any more, loading from source

filter_my_taxa <- function(myphyseq, 
                           name = "phyloseqobj",
                           prevfilter = "prevab",
                           PERFect_options = list(
                             method = "sim", # "sim" or "perm", perm takes longer
                             perf_alpha = 0.1, # this is the default
                             perf_k = 10000, # this is the number of permutation, default value
                             rollm = F, # this is the rolling mean argument
                             perf_alg = "fast", # algorithm for permutation filtering fast or full, full is slower
                             perf_hist = FALSE, # the hist option in PERFect_perm
                             vrb_out = TRUE # if extra output needs to be produced
                           ),
                           prevthreshold = 0.02,
                           passboth = T,
                           abthreshold = 0.001,
                           filenm = "prevfilter",
                           colvar_pattern = "[P|p]hylum",
                           gres = 150,
                           gtype = "tiff",
                           saveplot = T,
                           printplot = T,
                           outfolder = ".",
                           savepat = T){
  # creates list for the results, the object which will be returned by the function
  tax_filter_results <- list(phyobj = NULL,
                             phyobjname=NULL,
                             prevabplot=NULL,
                             prevabtable=NULL)
  # check if any taxon has 0 counts and removes it
  if(any(taxa_sums(myphyseq)==0)){
    ntaxa_pre <- ntaxa(myphyseq)
    myphyseq <- filter_taxa(myphyseq, function(x) sum(x)>0, prune = T)
    ntaxapost <- ntaxa(myphyseq)
    taxaremoved <- ntaxa_pre-ntaxapost
    if(verbose_output) red(cat("\n",taxaremoved," taxa were removed because of 0 sums","\n"))
  }
  # check if any samples has 0 counts and removes it
  if(any(sample_sums(myphyseq)==0)){
    nsamples_pre <- nsamples(myphyseq)
    myphyseq <- subset_samples(myphyseq, sample_sums(myphyseq)>0)
    nsamplespost <- nsamples(myphyseq)
    samplesremoved <- nsamplespre-nsamplespost
    if(verbose_output) red(cat("\n",samplesremoved," samples were removed because of 0 sums","\n"))
  }
  # get the option on the OTU table
  taxa_are_rows_flag <- myphyseq@otu_table@taxa_are_rows
  # obtain a prevalence and abundance plot
  OTUmatrixf <- as(otu_table(myphyseq), "matrix")
  if (taxa_are_rows_flag) {
    OTUmatrixf <- t(OTUmatrixf)
  }
  OTUmatrixf_relab <- OTUmatrixf / rowSums(OTUmatrixf)
  # calculate the prevalence on all columns
  prevdf <- apply(
    X = OTUmatrixf,
    MARGIN = 2,
    FUN = function(x) {
      sum(x > 0)
    }
  )
  # calculate minimum relative abundance
  min_rel_ab <- apply(
    X = OTUmatrixf_relab,
    MARGIN = 2,
    FUN = function(x) {
      min(x)
    }
  )
  # calculate max. rel abundance
  max_rel_ab <- apply(
    X = OTUmatrixf_relab,
    MARGIN = 2,
    FUN = function(x) {
      max(x)
    }
  )
  prevdf <- data.frame(
    Prevalence = prevdf,
    TotalAbundance = colSums(OTUmatrixf),
    min_rel_ab = min_rel_ab,
    max_rel_ab = max_rel_ab
  )
  # merge taxonomic information
  taxa_metadata <-
    as.data.frame(as(tax_table(myphyseq),"matrix")) %>% rownames_to_column(var = "label")
  tranks <- rank_names(myphyseq)
  prevdf <- prevdf %>%
    rownames_to_column(var = "label") %>%
    left_join(., select(taxa_metadata, label, any_of(tranks)))
  
  # apply the filter 
  # "prevab" is a prevalence and abundance filter
  # "PERFect" is a permutation filter, see https://doi.org/10.1093/biostatistics/kxy020.
  if (prevfilter == "prevab") {
    pass_prev_filter <-
      dplyr::filter(select(prevdf, label, Prevalence),
                    Prevalence > floor(nsamples(myphyseq) *
                                         prevthreshold)) %>%
      pull(label)
    if (passboth) {
      OTUtokeep <- intersect(names(which(max_rel_ab >= abthreshold)),
                             pass_prev_filter)
    } else {
      OTUtokeep <- union(names(which(max_rel_ab >= abthreshold)),
                         pass_prev_filter)
    }
  } else {
    # use PERFect
    # create the count data frame
    if(taxa_are_rows_flag){
      counts <- as.data.frame(t(otu_table(myphyseq)))
    } else {
      counts <- as.data.frame(otu_table(myphyseq))
    }
    if(PERFect_options$method == "sim"){
      tic("start PERFect_sim")
      res_sim <- PERFect_sim(X = counts, alpha = PERFect_options$perf_alpha)
      toc()
      if(PERFect_options$vrb_out){
        p <- pvals_Plots(PERFect = res_sim, X = counts, quantiles = c(0.25, 0.5, 0.8, 0.9), 
                         alpha=PERFect_options$perf_alpha)
        print(p$plot + ggtitle(str_c(filenm," Simultanenous Filtering")))
      }
      OTUtokeep <- colnames(res_sim$filtX)
    } else {
      tic("start PERFect_perm")
      res_sim <- PERFect_perm(X = counts, alpha = PERFect_options$perf_alpha, 
                              algorithm = PERFect_options$perf_alg, 
                              hist = PERFect_options$perf_hist,
                              rollmean = PERFect_options$rollm)
      toc()
      if(PERFect_options$vrb_out){
        p <- pvals_Plots(PERFect = res_sim, X = counts, quantiles = c(0.25, 0.5, 0.8, 0.9), 
                         alpha= PERFect_options$perf_alpha)
        print(p$plot + ggtitle(str_c(filenm,
                               "Permutation Filtering", 
                               PERFect_options$perf_alg, "algorithm", sep = " ")))
      }
      OTUtokeep <- colnames(res_sim$filtX)
    }
    
    
  }
  prevdf <- prevdf  %>%
    mutate(
      relAbundance = TotalAbundance / sum(TotalAbundance),
      pass_filters = ifelse(label %in% OTUtokeep, "T", "F")
    )
  # apply the filter and put the result in the return list
  tax_filter_results$phyobj <- prune_taxa(OTUtokeep, myphyseq)
  tax_filter_results$phyobjname <- name
  
  # prevalence vs abundance plot
  OTUmatrixf <-
    OTUmatrixf[, which(colnames(OTUmatrixf) %in% OTUtokeep)]
  # fraction of remaining sequences
  
  f_seq_ret <-
    round(sum(sample_sums(tax_filter_results$phyobj)) / sum(sample_sums(myphyseq)), 4)
  
  # make a plot
  title_text <- paste("Prevalence vs. abundance, by Phylum,", name, sep = " ")
  # original number of taxa
  ntaxa_prefilt <- ntaxa(myphyseq)
  subtitle_text <-
    paste(
      "using the filters you retain ",
      length(OTUtokeep),
      " taxa (triangles) out of ",
      ntaxa_prefilt,
      " (",
      f_seq_ret * 100,
      "% of init. seqs.)",
      sep = ""
    )
  # a prevalence and abundance plot change rank for color and facet as appropriate
  colorvar <- rank_names(myphyseq)[which(str_detect(rank_names(myphyseq), colvar_pattern))]
  prev_ab_plot <-
    ggplot(prevdf,
           aes(
             x = TotalAbundance,
             y = Prevalence / nrow(OTUmatrixf),
             shape = as.factor(pass_filters),
             color = .data[[colorvar]]
           )) +
    geom_point(size = 2, alpha = 0.7) +
    facet_wrap(~ .data[[colorvar]]) +
    geom_hline(
      yintercept = ifelse(prevfilter == "prevab", prevthreshold, 0),
      alpha = 0.5,
      linetype = 2
    ) +
    labs(
      x = "total abundance",
      y = "Prevalence [Frac. Samples]",
      shape = 'pass ab. threshold',
      title = title_text,
      subtitle = subtitle_text
    ) +
    scale_x_log10() +
    scale_y_continuous(minor_breaks = seq(0, 1, 0.05)) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 90)
    )
  if(printplot) print(prev_ab_plot)
  # add the plot to the list
  tax_filter_results$prevabplot <- prev_ab_plot
  
  if(saveplot) ggsave(prev_ab_plot,
                      file = paste(file.path(outfolder,filenm), 
                                   "_prevabg.", gtype, sep=""),
                      device = gtype,
                      dpi = gres, 
                      width = 7,
                      height = 7)
  # finish up with the table
  prevdf <- prevdf %>%
    mutate(relprev = Prevalence / nrow(OTUmatrixf)) %>%
    arrange(-relprev,-relAbundance)
  tax_filter_results$prevabtable <- prevdf
  
  if(savepat) {
    write_tsv(prevdf, 
              file = paste(file.path(outfolder,filenm), "_prevabt.txt",sep=""))
  }
  return(tax_filter_results)
}

# Credits and copyright notice

# The code for loading packages is derived from [this post](http://tinyurl.com/jjwyzph).  

# The code for PERFect functions is taken verbatim from 
# Smirnova E, Cao Q (2022). PERFect: Permutation filtration for microbiome data. R package version 1.12.0, https://github.com/cxquy91/PERFect.
# this is because most recent versions of Bioconductor do not support PERFect

# The code for the creation of the prevalence and abundance data frame is taken from here: 
# https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html  
# and significantly modified