# METAOlive_MAN v2.3
# 5/11/24

# a workflow for estimating microbial association networks from phyloseqs 
# of fungal and bacterial communities from the olive FMBN dataset

# recycle MAN_dairy_whole studies with modifications
# process bacteria and fungi separately, in two iterations
# do post-processing
# structure of networks and their properties
# identification of keystone taxa
# visualization and post-processing using scripts developed in other projects

# in a separate script
# differential network analysis when feasible: 
# starter/no starter
# spoiled/no spoiled for Penland studies

# a script designed to infer microbial association networks on whole studies
# on table olive products extracted from FoodMicrobionet
# Part of the 
# "METAolive" project, WP3 

# Install/load packages ---------------------------------------------------
# install packages, set general parameters and options
# THESE MUST BE RUN EVERY TIME YOU REOPEN THE PROJECT
# including the assignment of nc

.bioc_packages <- c("BiocManager", "phyloseq", "PERFect")
.cran_packages <- c("tidyverse", "tictoc", "beepr", "logr", "VennDiagram", 
                    "lobstr", "parallel", "tidygraph", "ggraph", "cowplot",
                    "RColorBrewer")

.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  if(!.inst[1]) {
    install.packages("BiocManager")
    .inst <- .bioc_packages %in% installed.packages()
  }
  if(any(!.inst[2:length(.inst)])) {
    BiocManager::install(.bioc_packages[!.inst], ask = F)
  }
}

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}

# lobstr is installed for obj_size
# cowplot for arranging plots created by ggplot in grids
# RColorBrewer to create the palette for phyla


# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

# NetCoMi needs to be installed from GitHub; the installation is somewhat
# complicated and some features would not work in the arm64 version of R
# see https://github.com/stefpeschel/NetCoMi#installation for details
require(NetCoMi)

# the following command detects the number of cores on UNIX/MacOS
ncores <- parallel::detectCores(logical = F) # to detect physical cores in MacOS
# other setup operations
opar <- par(no.readonly=TRUE) 
par(ask=F) 
# set.seed(1234) 


# further parameters ------------------------------------------------------
# play audio notifications
play_audio <- T
# time important steps
keep_time <- T
# verbose output: will print additional objects and messages
verbose_output <- T

# debug_mode: will save an error log with specific error messages
debug_mode <- T
if(debug_mode){
  if (!dir.exists("debugging_logs")) dir.create("debugging_logs")
}


# Processing and folders --------------------------------------------------
# the studies are ran in groups 0 is the trial run

iteration <- 1

# the name of the folder containing the phyloseq objects
data_folder <- "physeqs_F"
# the name of the folder containing additional functions + optional data
source_folder <- "source"
# the name of the output folder
output_folder <- file.path("output", "MAN_F")

# check if the output folder is there, if not, create it
if(sum(str_detect(list.dirs(), file.path(".","output")))==0) dir.create("output")
if(sum(str_detect(list.dirs(), file.path(".","output", "MAN_F")))==0) dir.create(file.path("output","MAN_F"))

# fn_pref is the prefix for all output filenames
fn_pref <- "whstudies_genus"
fn_pref_n <- str_c("whstudies_genus",iteration)


process_batch <- F
# if true, all phyloseqs in the data folder will be processed in batch (may take time)

# some graph options, shared everywhere
dpi_option <- 300
g_type <- "tiff" 


# Filtering options -------------------------------------------------------

min_samples <- 15
# the minimum number of samples in the data set: studies with < min_samples after filtering will be excluded from the analysis
min_seqs <- 1000
# the minimum number of sequences per sample

# removes uncharacterized taxa (i.e. taxa missing the identification at the phylum level or above)
# only needed if no tax glom is performed, set to F for fungi
rm_unchar <- F
# remove Eukaryota
rm_euk <- F
# remove chloroplasts and mitochondria
rm_chlmit <- F
# taxonomic agglomeration options: "none" means no aggregation; if the input data are ASVs or OTUs they will be used as such
# if the input file is from FMBN data in the original datasets are aggregated at different taxonomic levels.
# when comparing datasets obtained using different gene targets or pipelines should be set to "genus"
# "none", "genus", or "family"); do not use species, will most likely cause errors
taxglom <- "genus" # tax glom will be performed at the genus level
# remove taxa which are uncharacterized at the family, order and class level (should not be very many in most datasets)
# only needed if tax glom is not performed
above_genus_flag <- T

# filtering taxa
# flag to decide if OTU should be filtered, i.e. if the filter is to be applied at all
filterOTUs <- TRUE

# prev_filter is either "prevab" (for a prevalence and abundance filter)
# or "PERFect" for a permutation filter, otherwise will throw an error (sooner or later)
prev_filter <- "PERFect"
# options for the prevalence and abundance filter
prev_thr <- 0.02
ab_thr <- 0.001
pass_both <- TRUE 
# if true an OTU must pass both abundance and prevalence filter, otherwise
# passing either is enough

# options for the permutation filter
perf_opt <- list(
  method = "perm", # "sim" or "perm", perm takes longer
  perf_alpha = 0.1, # this is the default
  perf_k = 10000, # this is the number of permutations, default value
  rollm = F, # this is the rolling mean option of PERFect_perm
  perf_alg = "fast", # algorithm for permutation filtering fast or full, full is slower
  perf_hist = FALSE, # the hist option in PERFect_perm
  vrb_out = TRUE # printing extra graphs
)

# options for saving output of prevalence and abundance filter
save_prev_ab_plot <- T
# flag for saving the prev ab plot
print_prev_ab_plot <- F
# flag for printing the plot
save_prev_table <-T 
# flag for saving the prevalence and abundance table; should be saved because it can be used later
# to add info to the node stats
save_prev_ab_list <-T

# a vector with the methods used for inference, use the spelling used in `netConstruct()`
# seeNetCoMi documentation for a lsit of available methods and options;
# the first method is the one which will be used to do comparisons within sample
# if you want to run just one list only that
# I have removed "spring" which is too intensive in terms of calculations
inf_methods <- c("spieceasi", "sparcc", "ccrepe")

# the list of parameters to be passed in `netConstruct()` for each of the methods;
# there must be one entry for each of he parameters in the list
# some parameters are not really needed and will be ignored
inf_methods_param <- list(
  spieaceasi = list(
    sparMethod = "t-test",
    alpha = 0.001,
    measureParList = list(
      method = 'mb',
      lambda.min.ratio = 1e-2,
      nlambda = 25,
      pulsar.params = list(rep.num =
                             50)
    ),
    normmethodPar = "none",
    zeromethodPar = "none",
    dissFuncPar = "signed",
    verbosePar = 1
  ),
  # spring = list(
  #  sparMethod = "t-test",
  #  alpha = 0.001,
  #  measureParList = list(nlambda = 50, rep.num = 50),
  #  normmethodPar = "none",
  #  zeromethodPar = "none",
  #  dissFuncPar = "signed",
  #  verbosePar = 1
  #),
  sparcc = list(
    sparMethod = "t-test",
    alpha = 0.001,
    measureParList = list(iter = 100, inner_it = 20, th = 0.05),
    normmethodPar = "none",
    zeromethodPar = "none",
    dissFuncPar = "signed",
    verbosePar = 1
  ),
  ccrepe = list(
    sparMethod = "t-test",
    alpha = 0.001,
    measureParList = NULL,
    normmethodPar = "fractions",
    zeromethodPar = "none",
    dissFuncPar = "signed",
    verbosePar = 1
  )
)

# options for network analysis
# do Venn plot
do_Venn <- T
# calculates edge betweenness
calc_e_betw <-T
# merge node stats
merge_node_stats <- T

# functions ---------------------------------------------------------------

# load the filter_my_taxa_function
# allows to choose between a prevalence and abundance filter and a 
# permutation filter (package PERFect)
source(file.path("source","filter_my_taxa_function_1_1.R"))

# load functions for network inference and plotting
source(file.path("source","net_anal_funct_v1_1.R"))
# I have modified one of the funtions to allow species names in nodes, although
# I doubt is a good idea; same thing with node colors, now it is phylum, could be class

# phyloseq objects --------------------------------------------------------

if(keep_time) tic("complete batch operation")

# I need to remove all the studies which have a suffix like sample, environment
all_physeqs <- list.files(path = file.path(".", data_folder))
physeqs_to_use <- !str_detect(all_physeqs, "environment|sample")
# the phyloseqs I need to use
wh_studies_physeqs <- all_physeqs[physeqs_to_use]
study_ids <- str_split_i(wh_studies_physeqs, "_", 1)
rm(all_physeqs, physeqs_to_use)
if(!exists("study_groups")){
  study_groups <- data.frame(studyId = study_ids,
                             physeqs = wh_studies_physeqs) |> 
    separate(studyId, into = c("throwaway","stid"), sep = "T", remove = F) |>
    mutate(stid = as.numeric(stid)) |>
    arrange(stid) |>
    select(-throwaway)
  study_groups$group <- c(rep(1,4), rep(2,5))
  write_tsv(study_groups, file =str_c(fn_pref,"_studygroups.txt"))
} else {
  study_groups <- read.delim(file =file.path(getwd(),str_c(fn_pref,"_studygroups.txt")))
}

# this needs to be done only once: need a study_info file and a phylum
# palette file
if(!"whole_studies_info.txt" %in% list.files()){
  # create the study info and the palette from scratch and save them
  # I am using FMBN_plus, version 4.2.1, in a folder named FMBN
  FMBN_path <- file.path("FMBN", "FMBN_plus.rds") 
  FMBN <- readRDS(FMBN_path)
  study_info <- FMBN$studies
  study_info <- study_info |> 
    dplyr::filter(studyId %in% study_ids)
  # calculate cumulative sum of samples
  study_info <- study_info |> mutate(cuM_sum_samples=cumsum(samples))
  # there is 880 samples in the studies
  # need to divide in 5 groups approx the same size
  write_tsv(study_info, file = "whole_studies_info.txt")
  # get the 10 most abundant phyla and create and save a palette
  # get the samples 
  samples <- FMBN$samples_F |>
    dplyr::filter(studyId %in% study_ids)
  # get the edges
  edges <- FMBN$edges_F |>
    dplyr::filter(sampleId %in% samples$sampleId)
  # get taxa (and remove Chloroplast and Mitochondria and poorly identified stuff)
  taxa <- FMBN$taxa |>
    dplyr::filter(!is.na(phylum)) |>
    dplyr::filter(family != "Mitochondria") |>
    dplyr::filter(class != "Chloroplast" | order != "Chloroplast") |>
    dplyr::filter(domain == "Fungi")
  # merge phylum info into edges
  edges <- edges |> 
    left_join(dplyr::select(taxa, taxonId, phylum)) |>
    dplyr::filter(!is.na(phylum))
  # summarize by sampleId and phylum
  edges_sum <- edges |>
    summarise(phy_ab = sum(weight), .by = c(sampleId,phylum))
  # now the phyla with the highest average abundance
  phylum_mean_ab <- edges_sum |>
    summarise(ave_ab = mean(phy_ab), .by = phylum) |>
    arrange(desc(ave_ab)) 
  n_phyla <- n_distinct(phylum_mean_ab$phylum)
  n_phyla <- if_else(n_phyla>10, 10, n_phyla)
  phylum_palette <- data.frame(
    phylum = factor(c(phylum_mean_ab$phylum[1:n_phyla], "Other"), level = c(phylum_mean_ab$phylum[1:n_phyla], "Other")),
    phylum_color = brewer.pal(n_phyla+1, "Set3")
  )
  write_tsv(phylum_palette, file =str_c(fn_pref,"_phylumpalette.txt"))
  # clean up
  rm(FMBN, FMBN_path, samples, edges, taxa, edges_sum, phylum_mean_ab, n_phyla)
} else {
  # get the info from files, file.path needed if working on files in Dropbox
  study_info <- read.delim(file = file.path(".","whole_studies_info.txt"))
  phylum_palette <- read.delim(file= file.path(".",str_c(fn_pref,"_phylumpalette.txt")))
}


if(iteration == 0){
  # 0 is the trial run
  studies_to_filter <- str_c("ST182", "_DNA.rds")
} else {
  studies_to_filter <- study_groups |>
    dplyr::filter(group == iteration) |>
    dplyr::select(physeqs) |>
    pull()
}


study_info <- study_info %>%
  dplyr::filter(studyId %in% str_split_i(studies_to_filter, "_", 1))
n_physeqs <- length(studies_to_filter)

physeq_files <- map_chr(studies_to_filter, 
                        ~file.path(getwd(),"physeqs_F",.x))
names(physeq_files) <- studies_to_filter

physeq_list <- vector(mode = "list", length = n_physeqs)


physeq_list <- map(physeq_files, readRDS)

# need to remove suffix from names
names(physeq_list) <- str_split_i(names(physeq_list), "\\.",1)


# just in case (you should know what you are doing)

is_phyloseq <- sapply(physeq_list, function(x) class(x)=="phyloseq")
# clean up the list if necessary and give names
physeq_list <- physeq_list[is_phyloseq]

# should check if this is too large
# size_limit 
size_limit <- 500e6L
physeq_size <- obj_size(physeq_list)
size_warning <- ifelse(as.numeric(physeq_size)>size_limit,
                       "WARNING: too much data to process",
                       paste("object size", as.character(physeq_size), " B", sep = " ")) 



# prefilter and glom (genus level) ----------------------------------------

# remove mitochondria, chloroplasts
# quickly checking rank names
rank_names(physeq_list[[1]])

if(rm_chlmit){
  physeq_list_filt_1 <- physeq_list |>
    map(\(x) subset_taxa(x, !(family == "Mitochondria" | 
                                class == "Chloroplast" | order == "Chloroplast")))
} else {
  physeq_list_filt_1 <- physeq_list 
}

# remove Eukaryota
if(rm_euk){
  physeq_list_filt_2 <- physeq_list_filt_1 |> 
    map(\(x) subset_taxa(x, !(domain == "Eukaryota")))
} else {
  physeq_list_filt_2 <- physeq_list_filt_1  
}


# perform tax glom -------------------------------------------------------

# perfom tax_glom at the taxglom level (this will also remove anything
# whose taxonomic assignment is above taxonomic level stored in taxglom )

# tax glom at the taxglom level
if(taxglom != "none"){
  if(verbose_output) tic("starting tax glom")
  physeq_list_filt_3 <- physeq_list_filt_2 |> 
    map(\(x) tax_glom(x, taxglom))
  toc()
} else {
  physeq_list_filt_3 <- physeq_list_filt_2
}

# rename taxa if needed
if(taxglom != "none"){
  for(i in seq_along(physeq_list_filt_3)){
    my_temp_physeq <- physeq_list_filt_3[[i]]
    my_taxa_rank_glom <- which(rank_names(my_temp_physeq) == taxglom)
    # check for occurrences of duplicated taxanames (a rare occurrence)
    taxanames <- as.data.frame(tax_table(my_temp_physeq)) |>
      pull(taxglom)
    if(!any(duplicated(taxanames))){
      taxa_names(my_temp_physeq) <- taxanames
    } else {
      my_dupli_taxa <- taxanames[duplicated(taxanames)]
      my_temp_taxa_table <- as.data.frame(tax_table(my_temp_physeq))
      for(j in seq_along(my_dupli_taxa)){
        my_taxon <- my_dupli_taxa[j]
        my_taxon_index <- which(my_temp_taxa_table[,my_taxa_rank_glom] == my_taxon)
        # recover the taxonomy
        tax_replacement <- as.data.frame(tax_table(my_temp_physeq))[my_taxon_index[1],]
        my_temp_physeq <- merge_taxa(my_temp_physeq, eqtaxa = my_taxon_index, 
                                     archetype = 1L)
        my_temp_taxa_table <- as.data.frame(tax_table(my_temp_physeq))
        my_temp_taxa_table[my_taxon_index[1],] <- tax_replacement
        tax_table(my_temp_physeq) <- as.matrix(my_temp_taxa_table)
      }
      taxanames <- as.data.frame(tax_table(my_temp_physeq)) |>
        pull(taxglom)
      taxa_names(my_temp_physeq) <- taxanames
    }
    physeq_list_filt_3[[i]] <- my_temp_physeq
    rm(taxanames, my_temp_physeq)
  }
}

# Calculate diversity prior to filtering for prevalence and abundance --------
# this is done after preliminary filtering and tax_glom if any

if(keep_time) tic("calculate diversity pre-filter")
# may generate warnngs
div_est_prefilter <- map(physeq_list_filt_3, 
                         \(x) phyloseq::estimate_richness(x, 
                                                          split = F, 
                                                          measure=c("Observed","Chao1","Shannon"))) |>
  list_rbind(names_to = "dataset") |>
  remove_rownames() |>
  mutate(Pielou_J = Shannon/log(Observed))
# the last step adds Pielou_J evenness

# calculate and add ave Bray-Curtis dissimilarity
meanbcdist <- map(physeq_list_filt_3, 
                  \(x) phyloseq::distance(x, method="bray"))
div_est_prefilter$ave_BC <- unlist(map(meanbcdist, mean))


# save for further use
save(div_est_prefilter, 
     file = file.path(getwd(),
                      output_folder,
                      str_c(fn_pref_n, "_divprefilter.Rdata")))


# add this to study info
study_metadata <- data.frame(dataset = studies_to_filter) |>
  mutate(studyId = str_split_i(dataset, "_", 1)) |>
  mutate(dataset = str_split_i(dataset, "\\.", 1))

study_metadata <- left_join(study_metadata, div_est_prefilter) 

if(keep_time) toc()


# create a list with filtered phyloseqs -----------------------------------


if(keep_time) tic(msg = "perform taxa filtering")

# let's check
if(!(prev_filter %in% c("prevab", "PERFect"))) stop("wrong specification of OTU filter")


# neet to use a loop for handling names
physeq_list_filt_4 <- vector(mode = "list", length = n_physeqs)
names(physeq_list_filt_4) <- names(physeq_list_filt_3)
# I had to relax filtering criteria: use ov prevab,
# large prev_thr; this is likely due to low diversity and size
for(i in seq_along(physeq_list_filt_3)){
  fnm <- names(physeq_list_filt_4)[i]
  if(prev_filter == "prevab"){
    physeq_list_filt_4[[i]] <- filter_my_taxa(
      physeq_list_filt_3[[i]], 
      name = names(physeq_list_filt_4)[i],
      prevfilter = "prevab",
      prevthreshold = prev_thr,
      abthreshold = ab_thr,
      filenm = str_c(fnm,
                     "prevab",
                     sep="_"),
      outfolder = output_folder,
      savepat = T
    )
  } else {
    physeq_list_filt_4[[i]] <- filter_my_taxa(
      physeq_list_filt_3[[i]], 
      name = names(physeq_list_filt_4)[i],
      prevfilter = "PERFect",
      PERFect_options = perf_opt,
      filenm = str_c(fnm,"PERFect",sep="_"),
      outfolder = output_folder,
      saveplot = T,
      printplot = T,
      savepat = T
    )
  }
}

sapply(physeq_list_filt_4, class)
rm(fnm)

# need to remove further unidentified taxa
if(taxglom == "none" & above_genus_flag){
  for(i in seq_along(physeq_list_filt_4)){
    phy <- physeq_list_filt_4[[i]][[1]]
    phy <- subset_taxa(phy, !is.na(genus))
    physeq_list_filt_4[[i]][[1]] <- phy
  }
  rm(phy)
}

if(keep_time) toc()


# results of filtering ----------------------------------------------------

original_data <- data.frame(n_taxa = map_dbl(physeq_list_filt_1, ntaxa),
                            n_seqs = map_dbl(physeq_list_filt_1, \(x) sum(sample_sums(x)))) %>%
  mutate(step = "orig_data") %>%
  rownames_to_column("study") 
first_tax_filter <- data.frame(n_taxa = map_dbl(physeq_list_filt_2, ntaxa),
                               n_seqs = map_dbl(physeq_list_filt_2, \(x) sum(sample_sums(x)))) %>%
  mutate(step = "first_tax_filter") %>%
  rownames_to_column("study") 
after_tax_glom <- data.frame(n_taxa = map_dbl(physeq_list_filt_3, ntaxa),
                             n_seqs = map_dbl(physeq_list_filt_3, \(x) sum(sample_sums(x)))) %>%
  mutate(step = "after_tax_glom") %>%
  rownames_to_column("study")
n_taxa_summary <- bind_rows(original_data, first_tax_filter, after_tax_glom)

phy4taxa <- vector(mode="list", length = n_physeqs)

for(i in seq_along(physeq_list_filt_4)){
  phy4taxa[[i]] <- data.frame(n_taxa = ntaxa(physeq_list_filt_4[[i]][[1]]),
                              n_seqs = sum(sample_sums(physeq_list_filt_4[[i]][[1]])))
  names(phy4taxa)[i] <- names(physeq_list_filt_4)[i]
}
phy4 <- list_rbind(phy4taxa, names_to = "study") |>
  dplyr::mutate(step = "tax_filt") |>
  dplyr::select(study, n_taxa, n_seqs, step)
n_taxa_summary <- bind_rows(n_taxa_summary, phy4)
# cleanup
rm(i, phy4)

# need to standardize
n_taxa_summary_2 <- n_taxa_summary %>%
  group_by(study) %>%
  mutate(prop_taxa = n_taxa/max(n_taxa),
         prop_seqs = n_seqs/(max(n_seqs))) %>%
  ungroup()


ggplot(n_taxa_summary_2, mapping=aes(x=study, y =prop_taxa, color = step)) +
  geom_point() +
  ggtitle("Proportion of taxa remaining after each filter step") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(
  filename = file.path(getwd(), output_folder, str_c(fn_pref_n, "_taxaloss_perm.", g_type)), 
  dpi = dpi_option
  )

ggplot(n_taxa_summary_2, mapping=aes(x=study, y =prop_seqs, color = step)) +
  geom_point() +
  ggtitle("Proportion of sequences remaining after each filter step") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(
  filename = file.path(getwd(), output_folder, str_c(fn_pref_n, "seqloss_perm.", g_type)), 
  dpi = 300
  )


# let's save the list 
save(physeq_list_filt_4, 
     file = file.path(getwd(), output_folder, str_c(fn_pref_n, "phylist_perm_filtered.Rdata")))

# save the seq loss summary
write_tsv(
  n_taxa_summary_2,
  file = file.path(getwd(), output_folder, str_c(fn_pref_n, "filt_summary.txt"))
)

# make a list with only the physeq using a combiantion of map and pluck
physeq_list_filt_5 <- map(physeq_list_filt_4,
                          \(x) pluck(x,1))

# Inference of networks -------------------

# methods and parameters are in the options at the beginning of the script

if(keep_time) tic("Microbial association network inference")
if(verbose_output) cat("Please be patient, this will take a while...\n")
# list for results, 1 slot for each object
MAN_inf_results <- vector("list", length = length(physeq_list_filt_5))
# create list for methods, 1 slot for each method
inf_meth_list <- vector("list", length = length(inf_methods))

for(i in seq_along(physeq_list_filt_5)){
  name <- names(physeq_list_filt_5)[i]
  if(verbose_output) cat("\n", "Inferring network(s) for ", name, ", ", i, "of ", 
                         length(physeq_list_filt_4), "\n", sep=" ")
  if(keep_time) {
    ticmessage_object <- paste("microbial association network inference, ",
                               name, ", ", i, " of ", length(physeq_list_filt_4), sep = "")
    tic(ticmessage_object)
  }
  for(j in seq_along(inf_methods)){
    if(keep_time){
      ticmessage_method <- paste("\n", "inference with method ",
                                 inf_methods[j], ", ", j, " of ", length(inf_methods), sep = "")
      tic(ticmessage_method)
    }
    infmethod <- inf_methods[[j]]
    infparam <- inf_methods_param[[j]]
    
    # inference is carried out here using a user defined function loaded with the `source()` command 
    # you do not need to provide much detail, everything is taken care of in the options
    # of course could be customized in the function call

    inf_meth_list[[j]] <- infer_MAN(myphyseq = physeq_list_filt_5[[i]],
                                    inf_method = infmethod,
                                    method_parameters = infparam)
    names(inf_meth_list)[j] <- infmethod
    if(keep_time) toc()
  }
  MAN_inf_results[[i]] <- inf_meth_list
  names(MAN_inf_results)[i] <- name
  if(keep_time) toc()
}


# put together a report with classes of the objects
inference_report <- vector("list", length(MAN_inf_results))
for(i in seq_along(MAN_inf_results)){
  inference_report[[i]] <- map(MAN_inf_results[[i]], \(x) data.frame(obj_class = class(x))) |>
      list_rbind(names_to = "inf_method")
    names(inference_report)[i] <- names(MAN_inf_results)[i]
  }
  
inference_report_df <- bind_rows(inference_report, .id = "dataset")

# save the list and the data frame and do some clean-up

save(MAN_inf_results, 
     file = file.path(getwd(), output_folder, str_c(fn_pref_n,"MAN_inf_res.Rdata")))
write_tsv(inference_report_df,
          file = file.path(getwd(), output_folder, str_c(fn_pref_n,"MAN_inf_rep.txt")))

rm(inf_meth_list, inference_report)
if(play_audio) beep(sound = 6)
if(keep_time) toc()

# analyze networks and extract properties ---------------------------------

# calculate network statistics with netAnalyze

if(keep_time) tic("Calculate network statistics")

net_stats <- vector("list", length = length(MAN_inf_results))
# net_stat_results is a list, do the calculation for each of the datasets, all inference methods
for(i in seq_along(MAN_inf_results)){
  name <- names(MAN_inf_results)[i]
  if(verbose_output) cat("Calculating network stats for",name,"\n")
  net_stats[[i]] <- map(MAN_inf_results[[i]], calculate_net_stats)
  names(net_stats)[i] <- names(MAN_inf_results)[i]
}

netstat_report <- vector("list", length(net_stats))
for(i in seq_along(net_stats)){
  netstat_report[[i]] <- map(net_stats[[i]], \(x) data.frame(obj_class = class(x))) |>
    list_rbind(names_to = "inf_method")
  names(netstat_report)[i] <- names(net_stats)[i]
}

netstat_report_df <- bind_rows(netstat_report, .id = "dataset")
write_tsv(netstat_report_df,
          file = file.path(getwd(), output_folder, str_c(fn_pref_n,"netstat_rep.txt")))
rm(netstat_report)
if(keep_time) toc()


# extracting global network properties ------------------------------------

if(keep_time) tic("Extraction global network properties")
# extracting the global network stats 
global_props_list_a <-vector("list",length(net_stats))
global_props_list_l <-vector("list",length(net_stats)) 
for (i in seq_along(net_stats)){
  dataset <- names(net_stats)[i]
  mynetstats <- net_stats[[i]]
  # which are microNetProps
  which_microNet_props <- map_lgl(mynetstats, \(x) class(x) == "microNetProps")
  n_microNet_props <- sum(which_microNet_props)
  if(n_microNet_props==0){
    next
  } else {
    mynetstats <- mynetstats[which_microNet_props]
    global_props_list_all <-vector("list",length(mynetstats))
    global_props_list_lcc <-vector("list",length(mynetstats))
    for (j in seq_along(mynetstats)){
      nnodes <- sum(rowSums(mynetstats[[j]]$input$assoMat1 != 0)>1)
      ntaxa <- nrow(mynetstats[[j]]$input$assoMat1)
      nposedge <- sum(mynetstats[[j]]$input$assoMat1[lower.tri(mynetstats[[j]]$input$assoMat1)]>0)
      nnegedge <- sum(mynetstats[[j]]$input$assoMat1[lower.tri(mynetstats[[j]]$input$assoMat1)]<0)
      extra_prop_vector <- c(nnodes, ntaxa, nposedge, nnegedge)
      names(extra_prop_vector) <- c("nnodes", "ntaxa", "nposedge", "nnegedge")
      global_props_list_all[[j]] <- c(unlist(mynetstats[[j]]$globalProps),extra_prop_vector)
      global_props_list_lcc[[j]] <- c(unlist(mynetstats[[j]]$globalPropsLCC),extra_prop_vector)
      names(global_props_list_all)[j] <- names(global_props_list_lcc)[j]<- names(mynetstats[j])
    }
    # create data frame with results
    all_df <- bind_rows(global_props_list_all, .id = "method")
    lcc_df <- bind_rows(global_props_list_lcc, .id = "method") 
    global_props_list_a[[i]] <- all_df
    global_props_list_l[[i]] <- lcc_df
    names(global_props_list_a)[i] <- names(global_props_list_l)[i] <- names(net_stats)[i]
  }
}
# note that str_sub only works if you have <=9 inference methods in inf_methods
global_all_df <- bind_rows(global_props_list_a, .id = "dataset") 
global_lcc_df <- bind_rows(global_props_list_l, .id = "dataset") 

global_all_df <- left_join(global_all_df, 
                           study_metadata) |>
  mutate(across(where(is.numeric), ~na_if(.,Inf)))
global_lcc_df <- left_join(global_lcc_df, 
                           study_metadata) |>
  mutate(across(where(is.numeric), ~na_if(.,Inf)))

# both might contain Inf which are replaced by NA (using dplyr::na_if()) to be handled
# correctly if doing PCA by pairwise deletion.
# the problem only occurs in avPath1 and clustCoef1, but I am handling it with a scoped mutate
# a better solution might be the use of hablar::rationalize()
rm(mynetstats)


# save the data frames for further use
write_tsv(global_all_df, 
          file = file.path(getwd(), output_folder, str_c(fn_pref_n,"_netpropall.txt",sep="")))
write_tsv(global_lcc_df, 
          file = file.path(getwd(), output_folder, str_c(fn_pref_n, "_netproplcc.txt",sep="")))

# print a summary table
global_all_df

rm(global_props_list_a, global_props_list_l, global_props_list_all, 
   global_props_list_lcc, all_df, lcc_df, nnodes, ntaxa, nposedge, nnegedge, extra_prop_vector)
if(play_audio) beep(sound = 6)
if(keep_time) toc()

# extract node properties -------------------------------------------------

# extract node properties
# using a loop takes slightly longer that using functionals but handles names better
# note that when using ASVs comparing nodes between datasets does not make much sense

# get prev_ab info from the saved files, choosing only those
# belonging to this iteration

this_ite_studies <- names(physeq_list)

prev_ab_files <- str_c(this_ite_studies,"_", prev_filter,"_prevabt.txt")

prev_ab_file_studies <- str_split_i(prev_ab_files, str_c("_", prev_filter),1)
prev_ab_file_studies <- prev_ab_file_studies[which(prev_ab_file_studies %in% this_ite_studies)]
prev_ab_file_studies <- prev_ab_file_studies[order(match(prev_ab_file_studies, this_ite_studies))] 


# for some reason need to add the wd on the MacBook not on the iMac
# likely due to different location of the Dropbox files
prev_ab_files_path <- file.path(getwd(),
                                "output",
                                "MAN_F",
                                prev_ab_files)
node_stat_list <- map(prev_ab_files_path, read.delim)
names(node_stat_list) <- names(physeq_list)

if(keep_time) tic("Extracting node properties")

node_stats <- vector("list", length = length(net_stats))
for (i in seq_along(net_stats)) {
  node_properties <- node_stat_list[[i]]
  
  if(verbose_output) cat("extracting node stats for", names(net_stats)[i],"\n")
  for (j in seq_along(net_stats[[i]])) {
    if (class(net_stats[[i]][[j]]) == "microNetProps") {
      node_stats[[i]][[j]] <- extract_node_stats(net_stat_list = net_stats[[i]][[j]], 
                                                 nodestat = node_properties)
      method <- names(net_stats[[i]])[j]
      dataset <- names(net_stats)[i]
      nrows <- nrow(node_stats[[i]][[j]])
      node_stats[[i]][[j]] <- bind_cols(
        dataset = rep(dataset,nrows),
        method = rep(method,nrows),
        node_stats[[i]][[j]]
      )
    } else {
      cat("no node stats to return for",
          names(net_stats)[i],
          names(node_stats[[i]])[j],
          "\n")
      next
    }
  }
  node_stats[[i]]<-bind_rows(node_stats[[i]])
}
node_stats_df <- bind_rows(node_stats)
# perform some tidying
node_stats_df <- node_stats_df |>
  tidyr::separate(dataset, into = c("Study", "nTarg"), sep = "_", remove = F)

# merge colors for phyla

if(exists("phylum_palette")){
  node_stats_df <- node_stats_df |> 
    mutate(phylumC = ifelse(phylum %in% phylum_palette$phylum, phylum, "Other")) |>
    left_join(rename(phylum_palette, phylumC = phylum)) %>%
    mutate(phylumC = as.factor(phylumC)) 
}
node_stats_df$phylumC <- fct_relevel(node_stats_df$phylumC, "Other", after = Inf)

write_tsv(node_stats_df, 
          file.path(output_folder, str_c(fn_pref_n, "_nodestats_df.txt",sep="")))
rm(node_stats, method)
if(play_audio) beep(sound = 6)
if(keep_time) toc()


# build tidygraphs  --------------------------

if(keep_time) tic("List with tidygraph objects created")

# creating a tidygraph object for each net
tidygraph_list <- vector("list", length = length(MAN_inf_results))

for(i in seq_along(MAN_inf_results)){
  if(verbose_output) cat("Converting in tidygraphs for", names(MAN_inf_results)[i],"\n")
  # the warning, if any, is not very informative, should consider passing names of datasets and methods
  tidygraph_list[[i]] <- MAN_inf_results[[i]] %>% map(microNet_to_tidygraph, fail_w_err = F, use_asso_matrix = T)
}
names(tidygraph_list)<-names(MAN_inf_results)
# make a report
tidygraph_list_report <-vector(mode = "list", length =length(tidygraph_list))
for(i in seq_along(tidygraph_list)){
  tidygraph_list_report[[i]] <- purrr::map_dfr(tidygraph_list[[i]], \(x) class(x)[1])
  names(tidygraph_list_report)[i]<-names(tidygraph_list)[i]
}
tidygraph_list_report <- list_rbind(tidygraph_list_report, names_to = "dataset")
tidygraph_list_report <- tidygraph_list_report |>
  pivot_longer(cols = spieceasi:ccrepe, names_to = "method", values_to = "obj_returned")
# save the report
write_tsv(tidygraph_list_report, 
          file.path(str_c(fn_pref_n, "_tidygraphreport.txt",sep=""))) 

if(keep_time) toc()


# merge further node stats, extract edges ------------------------------------------------

# optionally merge further node stats (depends on merge_n_stats) 
# and calculate edge betweenness (depends on calc_e_betw)
# I am using a loop

if(keep_time) tic("Stats added to tidygraphs, edge dataframe created")
tidygraph_list_wstats <- vector("list", length = length(tidygraph_list))
# the list with the edge data frames
edge_list <- vector("list", length = length(tidygraph_list))

for(i in seq_along(tidygraph_list)){
  # need to be reset
  inner_tgl <- vector("list", length = length(tidygraph_list[[i]]))
  inner_el <- vector("list", length = length(tidygraph_list[[i]]))
  dtst <- names(tidygraph_list)[i]
  for(j in seq_along(inner_tgl)){
    if(is.tbl_graph(tidygraph_list[[i]][[j]])){
      mthd = names(tidygraph_list[[i]])[j]
      nstats <- node_stats_df %>% dplyr::filter(dataset == dtst & method == mthd)
      inner_tgl[[j]] <- merge_stats(tg = tidygraph_list[[i]][[j]],
                                    node_stats = nstats, 
                                    ebetw = calc_e_betw)
      # extract edge tibble
      inner_el[[j]] <- inner_tgl[[j]] %>% 
        activate(edges) %>% 
        as_tibble() %>%
        mutate(method = mthd)
      # do naming
      names(inner_tgl)[j] <- names(inner_el)[j] <- names(tidygraph_list[[i]])[j]
    } else {
      next
    }
  }
  tidygraph_list_wstats[[i]] <- inner_tgl
  edge_list[[i]] <- bind_rows(inner_el)
  if(nrow(edge_list[[i]])==0){
    edge_list[[i]] <- NULL
  } else {
    edge_list[[i]] <- edge_list[[i]] %>% 
      mutate(dataset = dtst) %>% dplyr::select(dataset, method, !(dataset:method))
  }
  names(tidygraph_list_wstats)[i] <- names(edge_list)[i] <- names(tidygraph_list)[i]
}
rm(nstats)

# build and fix the edge df
edge_list_df <- bind_rows(edge_list)


edge_list_df <- edge_list_df %>% 
  mutate(name_from_sorted = if_else(from_name < to_name, from_name, to_name),
         name_to_sorted = if_else(from_name < to_name, to_name, from_name)) %>%
  mutate(edge_name = str_c(name_from_sorted, name_to_sorted, sep = "--"))

# save the tidygraphs
saveRDS(tidygraph_list_wstats, 
        file =file.path(getwd(),
                        output_folder,
                        str_c(fn_pref_n, "tidygraphswstats.rds",sep="_")))

write_tsv(edge_list_df, 
          file = file.path(getwd(),
                           output_folder,
                           str_c(fn_pref_n, fn_pref_n, "edgelist_df.txt",sep="_")))
if(keep_time) toc()

# make Venn plots with edges ----------------------------------------------

# make Venn plots

# I am using a loop
if(do_Venn){
  if(keep_time) tic("Venn plots created and saved")
  Venn_list <- vector("list", length(tidygraph_list))
  for(i in seq_along(tidygraph_list)){
    # need to select only elements of class tbl_graph
    tblgrphs <- map_lgl(tidygraph_list[[i]], is.tbl_graph)
    names(Venn_list) <- names(tidygraph_list)
    if(length(tidygraph_list[[i]][tblgrphs])>1){
      inner_list <- map(tidygraph_list[[i]][tblgrphs], function(x) as_ids(E(x)))
      Venn_title <- names(tidygraph_list)[i]
      Venn_file <- file.path(getwd(), 
                             output_folder,
                             str_c(fn_pref_n, Venn_title, "venns.tiff",sep="_"))
      my_fill <- (2:5)[1:length(tidygraph_list[[i]][tblgrphs])]
      Venn_list[[i]] <- venn.diagram(inner_list, 
                                     fill = my_fill, 
                                     alpha = 0.3, 
                                     filename = Venn_file, 
                                     margin = 0.05, 
                                     main = Venn_title, 
                                     main.cex = 1.5, 
                                     main.fontface = "bold", 
                                     main.fontfamily = "sans", 
                                     main.pos = c(0.5, 1.05)
      )
    }
  }
  rm(inner_list, Venn_title, Venn_file)
  if(keep_time) toc()
}

if(play_audio) beep(sound = 6)


# plot tidygraphs ---------------------------------------------------------

if(keep_time) tic("Plotting the networks with ggraph")
# create a list of network plots

edge_col_scheme <- c(copres = "cyan", mut_ex = "magenta")

netplot_list <- vector("list", length = length(tidygraph_list_wstats))
for(i in seq_along(tidygraph_list_wstats)){
  netplot_list_2 <- vector("list", length = length(tidygraph_list_wstats[[i]]))
  dtst <- names(tidygraph_list_wstats)[i]
  for(j in seq_along(tidygraph_list_wstats[[i]])){
    if(!is.tbl_graph(tidygraph_list_wstats[[i]][[j]])){
      netplot_list_2[[j]] <- "no plot to return"
      next
    } else {
      tg <- tidygraph_list_wstats[[i]][[j]]
      mthd <- names(tidygraph_list_wstats[[i]])[j]
      # the default for the argument c0l0r is phylum, otherwise
      # use c0l0r = clust_memb or use c0l0r = phylumC for custom colors
      # argument lp = "bottom" is the default; "right" is an alternative
      netplot_list_2[[j]] <- plot_ggraph(tidy_graph = tg,
                                         method = mthd,
                                         stname = dtst,
                                         lp = "right",
                                         c0l0r = "phylumC",
                                         edge_cols = edge_col_scheme)
      names(netplot_list_2)[j] <- mthd
      print(netplot_list_2[[j]])
    }
  }
  netplot_list[[i]] <- netplot_list_2
  names(netplot_list)[i] <- dtst
}

# save tidygraph list
save(netplot_list, file = file.path(getwd(), 
                                    output_folder,
                                    str_c(fn_pref_n, "tdglist.Rdata",sep="_")))



rm(tg, dtst, mthd, netplot_list_2)
if(play_audio) beep(sound = 6)
if(keep_time) toc() # closing tidygraph plot

if(keep_time) toc() # closes batch operation


# This work was was carried out within the PRIN 2022 Project METAOlive 
# 2022NN28ZZ and received funding from the European Union Next-GenerationEU, 
# CUP C53D23007560001 (PIANO NAZIONALE DI RIPRESA E RESILIENZA (PNRR) 
# – MISSIONE 4 COMPONENTE 2,  INVESTIMENTO 1.4 – D.D. 1048 14/07/2023). 
# This script and its contents reflects only the authors’ views and opinions,  
# neither the European Union nor the European Commission can be considered
# responsible for them.

# MIT License

# Copyright (c) 2024 EugenioP

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
  
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
