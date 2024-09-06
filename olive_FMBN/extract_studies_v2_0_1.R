# Olive_extract_studies v2_0_1

# this version was run on FMBN 5_0_1 on 6/9/2024
# REMINDER, EACH TIME I NEED TO MAKE SURE THAT THE PHYLOSEQ LISTS INCLUDE ALL STUDIES

# a script designed to extract studies from FoodMicrobionet 5 for the
# "METAOlive" project


#  notes for self ---------------------------------------------------------

# I could modify the function to create separate phyloseqs for brine or fruit, but
# there is little point

# Install/load packages ---------------------------------------------------
# install packages, set general parameters and options

.bioc_packages <- c("BiocManager", "phyloseq")
.cran_packages <- c("tidyverse",  "data.table", "tictoc", "beepr", "logr", "crayon")

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

# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

# other setup operations
opar <- par(no.readonly=TRUE) 
par(ask=F) 
# set.seed(1234) 
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


# functions ---------------------------------------------------------------

# make_phyloseq -----------------------------------------------------------

# samples, edges and taxa are taken from FMBN 
# (possibly with modifications and/or filtering)
# samples_to_process is a vector of sampleIds to process in a phyloseq object
# study_name is the name of the study, i.e. the name of the slot in the list
# the function will produce between 1 and 5 (all, samples, cheese and/or milk 
# samples, environments, subject to the condition that there are at least 15 
# samples in each type) phyloseq objects
# taxa are pruned as appropriate

make_phyloseqs <- function(samples, edges, taxa, samples_to_process, study_name, output_folder){
  # get edges
  myedges <- edges %>% dplyr::filter(sampleId %in% samples_to_process) 
  mysamples <- samples %>% dplyr::filter(sampleId %in% samples_to_process)
  mytaxa <- taxa %>% dplyr::filter(taxonId %in% unique(myedges$taxonId))
  myedges <- left_join(myedges, dplyr::select(mysamples, sampleId, SRA_run, n_reads2)) %>%
    mutate(weight = round(weight/100*n_reads2))
  myedges <- left_join(myedges, dplyr::select(mytaxa, taxonId, label))
  myOTUtable <- myedges %>% dplyr::select(SRA_run, t_label = label, weight) %>%
    pivot_wider(names_from = SRA_run, values_from = weight, values_fill = 0) %>%
    column_to_rownames(var = "t_label")
  myOTUtable <- as.matrix(myOTUtable)
  mytaxa <- mytaxa %>% select(2:9) %>% column_to_rownames(var = "label") %>% as.matrix()
  mysamples <- mysamples %>% column_to_rownames(var = "SRA_run")
  myphseq <- phyloseq(otu_table(myOTUtable, taxa_are_rows = T), sample_data(mysamples),
                      tax_table(mytaxa))
  # save phyloseq objects for the food samples only, environmental samples only, 
  # whole phyloseq
  saveRDS(myphseq, file = file.path(output_folder,
                                    str_c(study_name,"rds", sep = ".")))
  L2_df <- data.frame(myphseq@sam_data) %>% dplyr::select(sampleId, s_type, L2, L3, nature, process)
  
  # food and environmental samples (I have tried a loop and it does not work
  # apparently the problem is due to subset_samples, which has a problem 
  # when working inside a loop inside a function)
  if(length(unique(mysamples$s_type))>1){
    if(nrow(dplyr::filter(mysamples, s_type == "Environment"))>=15){
      phy <- subset_samples(myphseq, s_type == "Environment")
      phy <- prune_taxa(taxa_sums(phy) > 0, phy)
      phyfn <- file.path(output_folder, str_c(study_name, "_environment.rds"))
      saveRDS(phy, file = phyfn)
    }
    if(nrow(dplyr::filter(mysamples, s_type == "Sample"))>=15){
      phy <- subset_samples(myphseq, s_type == "Sample")
      phy <- prune_taxa(taxa_sums(phy) > 0, phy)
      phyfn <- file.path(output_folder, str_c(study_name, "_sample.rds"))
      saveRDS(phy, file = phyfn)
    }
    }
}

# load FMBN ---------------------------------------------------------------

# I am using FMBN_plus, version 5_0_1, in a folder named FMBN
FMBN_path <- file.path("FMBN", "FMBN_plus.rds") 
FMBN <- readRDS(FMBN_path)


# filter studies ----------------------------------------------------------
studies <- FMBN$studies

studies_filt <- studies %>% 
  dplyr::filter(str_detect(bioinf_software, "R dada2")) %>%
  dplyr::filter(str_detect(short_descr, "olive|olives|Olive|Olives") & food_group == "Vegetables")

# 17 studies remaining
nrow(studies_filt)

# save the studies table
write_tsv(studies_filt, "FMBN_5_0_1_olive_studies.txt")


# get edges, samples, taxa for bacteria ------------------------------------------------

# the samples
samples_B <- FMBN$samples_B %>% 
  dplyr::filter(studyId %chin% studies_filt$studyId) %>%
  dplyr::filter(!is.na(SRA_run))

# remove any blank or mock
samples_B <- samples_B %>% 
  dplyr::filter(!str_detect(L1, "blank|Blank|mock|Mock"))
# add info from FoodEx2 + fix the issue due to mismatching L6
samples_B <- left_join(samples_B, dplyr::select(FMBN$foodex2exp, L2, L3, L6), by = "L6") 

# the edges
edges_B <- FMBN$edges_B %>% 
  dplyr::filter(sampleId %in% samples_B$sampleId)

# the taxa
taxa_B <- FMBN$taxa %>% 
  dplyr::filter(taxonId %in% unique(edges_B$taxonId))

# get edges, samples, taxa for fungi ------------------------------------------------

# the samples
samples_F <- FMBN$samples_F %>% 
  dplyr::filter(studyId %chin% studies_filt$studyId) %>%
  dplyr::filter(!is.na(SRA_run))

# remove any blank or mock
samples_F <- samples_F %>% 
  dplyr::filter(!str_detect(L1, "blank|Blank|mock|Mock"))
# add info from FoodEx2 + fix the issue due to mismatching L6
samples_F <- left_join(samples_F, dplyr::select(FMBN$foodex2exp, L2, L3, L6), by = "L6") 

# the edges
edges_F <- FMBN$edges_F %>% 
  dplyr::filter(sampleId %in% samples_F$sampleId)

# the taxa
taxa_F <- FMBN$taxa %>% 
  dplyr::filter(taxonId %in% unique(edges_F$taxonId))


# create phyloseqs --------------------------------------------------------


# check if the phyloseq folders exist, if not create it
if (!dir.exists("physeqs_B")) dir.create("physeqs_B")
if (!dir.exists("physeqs_F")) dir.create("physeqs_F")


# phyloseqs for bacteria --------------------------------------------------

# optionally start logging
if(debug_mode) selsampleslog_B <- log_open(file_name = file.path("debugging_logs","selsamples_B.log"))

# loop over the selected studies, putting the samples to process in a phyloseq in a list
samples_to_physeq_B <- vector(mode = "list")
# filter the studies
studies_filt_B <- studies_filt |>
  dplyr::filter(str_detect(target, "16S"))

for(i in seq_along(studies_filt_B$study)){
  if(debug_mode) log_print(str_c("processing study", i, sep = " "), console = F)
  # in a very few cases need to remove samples which do not have a SRA accession
  samples_temp <- samples_B %>% dplyr::filter(studyId == studies_filt_B$studyId[i])
  # if more than one target, split (only subgroups with ≥15 samples will be kept)
  samples_per_target <- samples_temp %>% 
    group_by(target1) %>% 
    summarise(n_samples = n())  
  if(debug_mode) log_print(str_c("detected", 
                                   length(unique(samples_temp$target1)),
                                   "targets for study", i, sep = " "), 
                             console = F)
    for(j in seq_along(samples_per_target$target1)){
      if(debug_mode) log_print(str_c("processing", 
                                     samples_per_target$target1[j],
                                     "for study", i, sep = " "), 
                               console = F)
      samples_to_add <- samples_temp %>% 
        dplyr::filter(target1 == samples_per_target$target1[j]) %>%
        pull(sampleId)
      list_pos <- length(samples_to_physeq_B)
      samples_to_physeq_B[[list_pos+1]] <- samples_to_add
      target_region <- str_sub(samples_per_target$target1[j], 5)
      names(samples_to_physeq_B)[list_pos+1] <- str_c(studies_filt_B$studyId[i],
                                                  target_region, sep = "_")
    }
  }
rm(samples_temp, samples_per_target, samples_to_add, list_pos, target_region, i, j)
# optionally, end logging
if(debug_mode){
  # Close the log
  log_close()
  # View log
  if(verbose_output) writeLines(readLines(selsampleslog_B, encoding = "UTF-8"))
}


# create and save phyloseqs for bacteria
if(keep_time) tic("\nCreating and saving phyloseqs for bacteria")
for(i in seq_along(samples_to_physeq_B)){
  make_phyloseqs(samples = samples_B, edges = edges_B, taxa = taxa_B,
                 samples_to_proces = samples_to_physeq_B[[i]],
                 study_name = names(samples_to_physeq_B)[i],
                 output_folder = "physeqs_B")
}
if(play_audio) beep(sound = 6)
if(keep_time) toc()

# phyloseqs for fungi --------------------------------------------------

# optionally start logging
if(debug_mode) selsampleslog_F <- log_open(file_name = file.path("debugging_logs","selsamples_F.log"))

# loop over the selected studies, putting the samples to process in a phyloseq in a list
samples_to_physeq_F <- vector(mode = "list")
# filter the studies
studies_filt_F <- studies_filt |>
  dplyr::filter(str_detect(target, "ITS"))

for(i in seq_along(studies_filt_F$study)){
  if(debug_mode) log_print(str_c("processing study", i, sep = " "), console = F)
  # in a very few cases need to remove samples which do not have a SRA accession
  samples_temp <- samples_F %>% dplyr::filter(studyId == studies_filt_F$studyId[i])
  # if more than one target, split (only subgroups with ≥15 samples will be kept)
  samples_per_target <- samples_temp %>% 
    group_by(target1) %>% 
    summarise(n_samples = n())  
  if(debug_mode) log_print(str_c("detected", 
                                 length(unique(samples_temp$target1)),
                                 "targets for study", i, sep = " "), 
                           console = F)
  for(j in seq_along(samples_per_target$target1)){
    if(debug_mode) log_print(str_c("processing", 
                                   samples_per_target$target1[j],
                                   "for study", i, sep = " "), 
                             console = F)
    samples_to_add <- samples_temp %>% 
      dplyr::filter(target1 == samples_per_target$target1[j]) %>%
      pull(sampleId)
    list_pos <- length(samples_to_physeq_F)
    samples_to_physeq_F[[list_pos+1]] <- samples_to_add
    target_region <- str_sub(samples_per_target$target1[j], 5)
    names(samples_to_physeq_F)[list_pos+1] <- str_c(studies_filt_F$studyId[i],
                                                    target_region, sep = "_")
  }
}
rm(samples_temp, samples_per_target, samples_to_add, list_pos, target_region, i, j)
# optionally, end logging
if(debug_mode){
  # Close the log
  log_close()
  # View log
  if(verbose_output) writeLines(readLines(selsampleslog_F, encoding = "UTF-8"))
}


# create and save phyloseqs for fungi
if(keep_time) tic("\nCreating and saving phyloseqs for fungi")
for(i in seq_along(samples_to_physeq_F)){
  make_phyloseqs(samples = samples_F, edges = edges_F, taxa = taxa_F,
                 samples_to_proces = samples_to_physeq_F[[i]],
                 study_name = names(samples_to_physeq_F)[i],
                 output_folder = "physeqs_F")
}
if(play_audio) beep(sound = 6)
if(keep_time) toc()




# load and merge phyloseqs ------------------------------------------------
# bacteria
# I need to remove all the studies which have a suffix like sample, environment
all_physeqs_B <- list.files(path = file.path(".", "physeqs_B"))
physeqs_to_use_B <- !str_detect(all_physeqs_B, "environment|sample")
# the phyloseqs I need to use
olive_physeqs_B <- all_physeqs_B[physeqs_to_use_B]
study_ids <- str_split_i(olive_physeqs_B, "_", 1)
physeq_files_B <- map_chr(olive_physeqs_B, 
                       ~file.path(getwd(),"physeqs_B",.x))
names(physeq_files_B) <- olive_physeqs_B

physeq_list_B <- map(physeq_files_B, readRDS)

# need to remove suffix from names
names(physeq_list_B) <- str_split_i(names(physeq_list_B), "\\.",1)
# manually
merged_physeq_B <- merge_phyloseq(physeq_list_B[[1]], physeq_list_B[[2]], 
                                physeq_list_B[[3]], physeq_list_B[[4]], 
                                physeq_list_B[[5]], physeq_list_B[[6]], 
                                physeq_list_B[[7]], physeq_list_B[[8]],
                                physeq_list_B[[9]], physeq_list_B[[10]],
                                physeq_list_B[[11]], physeq_list_B[[12]], 
                                physeq_list_B[[13]], physeq_list_B[[14]])

# quick check
View(as(otu_table(merged_physeq_B), "matrix"))
View(as(sample_data(merged_physeq_B), "data.frame"))

# save
saveRDS(merged_physeq_B, file = "olive_FMBN_5_B_phy.rds")
sample_data <- as(sample_data(merged_physeq_B), "data.frame")
write_tsv(sample_data, file = "olive_FMBN_5_B_sample_data.txt")

# fungi
# I need to remove all the studies which have a suffix like sample, environment
all_physeqs_F <- list.files(path = file.path(".", "physeqs_F"))
physeqs_to_use_F <- !str_detect(all_physeqs_F, "environment|sample")
# the phyloseqs I need to use
olive_physeqs_F <- all_physeqs_F[physeqs_to_use_F]
study_ids <- str_split_i(olive_physeqs_F, "_", 1)
physeq_files_F <- map_chr(olive_physeqs_F, 
                          ~file.path(getwd(),"physeqs_F",.x))
names(physeq_files_F) <- olive_physeqs_F

physeq_list_F <- map(physeq_files_F, readRDS)

# need to remove suffix from names
names(physeq_list_F) <- str_split_i(names(physeq_list_F), "\\.",1)
# manually
merged_physeq_F <- merge_phyloseq(physeq_list_F[[1]], physeq_list_F[[2]], 
                                  physeq_list_F[[3]], physeq_list_F[[4]], 
                                  physeq_list_F[[5]], physeq_list_F[[6]], 
                                  physeq_list_F[[7]], physeq_list_F[[8]],
                                  physeq_list_F[[9]])

# quick check
View(as(otu_table(merged_physeq_F), "matrix"))
View(as(sample_data(merged_physeq_F), "data.frame"))

# save
saveRDS(merged_physeq_F, file = "olive_FMBN_5_F_phy.rds")
sample_data <- as(sample_data(merged_physeq_F), "data.frame")
write_tsv(sample_data, file = "olive_FMBN_5_F_sample_data.txt")


# Session report ----------------------------------------------------------
# This script has been tested with R 4.3.3 or R 4.4.0, x86 build, on a MacBook Pro
# 14 inch, 2023, with chip Apple M3 Pro and 16 GB RAM and on a MacStudio 2023
# with Apple M2 Max chip, both running MacOs Sonoma 14.5
sessionInfo()


# Copyright notice --------------------------------------------------------

# Script created by Eugenio Parente (eugenio.parente@unibas.it), 
# Università degli Studi della Basilicata, 2023, 2024 https://github.com/ep142
# Permission is hereby granted, free of charge, to any person obtaining a copy 
# of this software and associated documentation files (the \"Software\"), to 
# deal in the Software without restriction, including without limitation the 
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or 
# sell copies of the Software, and to permit persons to whom the Software is 
# furnished to do so, subject to the following conditions:  
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
# IN THE SOFTWARE.

# Acknowledgements --------------------------------------------------------

# This work was was carried out within the PRIN 2022 Project METAOlive 2022NN28ZZ
# and received funding from the European Union Next-GenerationEU,
# CUP C53D23005460006
# (PIANO NAZIONALE DI RIPRESA E RESILIENZA (PNRR) – MISSIONE 4 COMPONENTE 2, 
# INVESTIMENTO 1.4 – D.D. 1048 14/07/2023). 
# This script and its contents reflects only the authors’ views and opinions, 
# neither the European Union nor the European Commission can be considered 
# responsible for them.
