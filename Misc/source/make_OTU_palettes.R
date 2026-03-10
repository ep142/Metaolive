# make_OTU_palettes
# well this is a misnomer really, should be feature, not OTU
# takes an OTU_matrix, taxa_metadata, the lowest level of aggregation of taxa
# and a flag (override_OTU, set to F if you want more than 15 categories in palettes)

make_OTU_palettes <- function(OTUmatrix, 
                              taxa_metadata = input_data$taxa_metadata,
                              tax_agg = input_data$tax_agg,
                              override_OTU = overrideOTU){
  ET3 <- taxa_metadata %>%
    dplyr::filter(label %in% colnames(OTUmatrix))
  # must fill the NAs with info from higher taxa
  ET3 <- ET3 %>% 
    mutate(domain = ifelse(is.na(domain),"Other", domain)) %>%
    mutate(phylum = ifelse(is.na(phylum), domain, phylum)) %>%
    mutate(class = ifelse(is.na(class), phylum, class)) %>%
    mutate(order = ifelse(is.na(order), class, order)) %>%
    mutate(family = ifelse(is.na(family), order, family))
  
  ET3$class <- as.factor(ET3$class)
  if (tax_agg == "class"){
    ET3$family <- ET3$class
  } else {
    ET3$family <- as.factor(ET3$family)
  }
  ncolsfamily <- nlevels(ET3$family)
  ncolsclass <- nlevels(ET3$class)
  familypalette <- if(ncolsfamily<=12) {
    brewer.pal(ncolsfamily, "Set3")[1:ncolsfamily]
  } else {
    distinctColorPalette(ncolsfamily)
  }
  classpalette <- if(ncolsclass<=12) {
    brewer.pal(ncolsclass, "Set3")[1:ncolsclass]
  } else {
    distinctColorPalette(ncolsclass)
  }
  OTUcolors <- ET3[, c("label", "class", "family")]
  fcolors <- data.frame(family = levels(OTUcolors$family),
                        fcolor = familypalette)
  ccolors <- data.frame(class = levels(OTUcolors$class),
                        ccolor = classpalette)
  
  # Merge the palettes with label and lineage info
  OTUcolors <- merge(OTUcolors, ccolors, by = "class")
  OTUcolors <- merge(OTUcolors, fcolors, by = "family")
  OTUcolors <- arrange(OTUcolors, label)
  # set the palette to use
  catOTUs <- "class"
  if (input_data$tax_agg != "class"){
    catOTUs <- "family"
    familystring <- paste(" and ", nlevels(OTUcolors$family), 
                          " families.",sep = "")
    if (!override_OTU) {
      catOTUs <- ifelse(nlevels(OTUcolors$family)<=15,
                        "family","class")
    }
    if (catOTUs == "class") {
      myOTUpalette <- as.character(OTUcolors$ccolor)
      myOTUtext <- as.character(ccolors$class)
      myOTUcolors <- as.character(ccolors$ccolor)
    } else {
      myOTUpalette <- as.character(OTUcolors$fcolor)
      myOTUtext <- as.character(fcolors$family)
      myOTUcolors <- as.character(fcolors$fcolor)
    }
  }
  OTU_palette_list <- list(
    palette = myOTUpalette,
    text = myOTUtext,
    colors = myOTUcolors,
    OTU_colors_df = OTUcolors,
    cat_OTU_flag = catOTUs
  )
  return(OTU_palette_list)
}