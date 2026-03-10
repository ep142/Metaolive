# Metaolive

This repository hosts data and scripts created by me (Eugenio Parente) for 
the "METAOLIVE: meta-omic approaches for traditional table olive fermentation" 
project (more [here](http://web.unibas.it/parente/?page_id=1945)). 
This project is co-funded by Ministero dell’Università e della Ricerca PRIN 2022, 
proposal 2022NN28ZZ, and received funding from the European Union 
Next-GenerationEU - PIANO NAZIONALE DI RIPRESA E RESILIENZA (PNRR).  

* olive_FMBN folder contains scripts and data with phyloseq objects (and their metadata; 
including an extra_metadata file with further metadata for bacterial and fungal samples) 
on olive fermentations extracted from [FoodMicrobionet v5](https://github.com/ep142/FoodMicrobionet).  

* the Misc folder contains miscellaneous scripts used to process data generated in the project:  

  + skan_it_growth_curve_v0_4_1.Rmd is a Rmarkdown file designed to process in 
  batch kinetics files generated with the [SkanIt](https://www.thermofisher.com/it/en/home/life-science/lab-equipment/microplate-instruments/plate-readers/software.html) software, fit a number of growth 
  curve models and select the best one in terms of AIC. The current version also fits a secondary model. An example growth curve 
  file is provided. I am pretty sure it can be adapted to other plate readers.  
  
  + FMBN_metaolive_stats R markdown documents contain the code for descriptive statistical analysis on OliveFMBN  
  
  + METAOlive_MAN_ite_1_F_genus_v2_3.R and METAOlive_postprocess_f_genus.Rmd, together with the Source folder containing accessory scripts are examples of the pipeline used for inference and analysis of Microbial Association Networks  

Data for bacteria have been used in a poster presented at the 
[Food Systems Microbiomes](https://foodsystemsmicrobiomes.org) conference,
an oral communication at the AJLAB meeting, on 9/8/2024, a poster presented at 
[MD2025](https://www.md25.simtrea.org) and for a recent publication:  
Ricciardi, A., López, F.N.A., Giavalisco, M., Pietrafesa, R., Parente, E., 2025. 
Determining the core bacterial and fungal genera in table olive fermentations. Int. J. Food Microbiol. 442, 111344. https://doi.org/10.1016/j.ijfoodmicro.2025.111344  

Please cite as:
Ricciardi, A., López, F.N.A., Giavalisco, M., Pietrafesa, R., Parente, E., 2025. 
Determining the core bacterial and fungal genera in table olive fermentations. Int. J. Food Microbiol. 442, 111344. https://doi.org/10.1016/j.ijfoodmicro.2025.111344  