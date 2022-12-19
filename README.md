# GM_trop_temp
This repository contains scripts for analyses related to the study "Gut microbiota variation between climatic zones and due to migration strategy in passerine birds"

01.demultiplexing.txt - Bash commands for demultiplexing fastq files provided by the sequencing facility.

02.dada2_denoising.R - Quality filtering and denoising of demultiplexed fastq files with data2 

03.preprare_phyloseq_database.Rmd - merges the dada2 outputs of the two sequencing runs into a single phyloseq object.

04.phylogeny_taxonomy_etc.R - assigns the taxonomy and creates a bacterial phylogeny (+ some other phyloseq customizations).

05.remove.contamianats_OTU.clustering.R - this script uses the decontam package to identify and eliminate contaminating ASVs. The ASVs are then clustered into OTUs using vsearch and 97% sequence similarity 

06.statistics.Rmd - scripts used for statistical analysis 

PH_final.R - the phyloseq database used in step 06.
