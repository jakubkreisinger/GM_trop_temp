#In this script, ASVs taxonomy is assigned using dada2 implementation of RDP classifier and Silva v. 138
#Chimeric ASVs are detected using uchime and eliminated
#Information on host diet is added
#Phylogenetic tree for bacterial ASVs is constructed using FastTree
#Consensus phylogeny for passerine hosts is generated

library(dada2)
library(phyloseq)
library(plyr)
library(ape)
library(phytools)
library(ShortRead)

load("/media/kreising/DATA/data/TROPICAL_TEMPERATE_Birds_2021/data/PHYLOSEQ_tropical12.R")
PHYLOSEQ_tropical12

###########
#Taxonomy##
###########

SEQS<-as.character(refseq(PHYLOSEQ_tropical12))
tax50<-assignTaxonomy(SEQS,refFasta = "~/DB/DADA2/silva_nr99_v138_train_set.fa.gz",minBoot = 50,multithread = TRUE,
                      verbose = TRUE)

save(tax50,file = "/media/kreising/DATA/data/TROPICAL_TEMPERATE_Birds_2021/data/tax50.R")
load("/media/kreising/DATA/data/TROPICAL_TEMPERATE_Birds_2021/data/tax50.R")

tax50<-tax_table(tax50)

PHYLOSEQ_tropical12_tax50<-merge_phyloseq(PHYLOSEQ_tropical12,tax50)

#########################################################
#Detection and elimination of chimeric ASVs with usearch#
#########################################################
writeFasta(refseq(PHYLOSEQ_tropical12_tax50),file = "/media/kreising/DATA/data/TROPICAL_TEMPERATE_Birds_2021/data/ref.fasta")
setwd("/media/kreising/DATA/data/TROPICAL_TEMPERATE_Birds_2021/data/")
system("~/software/bin/usearch7.0.1090_i86linux32 -uchime_ref ref.fasta -db ~/DB/gold.fasta -nonchimeras FASTA.uchime.fasta -strand plus")
NONCHIM<-readDNAStringSet("/media/kreising/DATA/data/TROPICAL_TEMPERATE_Birds_2021/data/FASTA.uchime.fasta")

FILTER<-taxa_names(PHYLOSEQ_tropical12_tax50)%in%names(NONCHIM)
PHYLOSEQ_tropical12_tax50<-prune_taxa(FILTER,PHYLOSEQ_tropical12_tax50)
save(PHYLOSEQ_tropical12_tax50,file = "/media/kreising/DATA/data/TROPICAL_TEMPERATE_Birds_2021/data/PHYLOSEQ_tropical12_tax50.R")

###########################
#Add diet data#############
###########################

POTRAVA<-read.delim("/media/kreising/DATA/data/BIRD_MICROBIOME/PASSERINES_COMPARATIVE/FIRST_RUN/COPHYLOGENY/SPECIES_DATA/BirdFuncDat.txt")

load("/media/kreising/DATA/data/TROPICAL_TEMPERATE_Birds_2021/data/PHYLOSEQ_tropical12_tax50.R")

sample_data(PHYLOSEQ_tropical12_tax50)$species<-gsub("Locustela","Locustella",sample_data(PHYLOSEQ_tropical12_tax50)$species)
sample_data(PHYLOSEQ_tropical12_tax50)$species<-gsub("Illadopsis_cleaveri_?","Illadopsis_cleaveri",sample_data(PHYLOSEQ_tropical12_tax50)$species)
sample_data(PHYLOSEQ_tropical12_tax50)$species<-gsub("[?]","",sample_data(PHYLOSEQ_tropical12_tax50)$species)
sample_data(PHYLOSEQ_tropical12_tax50)$species<-gsub("Hipolasis_icterina","Hippolais_icterina",sample_data(PHYLOSEQ_tropical12_tax50)$species)
sample_data(PHYLOSEQ_tropical12_tax50)$species<-gsub("Platysteira_cyanea","Dyaphorophyia_cyanea",sample_data(PHYLOSEQ_tropical12_tax50)$species)
sample_data(PHYLOSEQ_tropical12_tax50)$species<-gsub("Alethe_castanea","Alethe_diademata",sample_data(PHYLOSEQ_tropical12_tax50)$species)

sample_data(PHYLOSEQ_tropical12_tax80)$species<-gsub("Locustela","Locustella",sample_data(PHYLOSEQ_tropical12_tax50)$species)
sample_data(PHYLOSEQ_tropical12_tax80)$species<-gsub("Illadopsis_cleaveri_?","Illadopsis_cleaveri",sample_data(PHYLOSEQ_tropical12_tax50)$species)
sample_data(PHYLOSEQ_tropical12_tax80)$species<-gsub("[?]","",sample_data(PHYLOSEQ_tropical12_tax50)$species)
sample_data(PHYLOSEQ_tropical12_tax80)$species<-gsub("Hipolasis_icterina","Hippolais_icterina",sample_data(PHYLOSEQ_tropical12_tax50)$species)
sample_data(PHYLOSEQ_tropical12_tax80)$species<-gsub("Platysteira_cyanea","Dyaphorophyia_cyanea",sample_data(PHYLOSEQ_tropical12_tax50)$species)
sample_data(PHYLOSEQ_tropical12_tax80)$species<-gsub("Alethe_castanea","Alethe_diademata",sample_data(PHYLOSEQ_tropical12_tax50)$species)
SEQS_SPECIES<-unique(sample_data(PHYLOSEQ_tropical12_tax50)$species)


POTRAVA$Scientific<-gsub(" ","_",POTRAVA$Scientific)

POTRAVA$Scientific[POTRAVA$Scientific=="Parus_caeruleus"]<-"Cyanistes_caeruleus"
POTRAVA$Scientific[POTRAVA$Scientific=="Parus_palustris"]<-"Poecile_palustris"
POTRAVA$Scientific[POTRAVA$Scientific=="Parus_ater"]<-"Periparus_ater"
POTRAVA$Scientific[POTRAVA$Scientific=="Nectarinia_reichenowi"]<-"Cinnyris_reichenowi"
POTRAVA$Scientific[POTRAVA$Scientific=="Nectarinia_olivacea"]<-"Cyanomitra_olivacea"
POTRAVA$Scientific[POTRAVA$Scientific=="Platysteira_concreta"]<-"Dyaphorophyia_concreta"
POTRAVA$Scientific[POTRAVA$Scientific=="Platysteira_castanea"]<-"Dyaphorophyia_castanea"
POTRAVA$Scientific[POTRAVA$Scientific=="Platysteira_cyanea"]<-"Dyaphorophyia_cyanea"

#Proportion of carniv
CARNIV_HERBIV<-apply(POTRAVA[,10:14],1,sum)/100
POTRAVA.2<-data.frame(species=POTRAVA$Scientific,CARNIV_HERBIV)

SD<-data.frame(sample_data(PHYLOSEQ_tropical12_tax50))
SD<-join(SD,POTRAVA.2)
SD<-sample_data(SD)
sample_names(SD)<-sample_names(PHYLOSEQ_tropical12_tax50)
sample_data(PHYLOSEQ_tropical12_tax50)<-SD

#####################################
#Bacterial phylogeny#################
#####################################
library(phyloseq)
library(DECIPHER)
library(ShortRead)
library(ape)

# align ASVs sequences
DNA<-DNAStringSet(refseq(PHYLOSEQ_tropical12_tax50))
ALIGNED<-AlignSeqs(RNAStringSet(DNA),processors = 6,iterations = 10,refinements = 10)
writeFasta(DNAStringSet(ALIGNED)
           ,file = "/media/kreising/DATA/data/TROPICAL_TEMPERATE_Birds_2021/data/DECIPHER.aligned.fasta")

setwd("/media/kreising/DATA/data/TROPICAL_TEMPERATE_Birds_2021/data/")

# FastTree phylogeny
system("~/software/bin/FastTree -gtr -nt DECIPHER.aligned.fasta > FASTTREE2.tre")

# add phylogeny to the phyloseq
BACT_TREE<-read.tree("/media/kreising/DATA/data/TROPICAL_TEMPERATE_Birds_2021/data/FASTTREE2.tre")
#plot(BACT_TREE)
PHYLOSEQ_tropical12_tax50<-merge_phyloseq(PHYLOSEQ_tropical12_tax50,BACT_TREE)
PHYLOSEQ_tropical12_tax50.final<-PHYLOSEQ_tropical12_tax50

save(PHYLOSEQ_tropical12_tax50.final,
     file = "/media/kreising/DATA/data/TROPICAL_TEMPERATE_Birds_2021/data/PHYLOSEQ_tropical12_tax50.final.R")


###################################
#phylogenetic tree for hosts#######
###################################

library(phangorn)
TREE<-read.nexus("/media/kreising/DATA/data/TROPICAL_TEMPERATE_Birds_2021/data/output.nex")

# maximum clade credibility
TREE.mcc<-maxCladeCred(TREE,rooted=FALSE)

write.tree(TREE.mcc,file = "/media/kreising/DATA/data/TROPICAL_TEMPERATE_Birds_2021/data/TREE.mcc.tree")

# rename species in host phylogenetic tree (to make it compatible with phyloseq database)
TREE.mcc<-read.tree("/media/kreising/DATA/data/TROPICAL_TEMPERATE_Birds_2021/data/TREE.mcc.tree")
TREE.mcc.renamed<-TREE.mcc

TREE.mcc.renamed$tip.label[TREE.mcc.renamed$tip.label=="Parus_caeruleus"]<-"Cyanistes_caeruleus"
TREE.mcc.renamed$tip.label[TREE.mcc.renamed$tip.label=="Parus_ater"]<-"Periparus_ater"
TREE.mcc.renamed$tip.label[TREE.mcc.renamed$tip.label=="Parus_palustris"]<-"Poecile_palustris"
TREE.mcc.renamed$tip.label[TREE.mcc.renamed$tip.label=="Nectarinia_olivacea"]<-"Cyanomitra_olivacea"
TREE.mcc.renamed$tip.label[TREE.mcc.renamed$tip.label=="Platysteira_castanea"]<-"Dyaphorophyia_castanea"
TREE.mcc.renamed$tip.label[TREE.mcc.renamed$tip.label=="Platysteira_concreta"]<-"Dyaphorophyia_concreta"
TREE.mcc.renamed$tip.label[TREE.mcc.renamed$tip.label=="Platysteira_cyanea"]<-"Dyaphorophyia_cyanea"
TREE.mcc.renamed$tip.label[TREE.mcc.renamed$tip.label=="Nectarinia_reichenowi"]<-"Cinnyris_reichenowi"

write.tree(TREE.mcc.renamed,
           file = "/media/kreising/DATA/data/TROPICAL_TEMPERATE_Birds_2021/data/TREE.mcc.renamed.tree")

