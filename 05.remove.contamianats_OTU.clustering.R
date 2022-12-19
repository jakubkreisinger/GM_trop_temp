library(dada2)
library(phyloseq)
library(decontam)
library(ShortRead)

#in this section puttatively contaminating ASVs are identified and eliminated using decontam 
#decontam analyzes were conducted separately for the first and the second MiSeq run (assuming different sources of contamination)
#putatively contaminating ASVs identified based on [1] prevalence in controls vs. biological samples and [2] based in correlation with PCR output for biological samples
#In the last step, ASVs are cluster to OTUs based on 97% sequence similarity

####################################################
#[2] Load PHYLOSEQ including negative controls - First run##
####################################################
load("/media/kreising/DATA/data/BIRD_MICROBIOME/PASSERINES_COMPARATIVE/FIRST_RUN/DATA_RESULTS_NEW/PHYLOSEQ.dupl_samples.R")
PHYLOSEQ.dupl_samples_FIRST<-PHYLOSEQ.dupl_samples


######
#run2#
######
load("/media/kreising/DATA/data/Ptaci_RUN_24.10.2016/Run41_CBGP_JFM/03.DADA.res_NEW/PHYLOSEQ.dupl_samples.R")
PHYLOSEQ.dupl_samples_SECOND<-PHYLOSEQ.dupl_samples


##################################
#DECONTAM - FIRST RUN#############
##################################

#remove ASVs with low prevalences - cannot be statistically determined, if these represent contaminant or not
TAXSUMS<-taxa_sums(transform_sample_counts(PHYLOSEQ.dupl_samples_FIRST,function(x) ifelse(x>0,1,0)))
hist(log10(TAXSUMS))
PHYLOSEQ.dupl_samples_FIRST.trans<-transform_sample_counts(PHYLOSEQ.dupl_samples_FIRST, function(x) x/sum(x))
PHYLOSEQ.dupl_samples_FIRST.trans.sub<-prune_taxa(TAXSUMS>5,PHYLOSEQ.dupl_samples_FIRST.trans)

#identification based on PCR output (i.e. scores for gel band intensities)
summary(sample_data(PHYLOSEQ.dupl_samples_FIRST.trans.sub)$Koncentrace_na_gelu_1)
PHYLOSEQ.dupl_samples_FIRST.trans.sub<-prune_samples(!is.na(sample_data(PHYLOSEQ.dupl_samples_FIRST.trans.sub)$Koncentrace_na_gelu_1),PHYLOSEQ.dupl_samples_FIRST.trans.sub)
sample_data(PHYLOSEQ.dupl_samples_FIRST.trans.sub)$Koncentrace_na_gelu_1<-sample_data(PHYLOSEQ.dupl_samples_FIRST.trans.sub)$Koncentrace_na_gelu_1+0.01
contamdf.freq <- isContaminant(PHYLOSEQ.dupl_samples_FIRST.trans.sub, method="frequency", conc="Koncentrace_na_gelu_1")
contamdf.freq$FDR<-p.adjust(contamdf.freq$p,method = "fdr")
contamdf.freq_FIRST<-contamdf.freq

#identification based on prevalence in negative controls
sample_data(PHYLOSEQ.dupl_samples_FIRST.trans.sub)$is.neg <- regexpr("ontrol",sample_data(PHYLOSEQ.dupl_samples_FIRST.trans.sub)$species)>0
sum(sample_data(PHYLOSEQ.dupl_samples_FIRST.trans.sub)$is.neg)
contamdf.prev <- isContaminant(PHYLOSEQ.dupl_samples_FIRST.trans.sub, method="prevalence", neg="is.neg",threshold = 0.01)
contamdf.prev_FIRST<-contamdf.prev


###########################################
#DECONTAM - SECOND RUN#####################
###########################################

#remove ASVs with low prevalences - cannot be staticically determined, if these represent contaminant or not
PHYLOSEQ.dupl_samples_SECOND<-prune_taxa(taxa_sums(PHYLOSEQ.dupl_samples_SECOND)>0,PHYLOSEQ.dupl_samples_SECOND)
PHYLOSEQ.dupl_samples_SECOND<-prune_samples(sample_sums(PHYLOSEQ.dupl_samples_SECOND)>200,PHYLOSEQ.dupl_samples_SECOND)
TAXSUMS<-taxa_sums(transform_sample_counts(PHYLOSEQ.dupl_samples_SECOND,function(x) ifelse(x>0,1,0)))
hist(log10(TAXSUMS))
PHYLOSEQ.dupl_samples_SECOND.trans<-transform_sample_counts(PHYLOSEQ.dupl_samples_SECOND, function(x) x/sum(x))
PHYLOSEQ.dupl_samples_SECOND.trans.sub<-prune_taxa(TAXSUMS>5,PHYLOSEQ.dupl_samples_SECOND.trans)

#identification based on PCR output (i.e. scores for gel band intensities)
summary(sample_data(PHYLOSEQ.dupl_samples_SECOND.trans.sub)$koncentrace_gel)
PHYLOSEQ.dupl_samples_SECOND.trans.sub<-prune_samples(!is.na(sample_data(PHYLOSEQ.dupl_samples_SECOND.trans.sub)$koncentrace_gel),PHYLOSEQ.dupl_samples_SECOND.trans.sub)
sample_data(PHYLOSEQ.dupl_samples_SECOND.trans.sub)$Koncentrace_na_gelu_1<-sample_data(PHYLOSEQ.dupl_samples_SECOND.trans.sub)$koncentrace_gel+0.01
contamdf.freq <- isContaminant(PHYLOSEQ.dupl_samples_SECOND.trans.sub, method="frequency", conc="Koncentrace_na_gelu_1")
contamdf.freq_SECOND<-contamdf.freq

#identification based on prevalence in negative controls
sample_data(PHYLOSEQ.dupl_samples_SECOND.trans.sub)$is.neg <- regexpr("ontrol",sample_data(PHYLOSEQ.dupl_samples_SECOND.trans.sub)$druh)>0
sum(sample_data(PHYLOSEQ.dupl_samples_SECOND.trans.sub)$is.neg)
contamdf.prev <- isContaminant(PHYLOSEQ.dupl_samples_SECOND.trans.sub, method="prevalence", neg="is.neg",threshold = 0.01)
table(contamdf.prev$contaminant)
View(contamdf.prev)
contamdf.prev_SECOND<-contamdf.prev


#list of all contaminants for both runs
CONTAMINANTS<-c(rownames(contamdf.prev_FIRST)[contamdf.prev_FIRST$p<0.01],
                rownames(contamdf.prev_SECOND)[contamdf.prev_SECOND$p<0.01],
                rownames(contamdf.freq_FIRST)[contamdf.freq_FIRST$p<0.01],
                rownames(contamdf.freq_SECOND)[contamdf.freq_SECOND$p<0.01])

CONTAMINANTS<-unique(CONTAMINANTS)
length(CONTAMINANTS)

#load phyloseg with data from first and second run
#remove contaminating ASVs
load("/media/kreising/DATA/data/TROPICAL_TEMPERATE_Birds_2021/data/PHYLOSEQ_tropical12_tax50.final.R")
PHYLOSEQ_tropical12_tax50.final
sum(taxa_names(PHYLOSEQ_tropical12_tax50.final)%in%CONTAMINANTS)
FILTER<-!taxa_names(PHYLOSEQ_tropical12_tax50.final)%in%CONTAMINANTS
sum(FILTER==F)
PHYLOSEQ_tropical12_tax50.final.negfilt<-prune_taxa(FILTER,PHYLOSEQ_tropical12_tax50.final)


save(PHYLOSEQ_tropical12_tax80.final.negfilt,file = "/media/kreising/DATA/data/TROPICAL_TEMPERATE_Birds_2021/data/PHYLOSEQ_tropical12_tax80.final.negfilt")

######################################
#VSEARCH CLUSTERING - 97% similarity##
######################################
load("/media/kreising/DATA/data/TROPICAL_TEMPERATE_Birds_2021/data/PHYLOSEQ_tropical12_tax50.final.negfilt")

TO_CLUSTER.fasta<-refseq(PHYLOSEQ_tropical12_tax80.final.negfilt)
writeFasta(TO_CLUSTER.fasta,"/media/kreising/DATA/data/TROPICAL_TEMPERATE_Birds_2021/data/TO_CLUSTER.negfilt.fasta")
setwd("/media/kreising/DATA/data/TROPICAL_TEMPERATE_Birds_2021/data/")
system("~/software/bin/vsearch --cluster_size /media/kreising/DATA/data/TROPICAL_TEMPERATE_Birds_2021/data/TO_CLUSTER.negfilt.fasta --id 0.97 --uc UC.file.negfilt --centroids REF.centr.negfilt.fasta")


UF<-read.delim("/media/kreising/DATA/data/TROPICAL_TEMPERATE_Birds_2021/data/UC.file.negfilt",header = F,stringsAsFactors = F)
UF<-(UF[-grep("C",UF$V1),])
dim(UF)
UF$V10<-ifelse(UF$V10=="*",UF$V9,UF$V10)
UF<-UF[,9:10]
names(UF)<-c("OTU","CLUSTER")

#RDP 50% confidence
TAX50<-data.frame(tax_table(PHYLOSEQ_tropical12_tax50.final.negfilt))
TAX50<-data.frame(OTU=rownames(TAX50),TAX50)
SPECIES<-rownames(TAX50)
TAX50<-join(TAX50,UF)
RANKS<-colnames(TAX50)
TAX50<-tax_table(TAX50)
colnames(TAX50)<-RANKS
taxa_names(TAX50)<-as.character(TAX50[,colnames(TAX50)=="OTU"])
TAX50.red<-TAX50[,c(8,1)]
tax_table(PHYLOSEQ_tropical12_tax50.final.negfilt)<-TAX50.red
PHYLOSEQ_tropical12_tax50.final.negfilt_97OTU<-tax_glom(PHYLOSEQ_tropical12_tax50.final.negfilt,"CLUSTER")
PHYLOSEQ_tropical12_tax50.final.negfilt_97OTU<-merge_phyloseq(otu_table(PHYLOSEQ_tropical12_tax50.final.negfilt_97OTU),
                                                              TAX50[,-c(1,8)],
                                                              refseq(PHYLOSEQ_tropical12_tax50.final.negfilt_97OTU),
                                                              phy_tree(PHYLOSEQ_tropical12_tax50.final.negfilt_97OTU),
                                                              sample_data(PHYLOSEQ_tropical12_tax50.final.negfilt_97OTU))

save(PHYLOSEQ_tropical12_tax80.final.negfilt_97OTU,file = "/media/kreising/DATA/data/TROPICAL_TEMPERATE_Birds_2021/data/PHYLOSEQ_tropical12_tax80.final.negfilt_97OTU.R")
