#QUALITY FILTERING AND DENOISING USING DATA

library(dada2)

####################################################
# 1 Demultiplexed fastq files were quality filtered#
####################################################

# this step was done separately for the first and the second run
setwd("/media/kreising/DATA/data/Ptaci_RUN_24.10.2016/Run41_CBGP_JFM/00.RAW.FASTQ/")
WD<-"/media/kreising/DATA/data/Ptaci_RUN_24.10.2016/Run41_CBGP_JFM/00.RAW.FASTQ/"
LIST<-list.files()
LIST<-LIST[grep("Lib_",LIST)]
LIST

i<-1
for (i in 1:length(LIST))
{
  setwd(paste(LIST[i],"/Demultiplexed/",sep=""))
  
  
  print(paste("Filtering files for:", LIST[i]))
  #list of merged fastq files:
  LIST.fastq<-list.files()
  LIST.fastq<-LIST.fastq[grep("fastq.gz",LIST.fastq)]
  LIST.fastq
  fnFs <- LIST.fastq[grepl("pair1.fastq.gz", LIST.fastq)] 
  fnRs <- LIST.fastq[grepl("pair2.fastq.gz", LIST.fastq)] 
  sample.names<-gsub("pair[12].fastq.gz","",fnFs)
  
  #POD TEMITO JMENY SE BUDOU UKLADAT FILTROVANE SOUBORY
  filtFs <- paste0(sample.names,LIST[i], "_READ1_filt.fastq.gz")
  filtRs <- paste0(sample.names,LIST[i], "_READ2_filt.fastq.gz")
  
  #FILTROVANI
  #MAXIMALNI PREDPIKLADANY POCET CHYB - 1
  for(x in seq_along(fnFs)) {
    print(fnFs[x])
    fastqPairedFilter(c(fnFs[x], fnRs[x]), c(filtFs[x], filtRs[x]),
                      maxN=0, maxEE=2, minQ=2,truncQ=2,
                      compress=TRUE, verbose=TRUE,
                      minLen = c(240,220),truncLen = c(240,220))
  }
  
  setwd(WD)
  
}

#####################################################################################
#2 filtered files were denoised using dada2 and ASV abundance matrix was constucted #
#####################################################################################

#List of quality filtered fastq files
fns <- list.files()
fastqs <- fns[grepl(".fastq.gz$", fns)]
fastqs <- sort(fastqs) 

fnFs <- fastqs[grepl("_READ1_filt.fastq.gz", fastqs)] 
fnRs <- fastqs[grepl("_READ2_filt.fastq.gz", fastqs)] 
sample.names <- gsub("_READ1_filt.fastq.gz","",fnFs)

#fastq dereplication
derepFs <- derepFastq(fnFs,n = 1e+05, verbose=F)
derepRs <- derepFastq(fnRs,n = 1e+05, verbose=F)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

#deoising
dadaFs <- dada(derepFs, selfConsist = TRUE,MAX_CONSIST=20)
dadaRs <- dada(derepRs, selfConsist = TRUE,MAX_CONSIST=20)

#merge denoised forward and reverse ASVs
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, minOverlap = 10,maxMismatch=1,justConcatenate=F)

#abundance matrix
seqtab <- makeSequenceTable(mergers)
