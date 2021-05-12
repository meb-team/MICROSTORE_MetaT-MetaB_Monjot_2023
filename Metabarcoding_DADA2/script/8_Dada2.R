#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# !/usr/bin/env Rscript
# 15/04/2021
#
# Script to process metabarcoding data with DADA2
# Set directory, input and output -----------------------------------------------------------
# Detect R or Rstudio
se <- Sys.getenv()
if (length(se[grepl("rstudio",se,ignore.case=TRUE)]) > 0 ) { inputmode <- TRUE }
if (length(se[grepl("rstudio",se,ignore.case=TRUE)]) == 0 ) { inputmode <- FALSE }
#
# Input argument if using Rstudio
if (inputmode == TRUE) {
  input <- "V4-unified"
  dataBase <- "pr2_version_4.13.0_18S_dada2.fasta.gz"
  result <- "V4-unified-out"
  minLen <- 200
  maxLen <- 500
  maxN <- 0
  minOverlap <- 50
  maxMismatch <- 0
  Nthread <- 6
}
#
# Set argument if using R
args = commandArgs(trailingOnly=TRUE)
if ( inputmode == FALSE ) {
  input <- args[1]
  dataBase <- args[2]
  result <- args[3]
  minLen <- as.numeric(args[4])
  maxLen <- as.numeric(args[5])
  maxN <- as.numeric(args[6])
  minOverlap <- as.numeric(args[7])
  maxMismatch <- as.numeric(args[8])
  Nthread <- as.numeric(args[9])
}
#
# Import package and palette -----------------------------------------------------------
pkg <- c("ggplot2", "readxl","dplyr","tidyr","cowplot","FactoMineR","factoextra","reshape2","varhandle","ggrepel","ggpubr","ggsci","scales","hrbrthemes","GUniFrac","svglite","treemap", "VennDiagram","stringr","dada2")
lapply(pkg, require, character.only = TRUE)
palette <- c("#AD002ACC","#EEA236CC","#00468BCC","#0099B4CC","#D62768CC","#ADB6B6CC","#42B540CC","#357EBDCC","#1B1919CC","#ED0000CC","#FF7F0ECC","#925E9FCC","#D43F3ACC","#9632B8CC","#46B8DACC","#5CB85CCC","#EFC008CC")
  #show_col(palette)

# Set directory, create file result and download dataBase -----------------------------------------------------------
if (inputmode == TRUE) {
  current <- dirname(rstudioapi::getSourceEditorContext()$path)
  setwd(current)
  input <- paste(current,"../dataDADA2/reads",input, sep = "/")
  result <- paste(current,"../dataDADA2/result",result, sep = "/")
  dataBase <- paste(current,"../dataBase",dataBase, sep = "/")
}
if (inputmode == FALSE) {
  input <- paste("../dataDADA2/reads",input, sep = "/")
  result <- paste("../../result",result, sep = "/")
  dataBase <- paste("../../../dataBase",dataBase, sep = "/")
}
if (dir.exists(input) == FALSE) { stop("input doesn't exist") }
if (dir.exists(input) == TRUE) { setwd(input) }
# Download dataBase
if (dir.exists("../../../dataBase") == FALSE) { dir.create("../../../dataBase") }
version <- as.data.frame(strsplit(basename(dataBase), "_"))[3,]
version <- paste("v",version,sep="")
link <- paste("https://github.com/pr2database/pr2database/releases/download",version,basename(dataBase),sep ="/")
command <- paste("wget",link,"-P ../../../dataBase", sep = " ")
if (dir.exists("../../../dataBase") == TRUE) { if (length(list.files(path ="../../../dataBase",pattern = basename(dataBase))) == 0) { system(command = command)}}
# Create files result
if (dir.exists("../../result") == FALSE) { dir.create("../../result") }
if (dir.exists(result) == FALSE) { dir.create(result) }
if (dir.exists(paste(result,"/Quality", sep ="")) == FALSE) { dir.create(paste(result,"/Quality", sep ="")) }
if (dir.exists(paste(result,"/Table", sep ="")) == FALSE) { dir.create(paste(result,"/Table", sep ="")) }

#
# DADA2 pipeline ----------------------------------------------------------
# Forward and reverse fastq
fnFs <- sort(list.files(pattern="R1_clean.fastq", full.names = TRUE))
fnRs <- sort(list.files(pattern="R2_clean.fastq", full.names = TRUE))
# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Inspect read quality profiles
svglite(paste(result,"Quality/Quality_Fs_aggregate.svg",sep="/"),width = 6,height = 6)
plotQualityProfile(fnFs, aggregate = TRUE)
dev.off()

svglite(paste(result,"Quality/Quality_Rs_aggregate.svg",sep="/"),width = 6,height = 6)
plotQualityProfile(fnRs, aggregate = TRUE)
dev.off()

# Filter and trim
filtFs <- file.path("filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, minLen=minLen, maxLen=maxLen,
                     maxN=maxN, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=Nthread) # On Windows set multithread=FALSE
  #head(out)

# Learn the Error Rates
errF <- learnErrors(filtFs, multithread=Nthread)
svglite(paste(result,"Quality/Quality_Fs_score.svg",sep="/"),width = 6,height = 6)
plotErrors(errF, nominalQ=TRUE)
dev.off()

errR <- learnErrors(filtRs, multithread=Nthread)
svglite(paste(result,"Quality/Quality_Rs_score.svg",sep="/"),width = 6,height = 6)
plotErrors(errR, nominalQ=TRUE)
dev.off()

# Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=Nthread)
  #dadaFs[[1]]
dadaRs <- dada(filtRs, err=errR, multithread=Nthread)
  #dadaRs[[1]]

# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, minOverlap = minOverlap, maxMismatch = maxMismatch)
  #head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
 #dim(seqtab)
  #seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=Nthread, verbose=TRUE)
 #dim(seqtab.nochim)

#sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
 #head(track)
write.table(track, file = paste(result,"Quality/track.txt",sep="/"), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

# Assign taxonomy
## Prepare taxa.table
taxa <- assignTaxonomy(seqtab.nochim, dataBase, multithread=Nthread)
taxa.table <- as.data.frame(taxa)
taxa.table$Taxonomy <- paste(taxa.table[,"Kingdom"],taxa.table[,"Phylum"],taxa.table[,"Class"],taxa.table[,"Order"],taxa.table[,"Family"],taxa.table[,"Genus"],taxa.table[,"Species"],taxa.table[,8], sep = ";")
taxa.table$Sequences <- rownames(taxa.table)
##n Prepare ASV.table
seqtab.nochim.table <- t(seqtab.nochim)
seqtab.nochim.table <- as.data.frame(seqtab.nochim.table)
seqtab.nochim.table$Sequences <- rownames(seqtab.nochim.table)
## Combine
ASV.table <- merge(seqtab.nochim.table,taxa.table %>% select(Taxonomy, Sequences),by = "Sequences")
ASV.table$ASV <- paste("ASV",1:nrow(seqtab.nochim.table),sep ="_")
rownames(ASV.table) <- ASV.table$ASV ; ASV.table <- ASV.table %>% select(-ASV)
write.table(ASV.table, file = paste(result,"Table/ASV_table.csv",sep="/"), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

# Save Rdata
save.image(file = paste(result,"Table/my_work_space.RData",sep="/"))
