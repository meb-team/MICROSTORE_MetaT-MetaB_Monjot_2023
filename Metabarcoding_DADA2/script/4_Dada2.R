# !/usr/bin/env Rscript
#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
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
  input <- "V4-unified-correct-paired"
  dataBase <- "pr2_version_4.13.0_18S_dada2.fasta.gz"
  result <- "V4-unified-correct-paired-out"
  minLen <- 200
  maxLen <- 500
  maxN <- 0
  minOverlap <- 50
  maxMismatch <- 0
  Nthread <- 6
  FWD <- "GTGYCAGCMGCCGCGGTA"
  REV <- "TTGGYRAATGCTTTCGC"
  cutadapt <- "/Users/arthur/miniconda3/envs/cutadaptenv/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
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
  FWD <- args[10]
  REV <- args[11]
  cutadapt <- args[12]
}
#
# Import package and palette -----------------------------------------------------------
pkg <- c("readxl","dplyr","tidyr","svglite","stringr","dada2","ShortRead","Biostrings")
lapply(pkg, require, character.only = TRUE)

# Set directory, create file result and download dataBase -----------------------------------------------------------
if (inputmode == TRUE) {
  current <- dirname(rstudioapi::getSourceEditorContext()$path)
  setwd(current)
  input <- paste(current,"../dataDADA2/reads",input, sep = "/")
  result <- paste(current,"../dataDADA2/rawTable",result, sep = "/")
  dataBase <- paste(current,"../dataBase",dataBase, sep = "/")
}
if (inputmode == FALSE) {
  input <- paste("../dataDADA2/reads",input, sep = "/")
  result <- paste("../../rawTable",result, sep = "/")
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
if (dir.exists("../../rawTable") == FALSE) { dir.create("../../rawTable") }
if (dir.exists(result) == FALSE) { dir.create(result) }
if (dir.exists(paste(result,"/Quality", sep ="")) == FALSE) { dir.create(paste(result,"/Quality", sep ="")) }
if (dir.exists(paste(result,"/Table", sep ="")) == FALSE) { dir.create(paste(result,"/Table", sep ="")) }

#
# DADA2 pipeline ----------------------------------------------------------
# Forward and reverse fastq
fnFs <- sort(list.files(pattern="R1_correct_paired.fastq", full.names = TRUE))
fnRs <- sort(list.files(pattern="R2_correct_paired.fastq", full.names = TRUE))

# Identify primer
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

# Filter N base
fnFs.filtN <- file.path("filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path("filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

# Search primer 
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

# Remove primer
path.cut <- file.path("cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i], "--discard-untrimmed")) # input files
}
# Check trimming
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "R1_correct_paired.fastq", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2_correct_paired.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))

# Filter and trim
filtFs <- file.path("filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, minLen=minLen, maxLen=maxLen,
                     maxN=maxN, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=Nthread) # On Windows set multithread=FALSE

# Inspect read quality profiles
svglite(paste(result,"Quality/Quality_Fs_aggregate.svg",sep="/"),width = 6,height = 6)
plotQualityProfile(filtFs, aggregate = TRUE)
dev.off()

svglite(paste(result,"Quality/Quality_Rs_aggregate.svg",sep="/"),width = 6,height = 6)
plotQualityProfile(filtRs, aggregate = TRUE)
dev.off()

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
dadaRs <- dada(filtRs, err=errR, multithread=Nthread)

# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, minOverlap = minOverlap, maxMismatch = maxMismatch)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=Nthread, verbose=TRUE)

#sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.table(track, file = paste(result,"Quality/track.txt",sep="/"), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

# Assign taxonomy
## Prepare taxa.table
taxa <- assignTaxonomy(seqtab.nochim, dataBase, multithread=Nthread)
taxa.table <- as.data.frame(taxa)
taxa.table$Taxonomy <- paste(taxa.table[,"Kingdom"],taxa.table[,"Phylum"],taxa.table[,"Class"],taxa.table[,"Order"],taxa.table[,"Family"],taxa.table[,"Genus"],taxa.table[,"Species"],taxa.table[,8], sep = ";")
taxa.table$Sequences <- rownames(taxa.table)
# Prepare ASV.table
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
