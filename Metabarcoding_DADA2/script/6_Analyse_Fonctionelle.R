# !/usr/bin/env Rscript
#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 30/05/2022
#
# Script to generate taxonomic profile table
set.seed(1)
# Set directory, input and output -----------------------------------------------------------
# Detect R or Rstudio
  se <- Sys.getenv()
  if (length(se[grepl("rstudio",se,ignore.case=TRUE)]) > 0 ) { inputmode <- TRUE }
  if (length(se[grepl("rstudio",se,ignore.case=TRUE)]) == 0 ) { inputmode <- FALSE }
#
# Input argument if using Rstudio
  if (inputmode == TRUE) {
    result <- "V4-unified-correct-paired-out"
  }
#
  args = commandArgs(trailingOnly=TRUE)
  if ( inputmode == FALSE ) {
    result <- args[1]
  }
# Import package and set directory -----------------------------------------------------------
  pkg <- c("dplyr","tidyr","stringr")
  lapply(pkg, require, character.only = TRUE)
# Set directory, create file result -----------------------------------------------------------
  if (inputmode == TRUE) {
    current <- dirname(rstudioapi::getSourceEditorContext()$path)
    setwd(current)
    result <- paste(current,"../dataDADA2/result",result, sep = "/")
  }
  if (inputmode == FALSE) {
    result <- paste("../dataDADA2/result",result, sep = "/")
  }
## Set working directory
  if (dir.exists(result) == FALSE) { dir.create(result,recursive = TRUE) }
  if (dir.exists(result) == TRUE) { setwd(result) }
## Create result files
  system("mkdir Functional-Analyse")
# Input ASV Table & Definition table BTT ---------------------------------------------------------
  tableVinput <- read.csv(file = "ASV-Table/Seq-mat-rare.csv", sep = "\t")
  BTT_table <- read.csv(file = "../../../rawdata/BTTwellannoted.csv", sep = ";")
#
# Create taxonomix profiles ---------------------------------------------------------------
## Prepare taxonomic profiles with tableVinput 
  profile_tax_table <- as.data.frame(unique(tableVinput$Taxonomy)) ; colnames(profile_tax_table) <- "Taxonomy"
  profile_tax_table <- as.data.frame(profile_tax_table[order(profile_tax_table$Taxonomy),]) ; colnames(profile_tax_table) <- "Taxonomy"
  profile_tax_table <- separate(profile_tax_table, all_of(Taxonomy), c("Domain","Superphylum","Phylum","Class","Order","Family","Genus","Species"), sep =";")
  profile_tax_table[profile_tax_table=='NA'] <- NA
  profile_tax_table[,"Last"] <- NA
  profile_tax_table <- profile_tax_table %>% filter(Domain=="Eukaryota")
## Input exeption for "Last" Annotation
  exept <- c("environmental samples", "environvironmental samples", "unclassified", "unclassified","_X")
## If Last annotation doesn't exist with exeption, use prior annotation
  for (i in row.names(profile_tax_table)) {
   w <- length(profile_tax_table[i,][is.na(profile_tax_table[i,]) != TRUE])
   if (colnames(profile_tax_table[w])=="Species"){w=7}
   if ( length(grep(profile_tax_table[i,w],pattern = paste(exept,collapse="|"))) == 1 ) { 
     if ( length(grep(profile_tax_table[i,w-1],pattern = paste(exept,collapse="|"))) == 1 ) { 
       if (length(grep(profile_tax_table[i,w-2],pattern = paste(exept,collapse="|"))) == 1 ) { 
         if (length(grep(profile_tax_table[i,w-3],pattern = paste(exept,collapse="|"))) == 1 ) {
           if (length(grep(profile_tax_table[i,w-4],pattern = paste(exept,collapse="|"))) == 1 ) {profile_tax_table[i,"Last"] <- profile_tax_table[i,w-5]}
           else {profile_tax_table[i,"Last"] <- profile_tax_table[i,w-4] }}
         else {profile_tax_table[i,"Last"] <- profile_tax_table[i,w-3] }}
       else {profile_tax_table[i,"Last"] <- profile_tax_table[i,w-2] }}
     else { profile_tax_table[i,"Last"] <- profile_tax_table[i,w-1] }}
   else { profile_tax_table[i,"Last"] <- profile_tax_table[i,w] }
  }
## Join BTT table (Ramond et al., 2018)
  profile_tax_table <- left_join(x = profile_tax_table, y = BTT_table %>% select(Last,SizeMin,SizeMax,Cover,Shape,Spicule,Symmetry,Polarity,Colony,Motility,Chloroplast,Plast_Origin,Ingestion,Symbiontic,Zoosporic_Parasite,Cyst_Spore), by = "Last")
## Join size annotation for multiple Last annotation after BTT binding
  h <- 0
  while ( h < 3) { h <- h+1
  i <- 1
    while (i < nrow(profile_tax_table)-1) {
      i <- i+1
      j <- i+1
      if (profile_tax_table[i,"Last"] == profile_tax_table[j,"Last"]) { 
        profile_tax_table[i:j,"SizeMin"] <- min(c(profile_tax_table[i,"SizeMin"],profile_tax_table[j,"SizeMin"]))
        profile_tax_table[i:j,"SizeMax"] <- max(c(profile_tax_table[i,"SizeMax"],profile_tax_table[j,"SizeMax"]))}
    }
  }
## remove doublons
  profile_tax_table <- distinct(profile_tax_table)
  profile_tax_table_x <- distinct(profile_tax_table %>% select(-"Species"))
  y <- length(unique(profile_tax_table_x$Last))
  x <- 0
  for (i in unique(profile_tax_table_x$Last)) { if (length(grep(pattern = i, BTT_table$Last)) > 0) { x<- x+1}}
  z <- x*100/y
  le <- length(unique(profile_tax_table_x$Last))
## Save tables
  inftablecomp <- data.frame(z,le)
  colnames(inftablecomp)<-c("Proportion de références complétées grace à la table BTT","Références taxonomiques totales")
  write.table(inftablecomp, file = "Functional-Analyse/Taxonomic_Informations_Table.csv", append = FALSE, sep = "\t", dec = ".",row.names=FALSE)
  write.table(profile_tax_table_x, file = "Functional-Analyse/Taxonomic_Profil_Table_Final.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(profile_tax_table, file = "Functional-Analyse/Taxonomic_Profil_Table.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

