#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# !/usr/bin/env Rscript
# XX/XX/XXXX
#
# Script Composition
# Set directory, input and output -----------------------------------------------------------
# Detect R or Rstudio
se <- Sys.getenv()
if (length(se[grepl("rstudio",se,ignore.case=TRUE)]) > 0 ) { inputmode <- TRUE }
if (length(se[grepl("rstudio",se,ignore.case=TRUE)]) == 0 ) { inputmode <- FALSE }
#
# Input argument if using Rstudio
if (inputmode == TRUE) {
  input <- "Composition-V4-095-Vsearch2151-FilterOnlyOne-Unifyno-Rarefyyes-SortBHno-BH-Eukaryota"
  Mode <- "Superphylum"
  Group <- "Eukaryota"
}
#
args = commandArgs(trailingOnly=TRUE)
if ( inputmode == FALSE ) {
  input <- args[1]
}
# Import package and palette -----------------------------------------------------------
pkg <- c("ggplot2", "readxl","dplyr","tidyr","cowplot","FactoMineR","factoextra","reshape2","varhandle","ggrepel","ggpubr","ggsci","scales","hrbrthemes","GUniFrac","svglite","treemap", "VennDiagram","stringr")
lapply(pkg, require, character.only = TRUE)

#palette <-  sample(c(pal_locuszoom(alpha = 0.8)(7), pal_lancet(alpha = 0.8)(9)))
#show_col(palette)
#palette <- c("#D43F3ACC","#EEA236CC","#AD002ACC","#46B8DACC","#357EBDCC","#9632B8CC","#B8B8B8CC","#00468BCC","#ED0000CC","#42B540CC","#0099B4CC","#925E9FCC")
palette <- c("#AD002ACC","#EEA236CC","#00468BCC","#0099B4CC","#D62768CC","#ADB6B6CC","#42B540CC","#357EBDCC","#1B1919CC","#ED0000CC","#FF7F0ECC","#925E9FCC","#D43F3ACC","#9632B8CC","#46B8DACC","#5CB85CCC","#EFC008CC")
show_col(palette)

input <- paste("../result",input, sep = "/")
if (dir.exists(input) == FALSE) { stop("input doesn't exist") }
if (dir.exists(input) == TRUE) { setwd(input) }

system("mkdir Functional-Analyse")

# Input OTU Table & Definition table BTT ---------------------------------------------------------
tableVinput <- read.csv(file = "OTU-Table/Seq-mat-pool-rare.csv", sep = "\t")
BTT_table <- read.csv(file = "../../rawdata/BTTwellannoted.csv", sep = ";")

# Enter mode and Group to study -------------------------------------------

if (length(args)==2) {
  cat("Enter mode (Superphylum or Phylum or Class) : ");
  Mode <- readLines("stdin",n=1);
  cat("You entered")
  str(Mode);
  cat( "\n" )}
if (length(args)>2) {
  Mode <- args[9]
}

if (Mode == "Phylum") {
  if (length(args)==2) {
    cat("Enter Phylum (Fungi, Alveolata, etc) : ");
    Group <- readLines("stdin",n=1);
    cat("You entered")
    str(Group);
    cat( "\n" )}
  if (length(args)>2) {
    Group <- args[10]}
}
if (Mode == "Class") {
  if (length(args)==2) {
    cat("Enter Class (Bacillariophyta, Synurophyceae, etc) : ");
    Group <- readLines("stdin",n=1);
    cat("You entered")
    str(Group);
    cat( "\n" )}
  if (length(args)>2) {
    Group <- args[10]}
}
# Theme unique Dark perso -------------------------------------------------------
theme_unique_dark <- function (base_size = 12, base_family = "") {
  ret <- (theme_bw(base_size = base_size, base_family = base_family) +
            theme(text = element_text(colour = "black"),
                  title = element_text(color = "black", face = "bold"),
                  axis.ticks = element_blank(),
                  line = element_line(color = "black"),
                  rect = element_rect(fill = "white", color = "black"),
                  axis.title = element_text(color = "black", face = "bold"),
                  axis.text.y = element_blank(),
                  axis.text.x = element_blank(),
                  #axis.text.x = element_text(color = "black", size = 8, vjust = 2),
                  axis.line = element_line(color = "#969696", linetype = 1),
                  legend.background = element_rect(fill = NULL, color = NULL),
                  legend.position = "bottom",
                  legend.key = element_rect(fill = NA, color = NA, linetype = 0),
                  legend.text = element_text(size = 8),
                  legend.title = element_text(size = 8, face="bold"),
                  strip.background = element_rect(fill=NA,colour=NA,size=NA,linetype = NULL),
                  strip.text = element_text(color="black",face="bold",vjust=.5,hjust=.5),
                  panel.background = element_rect(fill = "white", color = NULL),
                  panel.border = element_blank(),
                  panel.grid = element_line(color = "#252525"),
                  panel.grid.major = element_line(color = "white"),
                  panel.grid.minor = element_line(color = "white"),
                  plot.background = element_rect(fill = "white", colour = "white", linetype = 0)))
  ret
} 
# Theme unique Dark perso 2 -------------------------------------------------------
theme_unique_darkbis <- function (base_size = 12, base_family = "") {
  ret <- (theme_bw(base_size = base_size, base_family = base_family) +
            theme(text = element_text(colour = "black"),
                  title = element_text(color = "black", face = "bold", vjust = 0.5, hjust = 0.5),
                  axis.ticks = element_blank(),
                  line = element_line(color = "black"),
                  rect = element_rect(fill = "white", color = "black"),
                  axis.title = element_text(color = "black", face = "bold"),
                  axis.text.y = element_blank(),
                  axis.text.x = element_blank(),
                  axis.line = element_line(color = "#969696", linetype = 1),
                  legend.background = element_rect(fill = NULL, color = NULL),
                  legend.position = "none",
                  legend.key = element_rect(fill = NA, color = NA, linetype = 0),
                  legend.text = element_text(size = 8),
                  legend.title = element_text(size = 8, face="bold"),
                  strip.background = element_rect(fill=NA,colour=NA,size=NA,linetype = NULL),
                  strip.text = element_text(color="black",face="bold",vjust=.5,hjust=.5),
                  panel.background = element_rect(fill = "white", color = NULL),
                  panel.border = element_blank(),
                  panel.grid = element_line(color = "#252525"),
                  panel.grid.major = element_line(color = "white"),
                  panel.grid.minor = element_line(color = "white"),
                  plot.background = element_rect(fill = "white", colour = "white", linetype = 0)))
  ret
} 




# Arrange taxonomy -------------------------------------------------------------
tax_table <- tableVinput %>% select(Taxonomy)
tax_table <- separate(tax_table, all_of(Taxonomy), c("Domain","Superphylum","Phylum","Class","Order","Family","Genus","Species","Taxo9","Taxo10","Taxo11"), sep =";")
tax_table[tax_table=='NA'] <- NA
# Create taxonomix profiles ---------------------------------------------------------------
profile_tax_table <- as.data.frame(unique(tableVinput$Taxonomy)) ; colnames(profile_tax_table) <- "Taxonomy"
profile_tax_table <- as.data.frame(profile_tax_table[order(profile_tax_table$Taxonomy),]) ; colnames(profile_tax_table) <- "Taxonomy"
profile_tax_table <- separate(profile_tax_table, all_of(Taxonomy), c("Domain","Superphylum","Phylum","Class","Order","Family","Genus","Species","Taxo9","Taxo10","Taxo11"), sep =";")
profile_tax_table[profile_tax_table=='NA'] <- NA
profile_tax_table[,"Last"] <- NA
exept <- c("environmental samples", "environvironmental samples")
for (i in row.names(profile_tax_table)) {
 w <- length(profile_tax_table[i,][is.na(profile_tax_table[i,]) != TRUE])
 if ( length(grep(profile_tax_table[i,w],all_of(exept))) == 1 ) { profile_tax_table[i,"Last"] <- profile_tax_table[i,w-1] }
 else { profile_tax_table[i,"Last"] <- profile_tax_table[i,w] }
}
profile_tax_table <- left_join(x = profile_tax_table, y = BTT_table %>% select(Last,SizeMin,SizeMax,Cover,Shape,Spicule,Symmetry,Polarity,Colony,Motility,Chloroplast,Plast_Origin,Ingestion,Symbiontic,Zoosporic_Parasite,Cyst_Spore), by = "Last")




