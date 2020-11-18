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
# Script Duplicat
# Set directory, input, output and import packages -----------------------------------------------------------
#
#output <- "Analyse-Composition-Rarefy-V9-097-199-NOfilter-Total"
#input <- "../dataPANAM/PANAM2/V9-result-097-199/OTU_distribution_tax.txt"
#region <- "V9"
#sortop <- "no"
#Taxonomy <- "NN"

#
args = commandArgs(trailingOnly=TRUE)

if (length(args)==2) {
cat("Enter ribosomal region (V4 or V9) : ");
region <- readLines("stdin",n=1);
cat("You entered")
str(region);
cat( "\n" )}
if (length(args)==4) {
region <- args[3]
}

output <- args[2]
input <- args[1]

# Import packages and create result directory ---------------------------------------------------------
pkg <- c("parallel","ggplot2", "readxl","dplyr","tidyr","cowplot","FactoMineR","factoextra","reshape2","varhandle","ggrepel","ggpubr","ggsci","scales","hrbrthemes","GUniFrac","svglite")
lapply(pkg, require, character.only = TRUE)
input <- paste("..",input, sep = "/")
output <- paste("../result",output, sep = "/")
if (dir.exists("../result") == FALSE) { dir.create("../result") }
if (dir.exists(output) == FALSE) { dir.create(output) }
if (dir.exists(output) == TRUE) { setwd(output) }

system("mkdir Figure-Sum")
system("mkdir AFC-Duplicat")
system("mkdir Composition")
system("mkdir HistOnly")
system("mkdir AFC-Distribution")
system("mkdir TableOnly")
system("mkdir Biplot")
system("mkdir Rarecurve")
system("mkdir Rarecurve/Nopool")
system("mkdir Rarecurve/Pool")
system("mkdir Diversity")
system("mkdir Diversity/Nopool")
system("mkdir Diversity/Pool")


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

# Input OTU Table ---------------------------------------------------------
tableVinput <- read.csv(file = input, sep = "\t")
##0.0005% filter
if (length(args)==2) {
  cat("Should I filter 0.0005% of total OTUs (yes or no) ? : ");
  sortop <- readLines("stdin",n=1);
  cat("You entered")
  str(sortop);
  cat( "\n" )}
if (length(args)==4) {
  sortop <- args[4]
}

if (sortop == "yes") {
  amplicon <- grep(pattern = "OSTA", colnames(tableVinput), value = TRUE)
  tableVinput$SUM <- rowSums(tableVinput %>% select(all_of(amplicon)))
  tableVinput <- tableVinput %>% filter(SUM > 0.0005*sum(tableVinput$SUM)/100)
  tableVinput <- tableVinput %>% select(-"SUM")
}
# Prepare data inf --------------------------------------------------------
infdataini <- read.table(file = "../../rawdata/data-inf.txt", sep = "\t", header = FALSE,row.names = "row", col.names = c("row","Id","Condition","Technologie","Region"))
infdataini$Rep <- rep("_1", each = nrow(infdataini))
infdataini$Variable <- paste(infdataini$Condition,infdataini$Rep, sep = "")
infdataini <- infdataini %>% select(-"Rep",-"Condition")
infdataini <- separate(infdataini, Variable, c("Condition","Date","Replicat"), sep = "_")
for (i in row.names(infdataini)) { 
  if (infdataini[i,"Replicat"] == "01") infdataini[i,"Replicat"] <- "1"
  if (infdataini[i,"Replicat"] == "02") infdataini[i,"Replicat"] <- "2"
}
infdataini$OSTA <- rep("OSTA", each = nrow(infdataini))
infdataini$Variable <- paste(infdataini$Id,infdataini$OSTA, sep = "")
infdataini <- infdataini %>% select(-"Id",-"OSTA")
infdataini <- separate(infdataini, Variable, c("Cin","Amplicon"), sep = "_")
infdataini <- infdataini %>% select(-"Cin")

## FR sample (V9 in reality and not V4)
for (i in row.names(infdataini)) { if (infdataini[i,"Amplicon"] == "FROSTA") { infdataini[i,"Region"] <- "V9"}}
for (i in row.names(infdataini)) { if (infdataini[i,"Amplicon"] == "FROSTA") { infdataini[i,"Technologie"] <- "Miseq"}}

## Prepare sample_df
infdataini <- separate(infdataini, Condition, c("empty","ADN","Cycle","Zone","Fraction"),sep = "")
infdataini <- infdataini %>% select(-"empty")
infdataini<-as.data.frame(infdataini)
col <- c("Amplicon","Technologie","Region","ADN","Cycle","Zone","Fraction","Date","Replicat")
samples_df <- infdataini[,all_of(col)]
samples_df$Condition <- paste(samples_df$ADN,samples_df$Cycle,samples_df$Zone,samples_df$Fraction, sep = "")
for (i in row.names(samples_df)) {
  if (samples_df[i,"Cycle"] == "J") { samples_df[i,"Cycle"] <- "Jour"}
  if (samples_df[i,"Cycle"] == "N") { samples_df[i,"Cycle"] <- "Nuit"}
  if (samples_df[i,"Zone"] == "O") { samples_df[i,"Zone"] <- "Oxique"}
  if (samples_df[i,"Zone"] == "A") { samples_df[i,"Zone"] <- "Anoxique"}
  if (samples_df[i,"Fraction"] == "G") { samples_df[i,"Fraction"] <- "Grande"}
  if (samples_df[i,"Fraction"] == "P") { samples_df[i,"Fraction"] <- "Petite"}}


# Select V4 or V9 ---------------------------------------------------------
samples_df<-filter(samples_df, Region == region)

# Prepare  Object -------------------------------------------------
  # Prepare seq_mat (seq_mat) ---------------------------------------------------------------------
## seq_mat
amplicon <- c(grep(pattern = "OSTA", colnames(tableVinput), value = TRUE))
seq_mat <- tableVinput %>% select(all_of(amplicon))
seq_mat[,all_of(amplicon)] <- lapply(seq_mat[,all_of(amplicon)], as.numeric)

  # Rarefy (seq_mat_rare) ---------------------------------------------------------------------
##Rarefy yes or no
if (length(args)==2) {
  cat("Should I rarefy global data (yes or no) ? : ");
  RarefyYoN <- readLines("stdin",n=1);
  cat("You entered")
  str(RarefyYoN);
  cat( "\n" )}
if (length(args)==5) {
  RarefyYoN <- args[5]
}
if (RarefyYoN == "yes") {
  seq_matt <- t(seq_mat)
  seq_mat_rare <- as.data.frame(t(Rarefy(seq_matt)$otu.tab.rff)) ## Ref : Jun Chen et al. (2012). Associating microbiome composition with environmental covariates using generalized UniFrac distances. 28(16): 2106–2113.
  rm(seq_matt)
}
if (RarefyYoN == "no") { seq_mat_rare <- seq_mat }

  # Prepare otu_mat (otu_mat and otu_mat_rare) ------------------------------------------------------------------
## otu_mat
otu_mat <- seq_mat
amplicon <- c(grep(pattern = "OSTA", colnames(seq_mat), value = TRUE))
otu_mat[,all_of(amplicon)][otu_mat[,all_of(amplicon)] != 0] <- 1
## otu_mat_rare
otu_mat_rare <- seq_mat_rare
amplicon <- c(grep(pattern = "OSTA", colnames(otu_mat_rare), value = TRUE))
otu_mat_rare[,all_of(amplicon)][otu_mat_rare[,all_of(amplicon)] != 0] <- 1

# OTU ---------------------------------------------------------------------
  # AFC Plot ----------------------------------------------------------------
  dt <- as.data.frame(otu_mat_rare)
  res.ca <- CA (dt, graph = FALSE,ncp = 2 )
  pdf("AFC-Duplicat/AFC_OTU.pdf",width = 16.00,height = 9.00)
  fviz_ca_col(res.ca, repel = TRUE, col.col = "lightgreen") + theme_unique_dark()
  dev.off()
  p <- get_ca_col(res.ca)
  coord <- p$coord
  coord <- as.data.frame(coord)
  coord$sample <- row.names(coord)
  data <- merge(x = coord, y = samples_df, by.x = "sample", by.y = "Amplicon")
  fviz_screeplot (res.ca, addlabels = TRUE, ylim = c(0, 50))
  
  Xseq <- fviz_screeplot (res.ca, addlabels = TRUE, ylim = c(0, 50))
  Xseq <- Xseq$data
  Dim1Seq <- paste("Dim 1 [",round(Xseq %>% filter(dim == 1) %>% select(eig),1),"%]", sep = "")
  Dim2Seq <- paste("Dim 2 [",round(Xseq %>% filter(dim == 2) %>% select(eig),1),"%]", sep = "")
  
  # Plot AFC ----------------------------------------------------------------
    # Overview* -------------------------------------------------------------------
  svglite("AFC-Duplicat/Analyse-OTU.svg",width = 10.00,height = 6.00)
  a <- ggplot(data, aes(y = `Dim 1`, x = `Dim 2`)) + geom_point(aes(color = Condition, shape = Date), size = 2) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    facet_grid(~Date, switch = "x") + theme_bw() +
    geom_label_repel(aes(label = sample,fill = Condition), color = 'white',size = 2.5, segment.size = 0,segment.color = "black", alpha = 0.8) +
    geom_polygon(aes(color = Condition)) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(title="Analyse Factorielle des Correspondances",x=Dim1Seq,y=Dim2Seq) +
    guides(color = FALSE)
  print(a)
  dev.off()
  
  svglite("AFC-Duplicat/Analyse-OTU-04.svg",width = 8.00,height = 6.00)
  a <- ggplot(data %>% filter(Date == "04"), aes(y = `Dim 1`, x = `Dim 2`)) + geom_point(aes(color = Condition, shape = Date), size = 3) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    geom_label_repel(aes(label = sample,fill = Condition), color = 'white',size = 4, segment.size = 0.5,segment.color = "black", alpha = 0.8) +
    #geom_polygon(aes(color = Condition)) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x=Dim1Seq,y=Dim2Seq) +
    guides(shape = FALSE, color = FALSE) + labs(color = "Condition")
  print(a)
  dev.off()
  
  svglite("AFC-Duplicat/Analyse-OTU-06.svg",width = 8.00,height = 6.00)
  a <- ggplot(data %>% filter(Date == "06"), aes(y = `Dim 1`, x = `Dim 2`)) + geom_point(aes(color = Condition, shape = Date), size = 3) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    geom_label_repel(aes(label = sample,fill = Condition), color = 'white',size = 4, segment.size = 0.5,segment.color = "black", alpha = 0.8) +
    #geom_polygon(aes(color = Condition)) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x=Dim1Seq,y=Dim2Seq) +
    guides(shape = FALSE, color = FALSE) + labs(color = "Condition")
  print(a)
  dev.off()
  
  svglite("AFC-Duplicat/Analyse-OTU-09.svg",width = 8.00,height = 6.00)
  a <- ggplot(data %>% filter(Date == "09"), aes(y = `Dim 1`, x = `Dim 2`)) + geom_point(aes(color = Condition, shape = Date), size = 3) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    geom_label_repel(aes(label = sample,fill = Condition), color = 'white',size = 4, segment.size = 0.5,segment.color = "black", alpha = 0.8) +
    #geom_polygon(aes(color = Condition)) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x=Dim1Seq,y=Dim2Seq) +
    guides(shape = FALSE, color = FALSE) + labs(color = "Condition")
  print(a)
  dev.off()
  
  svglite("AFC-Duplicat/Analyse-OTU-11.svg",width = 8.00,height = 6.00)
  a <- ggplot(data %>% filter(Date == "11"), aes(y = `Dim 1`, x = `Dim 2`)) + geom_point(aes(color = Condition, shape = Date), size = 3) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    geom_label_repel(aes(label = sample,fill = Condition), color = 'white',size = 4, segment.size = 0.5,segment.color = "black", alpha = 0.8) +
    #geom_polygon(aes(color = Condition)) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x=Dim1Seq,y=Dim2Seq) +
    guides(shape = FALSE, color = FALSE) + labs(color = "Condition")
  print(a)
  dev.off()
  

# Sequence ---------------------------------------------------------------------
  # AFC Plot ----------------------------------------------------------------
  dt <- as.data.frame(seq_mat_rare)
  res.ca <- CA (dt, graph = FALSE,ncp = 2 )
  pdf("AFC-Duplicat/AFC_Seq.pdf",width = 16.00,height = 9.00)
  fviz_ca_col(res.ca, repel = TRUE, col.col = "lightgreen") + theme_unique_dark()
  dev.off()
  p <- get_ca_col(res.ca)
  coord <- p$coord
  coord <- as.data.frame(coord)
  coord$sample <- row.names(coord)
  data <- merge(x = coord, y = samples_df, by.x = "sample", by.y = "Amplicon")
  fviz_screeplot (res.ca, addlabels = TRUE, ylim = c(0, 50))
  
  Xseq <- fviz_screeplot (res.ca, addlabels = TRUE, ylim = c(0, 50))
  Xseq <- Xseq$data
  Dim1Seq <- paste("Dim 1 [",round(Xseq %>% filter(dim == 1) %>% select(eig),1),"%]", sep = "")
  Dim2Seq <- paste("Dim 2 [",round(Xseq %>% filter(dim == 2) %>% select(eig),1),"%]", sep = "")
  
  # Plot AFC ----------------------------------------------------------------
    # Overview* -------------------------------------------------------------------
  svglite("AFC-Duplicat/Analyse-Seq.svg",width = 10.00,height = 6.00)
  a <- ggplot(data, aes(y = `Dim 1`, x = `Dim 2`)) + geom_point(aes(color = Condition, shape = Date), size = 2) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    facet_grid(~Date, switch = "x") + theme_bw() +
    geom_label_repel(aes(label = sample,fill = Condition), color = 'white',size = 2.5, segment.size = 0,segment.color = "black", alpha = 0.8) +
    geom_polygon(aes(color = Condition)) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(title="Analyse Factorielle des Correspondances",x=Dim1Seq,y=Dim2Seq) +
    guides(color = FALSE)
  print(a)
  dev.off()
  
  svglite("AFC-Duplicat/Analyse-Seq-04.svg",width = 8.00,height = 6.00)
  a <- ggplot(data %>% filter(Date == "04"), aes(y = `Dim 1`, x = `Dim 2`)) + geom_point(aes(color = Condition, shape = Date), size = 3) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    geom_label_repel(aes(label = sample,fill = Condition), color = 'white',size = 4, segment.size = 0.5,segment.color = "black", alpha = 0.8) +
    #geom_polygon(aes(color = Condition)) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x=Dim1Seq,y=Dim2Seq) +
    guides(shape = FALSE, color = FALSE) + labs(color = "Condition")
  print(a)
  dev.off()
  
  svglite("AFC-Duplicat/Analyse-Seq-06.svg",width = 8.00,height = 6.00)
  a <- ggplot(data %>% filter(Date == "06"), aes(y = `Dim 1`, x = `Dim 2`)) + geom_point(aes(color = Condition, shape = Date), size = 3) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    geom_label_repel(aes(label = sample,fill = Condition), color = 'white',size = 4, segment.size = 0.5,segment.color = "black", alpha = 0.8) +
    #geom_polygon(aes(color = Condition)) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x=Dim1Seq,y=Dim2Seq) +
    guides(shape = FALSE, color = FALSE) + labs(color = "Condition")
  print(a)
  dev.off()
  
  svglite("AFC-Duplicat/Analyse-Seq-09.svg",width = 8.00,height = 6.00)
  a <- ggplot(data %>% filter(Date == "09"), aes(y = `Dim 1`, x = `Dim 2`)) + geom_point(aes(color = Condition, shape = Date), size = 3) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    geom_label_repel(aes(label = sample,fill = Condition), color = 'white',size = 4, segment.size = 0.5,segment.color = "black", alpha = 0.8) +
    #geom_polygon(aes(color = Condition)) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x=Dim1Seq,y=Dim2Seq) +
    guides(shape = FALSE, color = FALSE) + labs(color = "Condition")
  print(a)
  dev.off()
  
  svglite("AFC-Duplicat/Analyse-Seq-11.svg",width = 8.00,height = 6.00)
  a <- ggplot(data %>% filter(Date == "11"), aes(y = `Dim 1`, x = `Dim 2`)) + geom_point(aes(color = Condition, shape = Date), size = 3) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    geom_label_repel(aes(label = sample,fill = Condition), color = 'white',size = 4, segment.size = 0.5,segment.color = "black", alpha = 0.8) +
    #geom_polygon(aes(color = Condition)) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x=Dim1Seq,y=Dim2Seq) +
    guides(shape = FALSE, color = FALSE) + labs(color = "Condition")
  print(a)
  dev.off()
  
  
