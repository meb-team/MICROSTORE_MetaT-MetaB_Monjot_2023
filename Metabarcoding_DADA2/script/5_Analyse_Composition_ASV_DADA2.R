# !/usr/bin/env Rscript
#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 25/05/2022
#
# Script Composition
set.seed(1)
# Set input and output -----------------------------------------------------------
## Detect R or Rstudio
  se <- Sys.getenv()
  if (length(se[grepl("rstudio",se,ignore.case=TRUE)]) > 0 ) { inputmode <- TRUE }
  if (length(se[grepl("rstudio",se,ignore.case=TRUE)]) == 0 ) { inputmode <- FALSE }
#
## Input argument if using Rstudio
  if (inputmode == TRUE) {
  input <- "V4-unified-correct-paired-out"
  result <- "V4-unified-correct-paired-out-compo"
  region <- "V4"
  sortop <- "Singleton"
  RarefyYoN <- "yes"
  }
#
## Input argument if using R
  args = commandArgs(trailingOnly=TRUE)
  if ( inputmode == FALSE ) {
    input <- args[1]
    result <- args[2]
    region <- args[3]
    sortop <- args[4]
    RarefyYoN <- args[5]
  }
#
# Import package and palette -----------------------------------------------------------
  pkg <- c("ggplot2","dplyr","tidyr","cowplot","FactoMineR","factoextra","reshape2","varhandle","ggrepel","ggpubr","ggsci","scales","hrbrthemes","GUniFrac","svglite","treemap", "VennDiagram","stringr","paletteer","elementalist","gtools","vegan","ape")
  lapply(pkg, require, character.only = TRUE)
  palette <- c("#AD002ACC","#EEA236CC","#00468BCC","#0099B4CC","#D62768CC","#ADB6B6CC","#42B540CC","#357EBDCC","#1B1919CC","#ED0000CC","#FF7F0ECC","#925E9FCC","#D43F3ACC","#9632B8CC","#46B8DACC","#5CB85CCC","#EFC008CC")
  palet_monimixo <- paletteer_d("ggthemes::Nuriel_Stone",n=2,direction=-1)
  paletspe <- rep(c("#e53a35","#f44336","#ef5350","#e57373","#ef9a9a",
                    "#7b1fa2","#8e24aa",
                    "#303f9f",
                    "#0288d1","#049be5","#64b5f6",
                    "#2e7d32","#689f38","#7cb342",
                    "#fbc02c","#fdd835","#ffeb3a","#ffee58","#fff176",
                    "#f57c00",
                    "#616161","#757575","#9e9e9e","#bdbdbd","#e0e0e0",
                    "#1876d2",
                    "#3897a7","#41acc1","#4cc6da","#80deea","#b2ebf2"),
                  4)
#
# Set directory, create file result -----------------------------------------------------------
  if (inputmode == TRUE) {
    current <- dirname(rstudioapi::getSourceEditorContext()$path)
    setwd(current)
    result <- paste(current,"../dataDADA2/result",result, sep = "/")
  }
  if (inputmode == FALSE) {
    result <- paste("../dataDADA2/result",result, sep = "/")
  }
  print(result)
  getwd()
## Set working directory
  if (dir.exists(result) == FALSE) { dir.create(result,recursive = TRUE) }
  if (dir.exists(result) == TRUE) { setwd(result) }
  getwd()
## Create result files
  dir.create("Stat-Analyse")
  dir.create("Composition")
  dir.create("Hist-Taxomy")
  dir.create("AFC-Distribution")
  dir.create("Venn")
  dir.create("ASV-Table")
  dir.create("Krona")
  dir.create("Betadisp")
  dir.create("dataTables")
#
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
                  axis.line.x = element_line(color = "#969696", linetype = 1),
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
  theme_unique_art <- function (base_size = 12, base_family = "") {
  ret <- (theme_minimal(base_size = base_size, base_family = base_family) +
            theme(text = element_text(colour = "black"),
                  title = element_text(color = "black", face = "bold", vjust = 0.5, hjust = 0.5),
                  axis.ticks = element_blank(),
                  line = element_line(color = "black"),
                  rect = element_rect(fill = "white", color = "black"),
                  axis.title = element_text(color = "black", face = "bold"),
                  #axis.text.y = element_blank(),
                  #axis.text.x = element_blank(),
                  #axis.line = element_line(color = "#969696", linetype = 1),
                  legend.background = element_blank(),
                  legend.key = element_rect(fill = NA, color = NA, linetype = 0),
                  legend.text = element_text(size = 8),
                  legend.title = element_text(size = 12, face="bold"),
                  panel.background = element_rect_round(radius = unit(0.5, "cm"),fill = "white", color = NULL,size = 1),
                  #legend.position = "none",
                  #strip.background = element_rect(fill=NA,colour=NA,size=NA,linetype = NULL),
                  strip.text = element_text(size = 12,color="black",face="bold",vjust=.5,hjust=.5),
                  #panel.background = element_rect(fill = "white", color = NULL),
                  panel.border = element_blank(),
                  panel.grid = element_line(color = "#252525"),
                  panel.grid.major = element_line(color = "grey86"),
                  panel.grid.minor = element_line(color = "grey86"),
                  plot.background = element_rect(fill = "white", colour = "white", linetype = 0)))
  ret
}
# Input ASV Table ---------------------------------------------------------
  tableVinput <- read.csv(file = paste("../../rawTable",input,"Table/ASV_table.csv",sep = "/"), sep = "\t")
## Raw_Krona
  raw_krona <- tableVinput
  amplicon <- grep(pattern = "OSTA", colnames(raw_krona), value = TRUE)
  raw_krona$SUM <- rowSums(raw_krona %>% select(all_of(amplicon)))
  raw_krona <- raw_krona %>% select(SUM,Taxonomy)
  write.table(raw_krona, file = "Krona/raw_seq_Totalx.csv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  system("cat Krona/raw_seq_Totalx.csv | tr ';' '\t' > Krona/raw_seq_Total.csv ; rm Krona/raw_seq_Totalx.csv")
## Filter
  if (sortop == "Bokulich") {
    amplicon <- grep(pattern = "OSTA", colnames(tableVinput), value = TRUE)
    tableVinput$SUM <- rowSums(tableVinput %>% select(all_of(amplicon)))
    tableVinput <- tableVinput %>% filter(SUM > 0.0005*sum(tableVinput$SUM)/100)
    tableVinput <- tableVinput %>% select(-"SUM")
  }
  if (sortop == "Singleton") {
    amplicon <- grep(pattern = "OSTA", colnames(tableVinput), value = TRUE)
    tableVinput$SUM <- rowSums(tableVinput %>% select(all_of(amplicon)))
    tableVinput <- tableVinput %>% filter(SUM > 1) 
    tableVinput <- tableVinput %>% select(-"SUM")
  }
  if (sortop == "Doubleton") {
    amplicon <- grep(pattern = "OSTA", colnames(tableVinput), value = TRUE)
    tableVinput$SUM <- rowSums(tableVinput %>% select(all_of(amplicon)))
    tableVinput <- tableVinput %>% filter(SUM > 2)
    tableVinput <- tableVinput %>% select(-"SUM")
  }
  if (sortop == "OnlyOne") {
    amplicon <- grep(pattern = "OSTA", colnames(tableVinput), value = TRUE)
    tableVinput$SUM <- rowSums(tableVinput %>% select(all_of(amplicon)))
    tableVinput$MAX <- rowMaxs(as.matrix(tableVinput %>% select(all_of(amplicon))))
    tableVinput <- tableVinput %>% filter(SUM != MAX)
    tableVinput <- tableVinput %>% select(-"SUM", -"MAX")
  }
#
# Select Eukaryota --------------------------------------------------------
  tableVinput$Domain <- 0
  for (i in row.names(tableVinput)) { tableVinput[i,"Domain"] <- str_split(tableVinput[i,"Taxonomy"],pattern = ";")[[1]][1]}
  tableVinput <- tableVinput %>% filter(Domain %in% c("Eukaryota")) %>% select(-Domain)
#
# Prepare data inf --------------------------------------------------------
  infdataini <- read.table(file = "../../../rawdata/data-inf.txt", sep = "\t", header = FALSE,row.names = "row", col.names = c("row","Id","Condition","Technologie","Region"))
  infdataini$Rep <- rep("_1", each = nrow(infdataini))
  infdataini$Variable <- paste(infdataini$Condition,infdataini$Rep, sep = "")
  infdataini <- infdataini %>% select(-"Rep",-"Condition")
  infdataini <- separate(infdataini, Variable, c("Condition","Period","Replicat"), sep = "_")
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
  col <- c("Amplicon","Technologie","Region","ADN","Cycle","Zone","Fraction","Period","Replicat")
  samples_df <- infdataini[,col]
  samples_df$Condition <- paste(samples_df$ADN,samples_df$Cycle,samples_df$Zone,samples_df$Fraction, sep = "")
  for (i in row.names(samples_df)) {
  if (samples_df[i,"Cycle"] == "J") { samples_df[i,"Cycle"] <- "Day"}
  if (samples_df[i,"Cycle"] == "N") { samples_df[i,"Cycle"] <- "Night"}
  if (samples_df[i,"Zone"] == "O") { samples_df[i,"Zone"] <- "Mixolimnion"}
  if (samples_df[i,"Zone"] == "A") { samples_df[i,"Zone"] <- "Monimolimnion"}
  if (samples_df[i,"Fraction"] == "G") { samples_df[i,"Fraction"] <- "Large"}
  if (samples_df[i,"Fraction"] == "P") { samples_df[i,"Fraction"] <- "Small"}
  }
#
# Select V4 or V9 ---------------------------------------------------------
  samples_df<-filter(samples_df, Region == region)
  write.table(samples_df, file = "ASV-Table/metadata.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Prepare  Object -------------------------------------------------
  # Prepare seq_mat ---------------------------------------------------------------------
  ## Prepare asv_mat
    amplicon <- c(grep(pattern = "OSTA", colnames(tableVinput), value = TRUE))
    seq_mat <- tableVinput %>% select(all_of(amplicon))
    seq_mat[,all_of(amplicon)] <- lapply(seq_mat[,amplicon], as.numeric)
    seq_mat$ASV_Id <- row.names(seq_mat)
  ## Remove amplicon with 0 seq
    amplicon <- grep(pattern = "OSTA", colnames(seq_mat), value = TRUE)
    seq_mat$SUM <- rowSums(seq_mat %>% select(all_of(amplicon)))
    seq_mat <- seq_mat %>% filter(SUM != 0)
    seq_mat <- seq_mat %>% select(-"SUM")
  ## Prepare correspondance table
    tblc <- samples_df %>% select("Amplicon", "Period", "Replicat", "Condition")
    tblc$Condition_Period <- paste(tblc$Condition,tblc$Period, sep = "_")
    tblc <- tblc %>% select(-"Period",-"Condition")
    tblc1 <- tblc %>% filter(tblc$Replicat == 1)
    tblc2 <- tblc %>% filter(tblc$Replicat == 2)
    tblcx <- merge(x = tblc1,y = tblc2, by = "Condition_Period", all = TRUE)
    Amplicon_without_duplicat <- c(tblcx[,"Amplicon.x"][which(is.na(tblcx[,"Amplicon.y"]) == TRUE)])
    for (i in row.names(tblcx)) { if (is.na(tblcx[i,"Amplicon.y"]) == TRUE) { tblcx[i,"Amplicon.y"] <- tblcx[i,"Amplicon.x"]}}
    for (i in row.names(tblcx)) { if (is.na(tblcx[i,"Replicat.y"]) == TRUE) { tblcx[i,"Replicat.y"] <- tblcx[i,"Replicat.x"]}}
    samples_df <- samples_df %>% filter(Replicat == 1)
#
  # Rarefy (seq_mat_rare) ---------------------------------------------------------------------
  ## Rarefy yes or no
    if (RarefyYoN == "yes") {
    seq_matt <- seq_mat %>% select(-"ASV_Id")
    seq_matt <- t(seq_matt)
    seq_mat_rare <- as.data.frame(t(Rarefy(seq_matt)$otu.tab.rff)) ## Ref : Jun Chen et al. (2012). Associating microbiome composition with environmental covariates using generalized UniFrac distances. 28(16): 2106–2113.
    rm(seq_matt)
    seq_mat_rare$ASV_Id <- rownames(seq_mat_rare)
    }
    if (RarefyYoN == "no") { seq_mat_rare <- seq_mat }
  ## Remove amplicon with 0 seq
    amplicon <- grep(pattern = "OSTA", colnames(seq_mat_rare), value = TRUE)
    seq_mat_rare$SUM <- rowSums(seq_mat_rare %>% select(all_of(amplicon)))
    seq_mat_rare <- seq_mat_rare %>% filter(SUM != 0)
    seq_mat_rare <- seq_mat_rare %>% select(-"SUM")
#
  # Prepare asv_mat (asv_mat and asv_mat_rare) ------------------------------------------------------------------
  ## asv_mat
    asv_mat <- seq_mat
    amplicon <- c(grep(pattern = "OSTA", colnames(seq_mat), value = TRUE))
    asv_mat[,all_of(amplicon)][asv_mat[,all_of(amplicon)] != 0] <- 1
  ## asv_mat_rare
    asv_mat_rare <- seq_mat_rare
    amplicon <- c(grep(pattern = "OSTA", colnames(asv_mat_rare), value = TRUE))
    asv_mat_rare[,all_of(amplicon)][asv_mat_rare[,all_of(amplicon)] != 0] <- 1
  ## Rarefaction
    a <- asv_mat_rare %>% select(-ASV_Id)
    b <- seq_mat_rare %>% select(-ASV_Id)
    a_dist <- vegdist(t(a), method="bray", diag = TRUE)
    b_dist <- vegdist(t(b), method="bray", diag = TRUE)
    capture.output(mantel(a_dist, b_dist, method="spearman", permutations=9999, strata = NULL,
           na.rm = FALSE, parallel = 10),file="Betadisp/Mantel_ASV_vs_SEQ.txt")
    
#
# Stat Rarefy -------------------------------------------------------------
## Sequence
  amplicon <- c(grep(pattern = "OSTA", colnames(seq_mat), value = TRUE))
  avRarefyS <- as.data.frame(colSums(seq_mat[,all_of(amplicon)]))
  amplicon <- c(grep(pattern = "OSTA", colnames(seq_mat_rare), value = TRUE))
  apRarefyS <- as.data.frame(colSums(seq_mat_rare[,all_of(amplicon)]))
  statRarefy <- cbind(avRarefyS,apRarefyS)
## ASV
  amplicon <- c(grep(pattern = "OSTA", colnames(asv_mat), value = TRUE))
  avRarefyO <- as.data.frame(colSums(asv_mat[,all_of(amplicon)]))
  amplicon <- c(grep(pattern = "OSTA", colnames(asv_mat_rare), value = TRUE))
  apRarefyO <- as.data.frame(colSums(asv_mat_rare[,all_of(amplicon)]))
  statRarefy <- cbind(statRarefy,avRarefyO,apRarefyO)
## Totaux
  statRarefy["Total",]<-colSums(statRarefy)
  colnames(statRarefy) <- c("avRarefy-Sequence","apRarefy-Sequence","avRarefy-ASV","apRarefy-ASV")
## Remove empty ASVs
### Avant raréfaction
  ASVav <- asv_mat
  ASVav$ASV_Id <- row.names(ASVav)
  amplicon <- c(grep(pattern = "OSTA", colnames(ASVav), value = TRUE))
  for ( i in row.names(ASVav)) { if (sum(ASVav[i,all_of(amplicon)]) == 0) { ASVav[i,"ASV_Id"] <- "Uniq"}}
  ASVav <- ASVav %>% filter(ASV_Id != "Uniq")
### Après rarefaction
  ASVap <- asv_mat_rare
  ASVap$ASV_Id <- row.names(ASVap)
  amplicon <- c(grep(pattern = "OSTA", colnames(ASVap), value = TRUE))
  for ( i in row.names(ASVap)) { if (sum(ASVap[i,all_of(amplicon)]) == 0) { ASVap[i,"ASV_Id"] <- "Uniq"}}
  ASVap <- ASVap %>% filter(ASV_Id != "Uniq")
### Final
  statRarefy["Total","avRarefy-ASV"] <- nrow(ASVav)
  statRarefy["Total","apRarefy-ASV"] <- nrow(ASVap)
  write.table(cbind(Amplicon=row.names(statRarefy),statRarefy), file = "Stat-Analyse/StatRarefy_withoutDuplicat.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  # Sum files ---------------------------------------------------------------
  ## Rarefaction stat Fig
    statRarefyASV <- statRarefy %>% select("avRarefy-ASV","apRarefy-ASV")
    colnames(statRarefyASV) <- c("No","Yes")
    statRarefyASV$Echantillon <- row.names(statRarefyASV)
    amplicon <-  c(grep(pattern = "OSTA", colnames(ASVap), value = TRUE))
    statRarefyASV <- melt(statRarefyASV[all_of(amplicon),], id = "Echantillon")
    my_comp <- list(c("Yes","No"))
  ## Plot  
    R_art <- ggplot(statRarefyASV, aes(x = variable, y = value)) + geom_violin(aes(fill = variable),trim = FALSE,width=0.3) + geom_boxplot(width=0.075) + geom_jitter(position=position_jitter(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      stat_compare_means(method="wilcox.test", label.y = max(statRarefyASV$value)-0.02*max(statRarefyASV$value),label.x = 1.3) +
      scale_fill_paletteer_d("ggthemes::Nuriel_Stone") + theme_unique_art() + theme(legend.position="none") +
      labs(x="Rarefaction",y="ASVs") #+ guides(color = FALSE)
    print(R_art)
    ggsave("Stat-Analyse/Stat-Rarefy-art.svg", device = "svg", width = 5, height = 5)
  ## Pooling stat Fig
    statRarefy <- statRarefy[!(row.names(statRarefy) %in% "Total"), ]
    statRarefy[,"Fusion"] <- "Yes"
    statRarefy[all_of(Amplicon_without_duplicat),"Fusion"] <- "No"
    statRarefy$Echantillon <- rownames(statRarefy)
    my_comp <- list(c("Yes","No"))
  ### Plot: ASV après Rarefy
  J_art <- ggplot(statRarefy, aes(x = Fusion, y = `apRarefy-ASV`)) + geom_violin(aes(fill = Fusion),trim = FALSE,width=0.3) + geom_boxplot(width=0.075) + geom_jitter(position=position_jitter(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    stat_compare_means(method="wilcox.test", label.y = max(statRarefy$`apRarefy-ASV`)-0.02*max(statRarefy$`apRarefy-ASV`),label.x = 1.3) +
    scale_fill_paletteer_d("ggthemes::Nuriel_Stone") + theme_unique_art() + theme(legend.position="none") +
    labs(x="Pooled",y="ASVs") #+ guides(color = FALSE)
  print(J_art)
  ggsave("Stat-Analyse/Analyse-SumApRare-ASV-art.svg", device = "svg", width = 5, height = 5)
  ### Plot: Séquence avant Rarefy
  G_art <- ggplot(statRarefy, aes(x = Fusion, y = `avRarefy-Sequence`)) + geom_violin(aes(fill = Fusion),trim = FALSE,width=0.3) + geom_boxplot(width=0.075) + geom_jitter(position=position_jitter(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    stat_compare_means(method="wilcox.test", label.y = max(statRarefy$`avRarefy-Sequence`)-0.02*max(statRarefy$`avRarefy-Sequence`),label.x = 1.3) +
    scale_fill_paletteer_d("ggthemes::Nuriel_Stone") + theme_unique_art() + theme(legend.position="none") +
    labs(x="Pooled",y="Sequences") #+ guides(color = FALSE)
  print(G_art)
  ggsave("Stat-Analyse/Analyse-SumAvRare-Sequences-art.svg", device = "svg", width = 5, height = 5)
  ### Plot: ASV avant Rarefy
  H_art <- ggplot(statRarefy, aes(x = Fusion, y = `avRarefy-ASV`)) + geom_violin(aes(fill = Fusion),trim = FALSE,width=0.3) + geom_boxplot(width=0.075) + geom_jitter(position=position_jitter(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    stat_compare_means(method="wilcox.test", label.y = max(statRarefy$`avRarefy-ASV`)-0.02*max(statRarefy$`avRarefy-ASV`),label.x = 1.3) +
    scale_fill_paletteer_d("ggthemes::Nuriel_Stone") + theme_unique_art() + theme(legend.position="none") +
    labs(x="Pooled",y="ASVs") #+ guides(color = FALSE)
  print(H_art)
  ggsave("Stat-Analyse/Analyse-SumAvRare-ASV-art.svg", device = "svg", width = 5, height = 5)
#
  # Taxonomy & statistics -----------------------------------------------------
  ## Create Tax stat table
    tax_table <- tableVinput %>% select(Taxonomy)
    tax_table$ASV_Id <- row.names(tax_table)
    tax_table <- separate(tax_table, Taxonomy, c("Domain","Supergroup","Division","Class","Order","Family","Genus","Species"), sep =";")
    tax_table[tax_table == ""] <- NA
    tax_table[tax_table == " "] <- NA
    tax_table[tax_table == "NA"] <- NA
  ## Continue Treatment tax_table
    for (i in c("Domain","Supergroup","Division","Class","Order","Family","Genus")) {
      j <- grep(i, colnames(tax_table)) + 1
      for (f in rownames(tax_table)) { 
        if (is.na(tax_table[f,i]) == TRUE) { tax_table[f,j] <- NA }}}
    tax_tablemix <- tax_table
    tax_tablemix$Taxonomy <- paste(tax_table[,"Domain"],tax_table[,"Supergroup"],tax_table[,"Division"],tax_table[,"Class"],tax_table[,"Order"],tax_table[,"Family"],tax_table[,"Genus"],tax_table[,"Species"], sep = ";")
    tax_tablemix <- tax_tablemix %>% select(ASV_Id, Taxonomy)
  ## Add taxonomy to seq and asv matrix
  ### Seq
    seq_mat <- merge(seq_mat,tax_tablemix,by = "ASV_Id")
    rownames(seq_mat) <- seq_mat$ASV_Id ; seq_mat <- seq_mat %>% select(-ASV_Id)
    seq_mat_rare <- merge(seq_mat_rare,tax_tablemix,by = "ASV_Id")
    rownames(seq_mat_rare) <- seq_mat_rare$ASV_Id ; seq_mat_rare <- seq_mat_rare %>% select(-ASV_Id)
  ### ASV
    asv_mat <- merge(asv_mat,tax_tablemix,by = "ASV_Id")
    rownames(asv_mat) <- asv_mat$ASV_Id ; asv_mat <- asv_mat %>% select(-ASV_Id)
    asv_mat_rare <- merge(asv_mat_rare,tax_tablemix,by = "ASV_Id")
    rownames(asv_mat_rare) <- asv_mat_rare$ASV_Id ; asv_mat_rare <- asv_mat_rare %>% select(-ASV_Id)
  ### Prepare table to statistic plot
    tax_table_rare <- tax_table %>% filter(row.names(tax_table) %in% row.names(seq_mat_rare))
    tax_table_TF <- as.data.frame(is.na(tax_table_rare))
    for (level in c("Domain","Supergroup","Division","Class","Order","Family","Genus","Species")) {
      assign(paste("count", level, sep = "_"), nrow(as.data.frame(tax_table_TF[,level][tax_table_TF[,level] == FALSE])))}
    tax_stat <- as.data.frame(c("Domain","Supergroup","Division","Class","Order","Family","Genus","Species")) ; colnames(tax_stat) <- "Level"
    tax_stat$Affiliated <- c(count_Domain,count_Supergroup,count_Division,count_Class,count_Order,count_Family,count_Genus,count_Species)
    tax_stat$`Not affiliated` <- nrow(seq_mat_rare) - tax_stat$Affiliated
    tax_stat <- melt(tax_stat, id = "Level")
  ## Plot
    tax_stat$Level <- factor(tax_stat$Level , levels=c("Domain","Supergroup","Division","Class","Order","Family","Genus","Species"))
    tax_stat$variable <- factor(tax_stat$variable , levels=c("Not affiliated","Affiliated"))
    tax_plot <- ggplot(tax_stat, mapping = aes(x= Level, y = value, fill = variable, color = variable ,linetype = variable), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
      scale_fill_manual(values = alpha(c("#1212ff","#1212ff"), c(0.25,0.5))) +
      geom_label(aes(y = value + 0.07*max(value),label = paste(value," ASVs","\n", round(value*100/nrow(seq_mat_rare),1)," %",sep =""), color = variable, alpha = variable),size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
      scale_color_manual(values = c("#0000FF00","black")) +
      scale_alpha_manual(values = c(0,0.5)) +
      scale_linetype_manual(values = c("blank","dashed")) +
      guides(fill=guide_legend("Affiliation"), color = "none", linetype = "none") +
      labs(x="Taxonomic level",y="ASVs") +
      theme(legend.title = element_text(face="bold"),
            axis.title = element_text(color = "black", face = "bold"))
    print(tax_plot)
    ggsave("Stat-Analyse/Stat-Taxo.svg", device = "svg", width = 9, height = 5)
#
  # ASV distribution statistics ---------------------------------------------
  ## Raw
    amplicon <- c(grep(pattern = "OSTA", colnames(seq_mat), value = TRUE))
    ASV_stat_raw <- as.data.frame(rowSums(seq_mat %>% select(all_of(amplicon)))) ; colnames(ASV_stat_raw) <- "Abundance"
    ASV_stat_raw$type <- "Raw"
  ## Rarefy
    amplicon <- c(grep(pattern = "OSTA", colnames(seq_mat_rare), value = TRUE))
    ASV_stat_rare <- as.data.frame(rowSums(seq_mat_rare %>% select(all_of(amplicon)))) ; colnames(ASV_stat_rare) <- "Abundance"
    ASV_stat_rare$type <- "Normalized"
  ## rbind
    ASV_stat <- rbind(ASV_stat_raw,ASV_stat_rare)
  ## x scale
    w <- 0
    i <- 1
    while ( i < max(ASV_stat$`Abundance`)) { print (i) 
      w <- c(w,i)
      i <- i*4}
  ## plot
    ASV_stat$type <- factor(ASV_stat$type , levels=c("Raw","Normalized"))
    ggplot(ASV_stat, aes(x = `Abundance`)) +facet_grid(.~type) + 
      stat_bin(binwidth = 1,color = "white", fill = "#1212ff", alpha = 0.5) +
      scale_y_continuous() +
      scale_x_continuous(trans = log2_trans(),breaks = w) +
      labs(y="ASVs") +
      theme(legend.title = element_text(face="bold"),
            axis.title = element_text(color = "black", face = "bold"),
            strip.text.x = element_text(color = "black", face = "bold", size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave("Stat-Analyse/ASV-Stat.svg", device = "svg", width = 10, height = 5)
#
# Save ASV Tables ------------------------------------------------------------
  write.table(seq_mat, file = "ASV-Table/Seq-mat.csv", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
  write.table(seq_mat_rare, file = "ASV-Table/Seq-mat-rare.csv", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
  write.table(asv_mat, file = "ASV-Table/ASV-mat.csv", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
  write.table(asv_mat_rare, file = "ASV-Table/ASV-mat-rare.csv", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
#
# Create sort Condition pattern -------------------------------------------
## Cycle
  CycleDay <- samples_df %>% filter(Cycle == "Day") %>% filter(Replicat == 1)
  CycleDay <- CycleDay$Amplicon
  CycleNight <- samples_df %>% filter(Cycle == "Night") %>% filter(Replicat == 1)
  CycleNight <- CycleNight$Amplicon
## Zone
  ZoneMixolimnion <- samples_df %>% filter(`Zone` == "Mixolimnion") %>% filter(Replicat == 1)
  ZoneMixolimnion <- ZoneMixolimnion$Amplicon
  ZoneMonimolimnion <- samples_df %>% filter(`Zone` == "Monimolimnion") %>% filter(Replicat == 1)
  ZoneMonimolimnion <- ZoneMonimolimnion$Amplicon
## Fraction
  FractionSmall <- samples_df %>% filter(`Fraction` == "Small") %>% filter(Replicat == 1)
  FractionSmall <- FractionSmall$Amplicon
  FractionLarge <- samples_df %>% filter(`Fraction` == "Large") %>% filter(Replicat == 1)
  FractionLarge <- FractionLarge$Amplicon
## Period
  Period04 <- samples_df %>% filter(Period == "04") %>% filter(Replicat == 1)
  Period04 <- Period04$Amplicon
  Period06 <- samples_df %>% filter(Period == "06") %>% filter(Replicat == 1)
  Period06 <- Period06$Amplicon
  Period09 <- samples_df %>% filter(Period == "09") %>% filter(Replicat == 1)
  Period09 <- Period09$Amplicon
  Period11 <- samples_df %>% filter(Period == "11") %>% filter(Replicat == 1)
  Period11 <- Period11$Amplicon
#
  # ASV disparity and conditions --------------------------------------------
  ## FractionXZone
    SmallXMixolimnion <- as.vector(t(samples_df %>% filter(Fraction == "Small") %>% filter(Zone == "Mixolimnion") %>% select(Amplicon)))
    SmallXMixolimnion_table <- as.data.frame(rowSums(seq_mat_rare %>% select(all_of(SmallXMixolimnion)))) ; colnames(SmallXMixolimnion_table) <- "Small & Mixolimnion"
    SmallXMonimolimnion <- as.vector(t(samples_df %>% filter(Fraction == "Small") %>% filter(Zone == "Monimolimnion") %>% select(Amplicon)))
    SmallXMonimolimnion_table <- as.data.frame(rowSums(seq_mat_rare %>% select(all_of(SmallXMonimolimnion)))) ; colnames(SmallXMonimolimnion_table) <- "Small & Monimolimnion"
    LargeXMixolimnion <- as.vector(t(samples_df %>% filter(Fraction == "Large") %>% filter(Zone == "Mixolimnion") %>% select(Amplicon)))
    LargeXMixolimnion_table <- as.data.frame(rowSums(seq_mat_rare %>% select(all_of(LargeXMixolimnion)))) ; colnames(LargeXMixolimnion_table) <- "Large & Mixolimnion"
    LargeXMonimolimnion <- as.vector(t(samples_df %>% filter(Fraction == "Large") %>% filter(Zone == "Monimolimnion") %>% select(Amplicon)))
    LargeXMonimolimnion_table <- as.data.frame(rowSums(seq_mat_rare %>% select(all_of(LargeXMonimolimnion)))) ; colnames(LargeXMonimolimnion_table) <- "Large & Monimolimnion"
    ASV_stat_rare_FractionXZone <- cbind(SmallXMixolimnion_table,SmallXMonimolimnion_table,LargeXMixolimnion_table,LargeXMonimolimnion_table)
    ASV_stat_rare_FractionXZone_melt <- melt(ASV_stat_rare_FractionXZone) ; colnames(ASV_stat_rare_FractionXZone_melt) <- c("Conditions","Abundance")
    ASV_stat_rare_FractionXZone_melt <- separate(ASV_stat_rare_FractionXZone_melt,Conditions,c("Fraction","Zone"),sep=" & ")
  ## x scale
    w <- 0
    i <- 1
    while ( i < max(ASV_stat_rare_FractionXZone_melt$Abundance)) { print (i) 
      w <- c(w,i)
      i <- i*4}
  ## plot
    ggplot(ASV_stat_rare_FractionXZone_melt, aes(x = Abundance)) + facet_grid(.~Fraction) + 
      stat_bin(binwidth = 1,color = "white", aes(fill = Zone), alpha = 1, position = "dodge2") +
      scale_y_continuous() +
      scale_x_continuous(trans = log2_trans(),breaks = w) +
      labs(y="ASVs") + geom_vline(xintercept = 182, linetype=2,color = "darkgrey", linewidth=0.5) +
      theme(legend.title = element_text(face="bold"),
            axis.title = element_text(color = "black", face = "bold"),
            strip.text.x = element_text(color = "black", face = "bold", size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1)) + theme_unique_art() + scale_fill_manual(values = palet_monimixo)
    ggsave("Stat-Analyse/ASV-Stat-FractionXZone.svg", device = "svg", width = 12, height = 5)
#
# Create dataframe Condition --------------------------------------------
  # Séquence ----------------------------------------------------------------
  ## Cycle
  ### Day
    Day_seq <- as.data.frame(rowSums(seq_mat_rare %>% select(all_of(CycleDay))))
    colnames(Day_seq) <- "TotalDay"
    Day_seq$ASV_Id <- row.names(Day_seq)
  ### Night
    Night_seq <- as.data.frame(rowSums(seq_mat_rare %>% select(all_of(CycleNight))))
    colnames(Night_seq) <- "TotalNight"
    Night_seq$ASV_Id <- row.names(Night_seq)
  ### Merge
    Cycle_seq <- merge(x = Day_seq,y = Night_seq, by = "ASV_Id")
  ## Zone
  ### Mixolimnion
    Mixolimnion_seq <- as.data.frame(rowSums(seq_mat_rare %>% select(all_of(ZoneMixolimnion))))
    colnames(Mixolimnion_seq) <- "TotalMixolimnion"
    Mixolimnion_seq$ASV_Id <- row.names(Mixolimnion_seq)
  ### Monimolimnion
    Monimolimnion_seq <- as.data.frame(rowSums(seq_mat_rare %>% select(all_of(ZoneMonimolimnion))))
    colnames(Monimolimnion_seq) <- "TotalMonimolimnion"
    Monimolimnion_seq$ASV_Id <- row.names(Monimolimnion_seq)
  ### Merge
    Zone_seq <- merge(x = Mixolimnion_seq,y = Monimolimnion_seq, by = "ASV_Id")
  ## Fraction
  ### Small
    Small_seq <- as.data.frame(rowSums(seq_mat_rare %>% select(all_of(FractionSmall))))
    colnames(Small_seq) <- "TotalSmall"
    Small_seq$ASV_Id <- row.names(Small_seq)
  ### Large
    Large_seq <- as.data.frame(rowSums(seq_mat_rare %>% select(all_of(FractionLarge))))
    colnames(Large_seq) <- "TotalLarge"
    Large_seq$ASV_Id <- row.names(Large_seq)
  ### Merge
    Fraction_seq <- merge(x = Small_seq,y = Large_seq, by = "ASV_Id")
  ## Period
  ### 04
    `04_seq` <- as.data.frame(rowSums(seq_mat_rare %>% select(all_of(Period04))))
    colnames(`04_seq`) <- "Total04"
    `04_seq`$ASV_Id <- row.names(`04_seq`)
  ### 06
    `06_seq` <- as.data.frame(rowSums(seq_mat_rare %>% select(all_of(Period06))))
    colnames(`06_seq`) <- "Total06"
    `06_seq`$ASV_Id <- row.names(`06_seq`)
  ### 09
    `09_seq` <- as.data.frame(rowSums(seq_mat_rare %>% select(all_of(Period09))))
    colnames(`09_seq`) <- "Total09"
    `09_seq`$ASV_Id <- row.names(`09_seq`)
  ### 11
    `11_seq` <- as.data.frame(rowSums(seq_mat_rare %>% select(all_of(Period11))))
    colnames(`11_seq`) <- "Total11"
    `11_seq`$ASV_Id <- row.names(`11_seq`)
  ### Merge
    Period_seq <- merge(x = `04_seq`,y = `06_seq`, by = "ASV_Id")
    Period_seq <- merge(x = Period_seq,y = `09_seq`, by = "ASV_Id")
    Period_seq <- merge(x = Period_seq,y = `11_seq`, by = "ASV_Id")
#
  # ASV ----------------------------------------------------------------
  ## Cycle
  ### Day
    Day_asv <- as.data.frame(rowSums(asv_mat_rare %>% select(all_of(CycleDay))))
    colnames(Day_asv) <- "TotalDay"
    Day_asv[,"TotalDay"][Day_asv[,"TotalDay"] > 0] <- 1
    Day_asv$ASV_Id <- row.names(Day_asv)
  ### Night
    Night_asv <- as.data.frame(rowSums(asv_mat_rare %>% select(all_of(CycleNight))))
    colnames(Night_asv) <- "TotalNight"
    Night_asv[,"TotalNight"][Night_asv[,"TotalNight"] > 0] <- 1
    Night_asv$ASV_Id <- row.names(Night_asv)
  ### Merge
    Cycle_asv <- merge(x = Day_asv,y = Night_asv, by = "ASV_Id")
  ## Zone
  ### Mixolimnion
    Mixolimnion_asv <- as.data.frame(rowSums(asv_mat_rare %>% select(all_of(ZoneMixolimnion))))
    colnames(Mixolimnion_asv) <- "TotalMixolimnion"
    Mixolimnion_asv[,"TotalMixolimnion"][Mixolimnion_asv[,"TotalMixolimnion"] > 0] <- 1
    Mixolimnion_asv$ASV_Id <- row.names(Mixolimnion_asv)
  ### Monimolimnion
    Monimolimnion_asv <- as.data.frame(rowSums(asv_mat_rare %>% select(all_of(ZoneMonimolimnion))))
    colnames(Monimolimnion_asv) <- "TotalMonimolimnion"
    Monimolimnion_asv[,"TotalMonimolimnion"][Monimolimnion_asv[,"TotalMonimolimnion"] > 0] <- 1
    Monimolimnion_asv$ASV_Id <- row.names(Monimolimnion_asv)
  ### Merge
    Zone_asv <- merge(x = Mixolimnion_asv,y = Monimolimnion_asv, by = "ASV_Id")
  ## Fraction
  ### Small
    Small_asv <- as.data.frame(rowSums(asv_mat_rare %>% select(all_of(FractionSmall))))
    colnames(Small_asv) <- "TotalSmall"
    Small_asv[,"TotalSmall"][Small_asv[,"TotalSmall"] > 0] <- 1
    Small_asv$ASV_Id <- row.names(Small_asv)
  ### Large
    Large_asv <- as.data.frame(rowSums(asv_mat_rare %>% select(all_of(FractionLarge))))
    colnames(Large_asv) <- "TotalLarge"
    Large_asv[,"TotalLarge"][Large_asv[,"TotalLarge"] > 0] <- 1
    Large_asv$ASV_Id <- row.names(Large_asv)
  ### Merge
    Fraction_asv <- merge(x = Small_asv,y = Large_asv, by = "ASV_Id")
  ## Period
  ### 04
    `04_asv` <- as.data.frame(rowSums(asv_mat_rare %>% select(all_of(Period04))))
    colnames(`04_asv`) <- "Total04"
    `04_asv`[,"Total04"][`04_asv`[,"Total04"] > 0] <- 1
    `04_asv`$ASV_Id <- row.names(`04_asv`)
  ### 06
    `06_asv` <- as.data.frame(rowSums(asv_mat_rare %>% select(all_of(Period06))))
    colnames(`06_asv`) <- "Total06"
    `06_asv`[,"Total06"][`06_asv`[,"Total06"] > 0] <- 1
    `06_asv`$ASV_Id <- row.names(`06_asv`)
  ### 09
    `09_asv` <- as.data.frame(rowSums(asv_mat_rare %>% select(all_of(Period09))))
    colnames(`09_asv`) <- "Total09"
    `09_asv`[,"Total09"][`09_asv`[,"Total09"] > 0] <- 1
    `09_asv`$ASV_Id <- row.names(`09_asv`)
  ### 11
    `11_asv` <- as.data.frame(rowSums(asv_mat_rare %>% select(all_of(Period11))))
    colnames(`11_asv`) <- "Total11"
    `11_asv`[,"Total11"][`11_asv`[,"Total11"] > 0] <- 1
    `11_asv`$ASV_Id <- row.names(`11_asv`)
  ### Merge
    Period_asv <- merge(x = `04_asv`,y = `06_asv`, by = "ASV_Id")
    Period_asv <- merge(x = Period_asv,y = `09_asv`, by = "ASV_Id")
    Period_asv <- merge(x = Period_asv,y = `11_asv`, by = "ASV_Id")
#
# AFC Plot BASE ----------------------------------------------------------------
  # Séquence ----------------------------------------------------------------
    dt <- as.data.frame(seq_mat_rare)
    dt <- dt %>% select(-"Taxonomy")
  ## CA
    res.ca <- CA(dt, graph = FALSE,ncp = 2)
    fviz_ca_row(res.ca, repel = FALSE, label = "none")
    p <- get_ca_row(res.ca)
    coord <- p$coord
    coord <- as.data.frame(coord)
    coord$ASV_Id <- row.names(coord)
  ### Cycle
    data_seq <- merge(x = coord, y = Cycle_seq, by = "ASV_Id")
  #### 100% - 0%
    data_seq$Day <- 0
    data_seq$Night <- 0
    data_seq[,"Day"][data_seq[,"TotalDay"] != 0] <-1
    data_seq[,"Night"][data_seq[,"TotalNight"] != 0] <-2
    data_seq$Cycle <- 0
    data_seq[,"Cycle"] <- data_seq[,"Day"] + data_seq[,"Night"]
    data_seq[,"Cycle"][data_seq[,"Cycle"] == 1] <- "Day"
    data_seq[,"Cycle"][data_seq[,"Cycle"] == 2] <- "Night"
    data_seq[,"Cycle"][data_seq[,"Cycle"] == 3] <- "Shared"
  #### 90% - 10%
    data_seq$Cycle90 <- "Shared"
    data_seq[,"Cycle90"][0.1*data_seq[,"TotalDay"] > data_seq[,"TotalNight"]] <- "Day"
    data_seq[,"Cycle90"][0.1*data_seq[,"TotalNight"] > data_seq[,"TotalDay"]] <- "Night"
  ### Zone
    data_seq <- merge(x = data_seq, y = Zone_seq, by = "ASV_Id")
  #### 100% - 0%
    data_seq$Mixolimnion <- 0
    data_seq$Monimolimnion <- 0
    data_seq[,"Mixolimnion"][data_seq[,"TotalMixolimnion"] != 0] <-1
    data_seq[,"Monimolimnion"][data_seq[,"TotalMonimolimnion"] != 0] <-2
    data_seq$Zone <- 0
    data_seq[,"Zone"] <- data_seq[,"Mixolimnion"] + data_seq[,"Monimolimnion"]
    data_seq[,"Zone"][data_seq[,"Zone"] == 1] <- "Mixolimnion"
    data_seq[,"Zone"][data_seq[,"Zone"] == 2] <- "Monimolimnion"
    data_seq[,"Zone"][data_seq[,"Zone"] == 3] <- "Shared"
  #### 90% - 10%
    data_seq$Zone90 <- "Shared"
    data_seq[,"Zone90"][0.1*data_seq[,"TotalMixolimnion"] > data_seq[,"TotalMonimolimnion"]] <- "Mixolimnion"
    data_seq[,"Zone90"][0.1*data_seq[,"TotalMonimolimnion"] > data_seq[,"TotalMixolimnion"]] <- "Monimolimnion"
### Fraction
    data_seq <- merge(x = data_seq, y = Fraction_seq, by = "ASV_Id")
    #### 100% - 0%
    data_seq$Small <- 0
    data_seq$Large <- 0
    data_seq[,"Small"][data_seq[,"TotalSmall"] != 0] <-1
    data_seq[,"Large"][data_seq[,"TotalLarge"] != 0] <-2
    data_seq$Fraction <- 0
    data_seq[,"Fraction"] <- data_seq[,"Small"] + data_seq[,"Large"]
    data_seq[,"Fraction"][data_seq[,"Fraction"] == 1] <- "Small"
    data_seq[,"Fraction"][data_seq[,"Fraction"] == 2] <- "Large"
    data_seq[,"Fraction"][data_seq[,"Fraction"] == 3] <- "Shared"
    #### 90% - 10%
    data_seq$Fraction90 <- "Shared"
    data_seq[,"Fraction90"][0.1*data_seq[,"TotalSmall"] > data_seq[,"TotalLarge"]] <- "Small"
    data_seq[,"Fraction90"][0.1*data_seq[,"TotalLarge"] > data_seq[,"TotalSmall"]] <- "Large"
  ### Period
    data_seq <- merge(x = data_seq, y = Period_seq, by = "ASV_Id")
    data_seq04 <- data_seq %>% filter(Total04 != 0) %>% select(-Total06,-Total09,-Total11)
    data_seq06 <- data_seq %>% filter(Total06 != 0) %>% select(-Total04,-Total09,-Total11)
    data_seq09 <- data_seq %>% filter(Total09 != 0) %>% select(-Total04,-Total06,-Total11)
    data_seq11 <- data_seq %>% filter(Total11 != 0) %>% select(-Total04,-Total06,-Total09)
  ## Dimension
    Xseq <- fviz_screeplot (res.ca, addlabels = TRUE, ylim = c(0, 50))
    Xseq <- Xseq$data
    Dim1seq <- paste("Dim 1 [",round(Xseq %>% filter(dim == 1) %>% select(eig),1),"%]", sep = "")
    Dim2seq <- paste("Dim 2 [",round(Xseq %>% filter(dim == 2) %>% select(eig),1),"%]", sep = "")
#
  # ASV ----------------------------------------------------------------
    dt <- as.data.frame(asv_mat_rare)
    dt <- dt %>% select(-"Taxonomy")
  ## CA
    res.ca <- CA(dt, graph = FALSE,ncp = 2 )
    fviz_ca_row(res.ca, repel = FALSE, label = "none")
    p <- get_ca_row(res.ca)
    coord <- p$coord
    coord <- as.data.frame(coord)
    coord$ASV_Id <- row.names(coord)
  ### Cycle
    data_asv <- merge(x = coord, y = Cycle_asv, by = "ASV_Id")
    data_asv$Day <- 0
    data_asv$Night <- 0
    data_asv[,"Day"][data_asv[,"TotalDay"] != 0] <- 1
    data_asv[,"Night"][data_asv[,"TotalNight"] != 0] <- 2
    data_asv$Cycle <- 0
    data_asv[,"Cycle"] <- data_asv[,"Day"] + data_asv[,"Night"]
    data_asv[,"Cycle"][data_asv[,"Cycle"] == 1] <- "Day"
    data_asv[,"Cycle"][data_asv[,"Cycle"] == 2] <- "Night"
    data_asv[,"Cycle"][data_asv[,"Cycle"] == 3] <- "Shared"
  ### Zone
    data_asv <- merge(x = data_asv, y = Zone_asv, by = "ASV_Id")
    data_asv$Mixolimnion <- 0
    data_asv$Monimolimnion <- 0
    data_asv[,"Mixolimnion"][data_asv[,"TotalMixolimnion"] != 0] <- 1
    data_asv[,"Monimolimnion"][data_asv[,"TotalMonimolimnion"] != 0] <- 2
    data_asv$Zone <- 0
    data_asv[,"Zone"] <- data_asv[,"Mixolimnion"] + data_asv[,"Monimolimnion"]
    data_asv[,"Zone"][data_asv[,"Zone"] == 1] <- "Mixolimnion"
    data_asv[,"Zone"][data_asv[,"Zone"] == 2] <- "Monimolimnion"
    data_asv[,"Zone"][data_asv[,"Zone"] == 3] <- "Shared"
  ### Fraction
    data_asv <- merge(x = data_asv, y = Fraction_asv, by = "ASV_Id")
    data_asv$Small <- 0
    data_asv$Large <- 0
    data_asv[,"Small"][data_asv[,"TotalSmall"] != 0] <- 1
    data_asv[,"Large"][data_asv[,"TotalLarge"] != 0] <- 2
    data_asv$Fraction <- 0
    data_asv[,"Fraction"] <- data_asv[,"Small"] + data_asv[,"Large"]
    data_asv[,"Fraction"][data_asv[,"Fraction"] == 1] <- "Small"
    data_asv[,"Fraction"][data_asv[,"Fraction"] == 2] <- "Large"
    data_asv[,"Fraction"][data_asv[,"Fraction"] == 3] <- "Shared"
  ### Period
    data_asv <- merge(x = data_asv, y = Period_asv, by = "ASV_Id")
    data_asv04 <- data_asv %>% filter(Total04 != 0) %>% select(-Total06,-Total09,-Total11)
    data_asv06 <- data_asv %>% filter(Total06 != 0) %>% select(-Total04,-Total09,-Total11)
    data_asv09 <- data_asv %>% filter(Total09 != 0) %>% select(-Total04,-Total06,-Total11)
    data_asv11 <- data_asv %>% filter(Total11 != 0) %>% select(-Total04,-Total06,-Total09)
  ## Dimension
    Xasv <- fviz_screeplot (res.ca, addlabels = TRUE, ylim = c(0, 50))
    Xasv <- Xasv$data
    Dim1asv <- paste("Dim 1 [",round(Xasv %>% filter(dim == 1) %>% select(eig),1),"%]", sep = "")
    Dim2asv <- paste("Dim 2 [",round(Xasv %>% filter(dim == 2) %>% select(eig),1),"%]", sep = "")
#
# Beta diversity script ---------------------------------------------------
  # Séquence ----------------------------------------------------------------
  dt <- as.data.frame(seq_mat_rare)
  dt <- t(dt %>% select(-"Taxonomy"))
  dt <- merge(dt, samples_df %>% select(Amplicon, Cycle, Period, Fraction, Zone), by.x = "row.names", by.y = "Amplicon")
  row.names(dt) <- dt$Row.names
  dtx <- dt %>% select(-Cycle,-Period, -Fraction, -Zone, -Row.names)
  dist <- vegdist(dtx, method ="bray")
  NMDS_tab=metaMDS(dist, parralel = 10)
## Plot
  Condition <-dt$Zone
  plot(NMDS_tab)
  ordiellipse(NMDS_tab,groups = Condition,label = T)
  orditorp(NMDS_tab,display="sites")
## ggplot
  scores(NMDS_tab,display="sites") %>%
    cbind(dt %>% select(Row.names, Cycle, Period, Fraction, Zone)) %>%
    cbind(statRarefy) %>%
    ggplot(aes(x = NMDS1, y = NMDS2)) +
    geom_point(aes(size = `apRarefy-ASV`, color = Zone)) +
    #facet_grid(~Fraction) +
    stat_ellipse(geom = "polygon", aes(group = Zone, color = Zone), alpha = 0.1, level = 0.95) +
    annotate("text", x = min(scores(NMDS_tab,display="sites")), y = max(scores(NMDS_tab,display="sites"))+1.2, label = paste0("stress: ", format(NMDS_tab$stress, digits = 4)), hjust = 0) +
    annotate("text", x = min(scores(NMDS_tab,display="sites")), y = max(scores(NMDS_tab,display="sites"))+1.2, label = paste0("stress: ", format(NMDS_tab$stress, digits = 4)), hjust = 0) +
    theme_bw() + scale_color_jco()  + geom_text(aes(label = Echantillon))
  ggsave("Betadisp/NMDS_Zone_seq.svg", device = "svg", width = 5, height = 5)
## Adonis
  capture.output(adonis2(formula = dist~dt$Zone+dt$Fraction+dt$Cycle+dt$Period,parralel = 10),file="Betadisp/Adonis_seq.txt")
## ANOSIM
  #
  Zone_anosim <- anosim(dist, dt$Zone, permutations = 9999, parallel = 10)
  capture.output(Zone_anosim,file="Betadisp/anosim_Zone_seq.txt")
  #
  Fraction_anosim <- anosim(dist, dt$Fraction, permutations = 9999, parallel = 10)
  capture.output(Fraction_anosim,file="Betadisp/anosim_Fraction_seq.txt")
  #
  Periods_anosim <- anosim(dist, dt$Period, permutations = 9999, parallel = 10)
  capture.output(Periods_anosim,file="Betadisp/anosim_Periods_seq.txt")
  #
  Cycle_anosim <- anosim(dist, dt$Cycle, permutations = 9999, parallel = 10)
  capture.output(Cycle_anosim,file="Betadisp/anosim_Cycle_seq.txt")
#
## Betadisper
  # Betadisper
  mod <- betadisper(dist,dt$Zone)
  labs <- paste("Dimension", 1:4, "(", 
                round(100*mod$eig / sum(mod$eig), 2), "%)")
  svglite("Betadisp/PCoA_Zone_seq.svg", width = 5.00, height = 5.00)
  plot(mod, cex=1, pch=15:17, cex.lab=1.25,
       xlab=labs[1], ylab=labs[2], main = "Zone",
       hull=TRUE, ellipse=FALSE, conf=0.95, lwd=2)
  dev.off()
  # Period
  mod <- betadisper(dist,dt$Period)
  labs <- paste("Dimension", 1:4, "(", 
                round(100*mod$eig / sum(mod$eig), 2), "%)")
  svglite("Betadisp/PCoA_Periods_seq.svg", width = 5.00, height = 5.00)
  plot(mod, cex=1, pch=15:17, cex.lab=1.25,
       xlab=labs[1], ylab=labs[2], main = "Periods",
       hull=TRUE, ellipse=FALSE, conf=0.95, lwd=2)
  dev.off()
  # Fraction
  mod <- betadisper(dist,dt$Fraction)
  labs <- paste("Dimension", 1:4, "(", 
                round(100*mod$eig / sum(mod$eig), 2), "%)")
  svglite("Betadisp/PCoA_Fraction_seq.svg", width = 5.00, height = 5.00)
  plot(mod, cex=1, pch=15:17, cex.lab=1.25,
       xlab=labs[1], ylab=labs[2], main = "Fraction",
       hull=TRUE, ellipse=FALSE, conf=0.95, lwd=2)
  dev.off()
  # Cycle
  mod <- betadisper(dist,dt$Cycle)
  labs <- paste("Dimension", 1:4, "(", 
                round(100*mod$eig / sum(mod$eig), 2), "%)")
  svglite("Betadisp/PCoA_Cycle_seq.svg", width = 5.00, height = 5.00)
  plot(mod, cex=1, pch=15:17, cex.lab=1.25,
       xlab=labs[1], ylab=labs[2], main = "Cycle",
       hull=TRUE, ellipse=FALSE, conf=0.95, lwd=2)
  dev.off()
#
  # ASV ----------------------------------------------------------------
  dt <- as.data.frame(asv_mat_rare)
  dt <- t(dt %>% select(-"Taxonomy"))
  dt <- merge(dt, samples_df %>% select(Amplicon, Cycle, Period, Fraction, Zone), by.x = "row.names", by.y = "Amplicon")
  row.names(dt) <- dt$Row.names
  dtx <- dt %>% select(-Cycle,-Period, -Fraction, -Zone, -Row.names)
  dist <- vegdist(dtx, method ="bray")
  NMDS_tab=metaMDS(dist, parralel = 10)
  ## Plot
  Condition <-dt$Zone
  plot(NMDS_tab)
  ordiellipse(NMDS_tab,groups = Condition,label = T)
  orditorp(NMDS_tab,display="sites")
  ## ggplot
  scores(NMDS_tab,display="sites") %>%
    cbind(dt %>% select(Row.names, Cycle, Period, Fraction, Zone)) %>%
    cbind(statRarefy) %>%
    ggplot(aes(x = NMDS1, y = NMDS2)) +
    geom_point(aes(size = `apRarefy-ASV`, color = Zone)) +
    #facet_grid(~Fraction) +
    stat_ellipse(geom = "polygon", aes(group = Zone, color = Zone), alpha = 0.1, level = 0.95) +
    annotate("text", x = min(scores(NMDS_tab,display="sites")), y = max(scores(NMDS_tab,display="sites"))+1.2, label = paste0("stress: ", format(NMDS_tab$stress, digits = 4)), hjust = 0) +
    annotate("text", x = min(scores(NMDS_tab,display="sites")), y = max(scores(NMDS_tab,display="sites"))+1.2, label = paste0("stress: ", format(NMDS_tab$stress, digits = 4)), hjust = 0) +
    theme_bw() + scale_color_jco()  + geom_text(aes(label = Echantillon))
  ggsave("Betadisp/NMDS_Zone_asv.svg", device = "svg", width = 5, height = 5)
  ## Adonis
  capture.output(adonis2(formula = dist~dt$Zone+dt$Fraction+dt$Cycle+dt$Period,parralel = 10),file="Betadisp/Adonis_asv.txt")
  ## ANOSIM
  #
  Zone_anosim_asv <- anosim(dist, dt$Zone, permutations = 9999, parallel = 10)
  capture.output(Zone_anosim_asv,file="Betadisp/anosim_Zone_asv.txt")
  #
  Fraction_anosim_asv <- anosim(dist, dt$Fraction, permutations = 9999, parallel = 10)
  capture.output(Fraction_anosim_asv,file="Betadisp/anosim_Fraction_asv.txt")
  #
  Periods_anosim_asv <- anosim(dist, dt$Period, permutations = 9999, parallel = 10)
  capture.output(Periods_anosim_asv,file="Betadisp/anosim_Periods_asv.txt")
  #
  Cycle_anosim_asv <- anosim(dist, dt$Cycle, permutations = 9999, parallel = 10)
  capture.output(Cycle_anosim_asv,file="Betadisp/anosim_Cycle_asv.txt")
  #
  ## Betadisper
  # Betadisper
  mod <- betadisper(dist,dt$Zone)
  labs <- paste("Dimension", 1:4, "(", 
                round(100*mod$eig / sum(mod$eig), 2), "%)")
  svglite("Betadisp/PCoA_Zone_asv.svg", width = 5.00, height = 5.00)
  plot(mod, cex=1, pch=15:17, cex.lab=1.25,
       xlab=labs[1], ylab=labs[2], main = "Zone",
       hull=TRUE, ellipse=FALSE, conf=0.95, lwd=2)
  dev.off()
  # Period
  mod <- betadisper(dist,dt$Period)
  labs <- paste("Dimension", 1:4, "(", 
                round(100*mod$eig / sum(mod$eig), 2), "%)")
  svglite("Betadisp/PCoA_Periods_asv.svg", width = 5.00, height = 5.00)
  plot(mod, cex=1, pch=15:17, cex.lab=1.25,
       xlab=labs[1], ylab=labs[2], main = "Periods",
       hull=TRUE, ellipse=FALSE, conf=0.95, lwd=2)
  dev.off()
  # Fraction
  mod <- betadisper(dist,dt$Fraction)
  labs <- paste("Dimension", 1:4, "(", 
                round(100*mod$eig / sum(mod$eig), 2), "%)")
  svglite("Betadisp/PCoA_Fraction_asv.svg", width = 5.00, height = 5.00)
  plot(mod, cex=1, pch=15:17, cex.lab=1.25,
       xlab=labs[1], ylab=labs[2], main = "Fraction",
       hull=TRUE, ellipse=FALSE, conf=0.95, lwd=2)
  dev.off()
  # Cycle
  mod <- betadisper(dist,dt$Cycle)
  labs <- paste("Dimension", 1:4, "(", 
                round(100*mod$eig / sum(mod$eig), 2), "%)")
  svglite("Betadisp/PCoA_Cycle_asv.svg", width = 5.00, height = 5.00)
  plot(mod, cex=1, pch=15:17, cex.lab=1.25,
       xlab=labs[1], ylab=labs[2], main = "Cycle",
       hull=TRUE, ellipse=FALSE, conf=0.95, lwd=2)
  dev.off()
  #
  
# AFC Taxonomy ------------------------------------------------------------
  # Séquence ----------------------------------------------------------------
    data_seq_tax <- data_seq
    colnames(tax_tablemix) <- c("ASV_Id","Taxonomy")
    data_seq_tax <- merge(data_seq_tax,tax_tablemix, by = "ASV_Id")
    data_seq_tax <- separate(data_seq_tax, Taxonomy,c("Domain","Supergroup","Division","Class","Order","Family","Genus","Species"), sep =";")
    data_seq_tax[data_seq_tax == ""] <- NA
    data_seq_tax[data_seq_tax == " "] <- NA
    data_seq_tax[data_seq_tax == "NA"] <- NA
    data_seq_tax[is.na(data_seq_tax) == TRUE] <- "Not Affiliated"
    sortDomain <- c("Eukaryota","Bacteria","Not Affiliated")
    data_seq_tax <- data_seq_tax %>% filter(Domain %in% sortDomain)
    for (i in row.names(data_seq_tax)) {
      if (data_seq_tax[i,"Domain"] == "Bacteria") { 
        data_seq_tax[i,"Supergroup"] <- "Bacteria_X"
        data_seq_tax[i,"Division"] <- "Bacteria_XX"
        data_seq_tax[i,"Class"] <- "Bacteria_XXX"
        data_seq_tax[i,"Order"] <- "Bacteria_XXXX"
        data_seq_tax[i,"Family"] <- "Bacteria_XXXXX"
        data_seq_tax[i,"Genus"] <- "Bacteria_XXXXXX"
        data_seq_tax[i,"Species"] <- "Bacteria_XXXXXXX"
      }
    }
#
    # Export data_seq_tax for Krona -------------------------------------------
    ## Day
      write.table(data_seq_tax %>% select(TotalDay, Domain, Supergroup, Division, Class, Order, Family, Genus, Species), file = "Krona/data_seq_Day.csv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    ## Night
      write.table(data_seq_tax %>% select(TotalNight, Domain, Supergroup, Division, Class, Order, Family, Genus, Species), file = "Krona/data_seq_Night.csv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    ## Mixolimnion
      write.table(data_seq_tax %>% select(TotalMixolimnion, Domain, Supergroup, Division, Class, Order, Family, Genus, Species), file = "Krona/data_seq_Mixolimnion.csv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    ## Monimolimnion
      write.table(data_seq_tax %>% select(TotalMonimolimnion, Domain, Supergroup, Division, Class, Order, Family, Genus, Species), file = "Krona/data_seq_Monimolimnion.csv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    ## Large
      write.table(data_seq_tax %>% select(TotalLarge, Domain, Supergroup, Division, Class, Order, Family, Genus, Species), file = "Krona/data_seq_Large.csv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    ## Small
      write.table(data_seq_tax %>% select(TotalSmall, Domain, Supergroup, Division, Class, Order, Family, Genus, Species), file = "Krona/data_seq_Small.csv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    ## Avril
      write.table(data_seq_tax %>% select(Total04, Domain, Supergroup, Division, Class, Order, Family, Genus, Species), file = "Krona/data_seq_Avril.csv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    ## Juin
      write.table(data_seq_tax %>% select(Total06, Domain, Supergroup, Division, Class, Order, Family, Genus, Species), file = "Krona/data_seq_Juin.csv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    ## Septembre
      write.table(data_seq_tax %>% select(Total09, Domain, Supergroup, Division, Class, Order, Family, Genus, Species), file = "Krona/data_seq_Septembre.csv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    ## Novembre
      write.table(data_seq_tax %>% select(Total11, Domain, Supergroup, Division, Class, Order, Family, Genus, Species), file = "Krona/data_seq_Novembre.csv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    ## Total
      data_seq_tax_Sum <- data_seq_tax %>% select(TotalDay, TotalNight, Domain, Supergroup, Division, Class, Order, Family, Genus, Species)
      data_seq_tax_Sum$Total <- data_seq_tax_Sum$TotalDay + data_seq_tax_Sum$TotalNight
      data_seq_tax_Sum <- data_seq_tax_Sum %>% select(Total, Domain, Supergroup, Division, Class, Order, Family, Genus, Species)
      write.table(data_seq_tax_Sum %>% select(Total, Domain, Supergroup, Division, Class, Order, Family, Genus, Species), file = "Krona/data_seq_Total.csv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    ## Launch Krona
      system("for t in $(ls Krona/ | grep '.csv$')
      do
      perl ../../../script/KronaTools-2.8/scripts/ImportText.pl Krona/$t -o Krona/$t.html
      done")
#
    # Sort interest RangInterest ------------------------------------------------------
      data_seq_tax$RangInterest <- data_seq_tax$Supergroup 
#
    # 100% - 0% -------------------------------------------------------------
      # Special Zone ------------------------------------------------------------
        data_seq_tax_specialZ <- data_seq_tax
        data_seq_tax_specialZ$Period<-"Shared"
        for (i in row.names(data_seq_tax_specialZ)) { 
          if (data_seq_tax_specialZ[i,"Total04"]>0 && data_seq_tax_specialZ[i,"Total06"]==0 && data_seq_tax_specialZ[i,"Total09"]==0 && data_seq_tax_specialZ[i,"Total11"]==0) { data_seq_tax_specialZ[i,"Period"] <- "Avril"}
          if (data_seq_tax_specialZ[i,"Total06"]>0 && data_seq_tax_specialZ[i,"Total04"]==0 && data_seq_tax_specialZ[i,"Total09"]==0 && data_seq_tax_specialZ[i,"Total11"]==0) { data_seq_tax_specialZ[i,"Period"] <- "Juin"}
          if (data_seq_tax_specialZ[i,"Total09"]>0 && data_seq_tax_specialZ[i,"Total06"]==0 && data_seq_tax_specialZ[i,"Total04"]==0 && data_seq_tax_specialZ[i,"Total11"]==0) { data_seq_tax_specialZ[i,"Period"] <- "Septembre"}
          if (data_seq_tax_specialZ[i,"Total11"]>0 && data_seq_tax_specialZ[i,"Total06"]==0 && data_seq_tax_specialZ[i,"Total09"]==0 && data_seq_tax_specialZ[i,"Total04"]==0) { data_seq_tax_specialZ[i,"Period"] <- "Novembre"}
        }
        data_seq_tax_specialZ$Condition<-NULL
        for (i in row.names(data_seq_tax_specialZ)) { if (is.na(data_seq_tax_specialZ[i,"Period"])==FALSE ) {data_seq_tax_specialZ[i,"Condition"] <- paste(data_seq_tax_specialZ[i,"Zone90"],data_seq_tax_specialZ[i,"Period"],sep=" & ") }}
        data_seq_tax_specialZ <- data_seq_tax_specialZ %>% filter(is.na(Condition)==FALSE)
      ##Figure Zone
        multiZone100 <- ggplot(data_seq_tax_specialZ, aes(y = `Dim 2`, x = `Dim 1`, color = Period)) + geom_point(size = 2) +
          geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
          geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
          stat_ellipse(aes(alpha = Period,linetype =Period),geom = "polygon",type = "norm") +
          scale_linetype_manual(values = c(2,2,2,2,0)) +
          theme(axis.title = element_text(face="bold", size=12), 
                axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
                title = element_text(face="bold", size=14),
                legend.title = element_text(face="bold"),
                legend.position = "right",
                legend.text = element_text(size=10)) +
          scale_alpha_manual(values = c(0.1,0.1,0.1,0.1,0)) +
          scale_color_manual(values = c("#F8766D","#80B69B","#7F76A7","#00A5FF","lightgrey")) +
          labs(x=Dim1seq,y=Dim2seq,color = "Period", alpha = "Period", linetype = "Period")
        print(multiZone100)
        ggsave("Betadisp/Period100.svg", device = "svg", width = 5, height = 5)
#
      # All_2018 ----------------------------------------------------------
        Var_Condition <- c("Zone","Cycle","Fraction")
        myplots100_All_2018 = lapply(Var_Condition,function(Var_Condition) {
          if (Var_Condition == "Zone") { order <- c("Monimolimnion","Mixolimnion","Shared")}
          if (Var_Condition == "Cycle") { order <- c("Day","Night","Shared")}
          if (Var_Condition == "Fraction") { order <- c("Large","Small","Shared")}
          ggplot(data_seq_tax, aes(y = `Dim 2`, x = `Dim 1`, color = factor(get(Var_Condition),level=order))) + geom_point(size = 2) +
            geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
            geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
            stat_ellipse(aes(linetype = factor(get(Var_Condition),level=order), alpha = factor(get(Var_Condition),level=order)),geom = "polygon",type = "norm") +
            scale_linetype_manual(values = c(2,2,0)) +
            theme(axis.title = element_text(face="bold", size=12), 
                  axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
                  title = element_text(face="bold", size=14),
                  legend.title = element_text(face="bold"),
                  legend.position = "right",
                  legend.text = element_text(size=10)) +
            scale_alpha_manual(values = c(0.1,0.1,0)) +
            scale_color_manual(values = c("#F8766D","#00A5FF","lightgrey")) +
            labs(x=Dim1seq,y=Dim2seq,color = Var_Condition, alpha = Var_Condition, linetype = Var_Condition)
          }
        )
        svglite("AFC-Distribution/AFC-Sequence-100.svg",width = 12.00,height = 8.00)
        print(plot_grid(plot_grid(plotlist = myplots100_All_2018), ncol = 1, nrow = 1))
        dev.off()
#
      # Each Periods ------------------------------------------------------------
        for (Var_Periods in c("04","06","09","11")){ print(Var_Periods)
          assign(paste0("data_seq",Var_Periods,"_tax"),merge(x=get(paste0("data_seq",Var_Periods)),y=tax_tablemix,by = "ASV_Id"))
          Var_Condition <- c("Zone","Cycle","Fraction")
          myplots_Each_Periods = lapply(Var_Condition,function(Var_Condition) {
            if (Var_Condition == "Zone") { order <- c("Monimolimnion","Mixolimnion","Shared")}
            if (Var_Condition == "Cycle") { order <- c("Day","Night","Shared")}
            if (Var_Condition == "Fraction") { order <- c("Large","Small","Shared")}
            ggplot(get(paste0("data_seq",Var_Periods,"_tax")), aes(y = `Dim 2`, x = `Dim 1`, color = factor(get(Var_Condition),level=order))) + geom_point(size = 2) +
              geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
              geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
              stat_ellipse(aes(linetype = factor(get(Var_Condition),level=order), alpha = factor(get(Var_Condition),level=order)),geom = "polygon",type = "norm") +
              scale_linetype_manual(values = c(2,2,0)) +
              theme(axis.title = element_text(face="bold", size=12), 
                    axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
                    title = element_text(face="bold", size=14),
                    legend.title = element_text(face="bold"),
                    legend.position = "right",
                    legend.text = element_text(size=10)) +
              scale_alpha_manual(values = c(0.1,0.1,0)) +
              scale_color_manual(values = c("#F8766D","#00A5FF","lightgrey")) +
              labs(x=Dim1seq,y=Dim2seq,color = Var_Condition, alpha = Var_Condition, linetype = Var_Condition)
            }
          )
          title <- paste0("AFC-Distribution/AFC-Sequence-",Var_Periods,"-100.svg")
          svglite(title,width = 12.00,height = 8.00)
          print(plot_grid(plot_grid(plotlist = myplots_Each_Periods), ncol = 1, nrow = 1, rel_widths = c(3,3),rel_heights = c(3,3)))
          dev.off()
        }
#
    # 90% - 10% -------------------------------------------------------------
      # Special Zone ------------------------------------------------------------
        data_seq_tax_specialZ <- data_seq_tax
        data_seq_tax_specialZ$Period<-"Shared"
        for (i in row.names(data_seq_tax_specialZ)) { 
          if (data_seq_tax_specialZ[i,"Total04"]*0.10 > data_seq_tax_specialZ[i,"Total06"] + data_seq_tax_specialZ[i,"Total09"] + data_seq_tax_specialZ[i,"Total11"]) { data_seq_tax_specialZ[i,"Period"] <- "April"}
          if (data_seq_tax_specialZ[i,"Total06"]*0.10 > data_seq_tax_specialZ[i,"Total04"] + data_seq_tax_specialZ[i,"Total09"] + data_seq_tax_specialZ[i,"Total11"]) { data_seq_tax_specialZ[i,"Period"] <- "June"}
          if (data_seq_tax_specialZ[i,"Total09"]*0.10 > data_seq_tax_specialZ[i,"Total06"] + data_seq_tax_specialZ[i,"Total04"] + data_seq_tax_specialZ[i,"Total11"]) { data_seq_tax_specialZ[i,"Period"] <- "September"}
          if (data_seq_tax_specialZ[i,"Total11"]*0.10 > data_seq_tax_specialZ[i,"Total06"] + data_seq_tax_specialZ[i,"Total09"] + data_seq_tax_specialZ[i,"Total04"]) { data_seq_tax_specialZ[i,"Period"] <- "November"}
        }
        data_seq_tax_specialZ$Condition<-NULL
        for (i in row.names(data_seq_tax_specialZ)) { if (is.na(data_seq_tax_specialZ[i,"Period"])==FALSE ) {data_seq_tax_specialZ[i,"Condition"] <- paste(data_seq_tax_specialZ[i,"Zone90"],data_seq_tax_specialZ[i,"Period"],sep=" & ") }}
        data_seq_tax_specialZ <- data_seq_tax_specialZ %>% filter(is.na(Condition)==FALSE)
      ##Figure Zone
        multiZone90 <- ggplot(data_seq_tax_specialZ, aes(y = `Dim 2`, x = `Dim 1`, color = factor(Period,level=c("April","June","September","November","Shared")))) + geom_point(size = 2) +
          geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
          geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
          stat_ellipse(aes(alpha = factor(Period,level=c("April","June","September","November","Shared")),linetype =factor(Period,level=c("April","June","September","November","Shared"))),geom = "polygon",type = "norm") +
          scale_linetype_manual(values = c(2,2,2,2,0)) +
          theme(axis.title = element_text(face="bold", size=12), 
                axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
                title = element_text(face="bold", size=14),
                legend.title = element_text(face="bold"),
                legend.position = "right",
                legend.text = element_text(size=10)) +
          scale_alpha_manual(values = c(0.1,0.1,0.1,0.1,0)) +
          scale_color_manual(values = c("#F8766D","#80B69B","#7F76A7","#00A5FF","lightgrey")) +
          labs(x=Dim1seq,y=Dim2seq,color = "Period", alpha = "Period", linetype = "Period")
        print(multiZone90)
        ggsave("Betadisp/Period90.svg", device = "svg", width = 5, height = 5)
#
      # All_2018 ----------------------------------------------------------
        Var_Condition <- c("Zone90","Cycle90","Fraction90")
        myplots90_All_2018 = lapply(Var_Condition,function(Var_Condition) {
          if (Var_Condition == "Zone90") { order <- c("Monimolimnion","Mixolimnion","Shared")}
          if (Var_Condition == "Cycle90") { order <- c("Day","Night","Shared")}
          if (Var_Condition == "Fraction90") { order <- c("Large","Small","Shared")}
          ggplot(data_seq_tax, aes(y = `Dim 2`, x = `Dim 1`, color = factor(get(Var_Condition),level=order))) + geom_point(size = 2) +
            geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
            geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
            stat_ellipse(aes(linetype = factor(get(Var_Condition),level=order), alpha = factor(get(Var_Condition),level=order)),geom = "polygon",type = "norm") +
            scale_linetype_manual(values = c(2,2,0)) +
            theme(axis.title = element_text(face="bold", size=12), 
                  axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
                  title = element_text(face="bold", size=14),
                  legend.title = element_text(face="bold"),
                  legend.position = "right",
                  legend.text = element_text(size=10)) +
            scale_alpha_manual(values = c(0.1,0.1,0)) +
            scale_color_manual(values = c("#F8766D","#00A5FF","lightgrey")) +
            labs(x=Dim1seq,y=Dim2seq,color = Var_Condition, alpha = Var_Condition, linetype = Var_Condition)
        }
        )
        svglite("AFC-Distribution/AFC-Sequence-90.svg",width = 12.00,height = 8.00)
        print(plot_grid(plot_grid(plotlist = myplots90_All_2018), ncol = 1, nrow = 1))
        dev.off()
#
      # Final_Fig ---------------------------------------------------------------
      ## Cycle
        aCycle90 <- myplots90_All_2018[[2]] +
          scale_y_continuous(trans = "reverse", limits = c(3.1,-3.1),position = 'right') +
          theme(legend.position = c(0.87, 0.2),
                legend.background = element_rect(fill = "white", color = "black")) +
          xlim(-3.3,3.2) +
          scale_color_manual(values = c("#ad494aFF","#7375b5FF","lightgrey")) +
          labs(x=Dim1seq,y=Dim2seq,color = "Cycles", alpha = "Cycles", linetype = "Cycles") +
          annotate(geom="text", x=-3.3, y=2.7, col="black", fontface = "bold", size =3.5, hjust = 0,
                   label=paste("ANOSIM test: ",round(Cycle_anosim$statistic,4),"\np-value: ",round(Cycle_anosim$signif,4)))
        print(aCycle90)
      ## Fraction
        cFraction90 <- myplots90_All_2018[[3]] +
          scale_x_continuous(limits = c(-3.3,3.2), position = 'top') + scale_y_continuous(trans = "reverse", limits = c(3.1,-3.1),position = 'right') +
          theme(legend.position = c(0.87, 0.2),
                legend.background = element_rect(fill = "white", color = "black")) +
          scale_color_manual(values = c("#ad494aFF","#7375b5FF","lightgrey")) +
          labs(x=Dim1seq,y=Dim2seq,color = "Fractions", alpha = "Fractions", linetype = "Fractions") +
          annotate(geom="text", x=-3.3, y=2.7, col="black", fontface = "bold", size =3.5, hjust = 0,
                   label=paste("ANOSIM test: ",round(Fraction_anosim$statistic,4),"\np-value: ",round(Fraction_anosim$signif,4)))
        print(cFraction90)
      ## Zone
        bZone90 <- myplots90_All_2018[[1]] +
          scale_x_continuous(limits = c(-3.3,3.2), position = 'top') +
          theme(legend.position = c(0.85, 0.2),
                legend.background = element_rect(fill = "white", color = "black")) +
          ylim(3.1,-3.1) +
          scale_color_manual(values = c("#ad494aFF","#7375b5FF","lightgrey")) +
          labs(x=Dim1seq,y=Dim2seq,color = "Zones", alpha = "Zones", linetype = "Zones") +
          annotate(geom="text", x=-3.3, y=2.7, col="black", fontface = "bold", size =3.5, hjust = 0,
                   label=paste("ANOSIM test: ",round(Zone_anosim$statistic,4),"\np-value: ",round(Zone_anosim$signif,4)))
        print(bZone90)
      ## Periods
        multiZone90x <- multiZone90 +
          theme(legend.position = c(0.87, 0.27),
                legend.background = element_rect(fill = "white", color = "black")) +
          xlim(-3.3,3.2) +
          ylim(3.1,-3.1) +
          scale_color_manual(values = c("#7375b5FF","#8ca252FF","#f6ae34FF","#ad494aFF","lightgrey")) +
          labs(x=Dim1seq,y=Dim2seq,color = "Periods", alpha = "Periods", linetype = "Periods") +
          annotate(geom="text", x=-3.3, y=2.7, col="black", fontface = "bold", size =3.5, hjust = 0,
                   label=paste("ANOSIM test: ",round(Periods_anosim$statistic,4),"\np-value: ",round(Periods_anosim$signif,4)))
        print(multiZone90x)
        ## Coplot
        svglite("AFC-Distribution/AFC-Sequence-90-article.svg",width = 14.00,height = 8.00)
        All_plot90 <- plot_grid( bZone90,cFraction90,multiZone90x,aCycle90, ncol = 2, nrow = 2, rel_widths = c(4,4), rel_heights = c(3),labels=c("A","B","C","D"),align="v")
        print(All_plot90)
        dev.off()
#
      # Each Periods ------------------------------------------------------------
        for (Var_Periods in c("04","06","09","11")){ print(Var_Periods)
          assign(paste0("data_seq",Var_Periods,"_tax"),merge(x=get(paste0("data_seq",Var_Periods)),y=tax_tablemix,by = "ASV_Id"))
          Var_Condition <- c("Zone90","Cycle90","Fraction90")
          myplots_Each_Periods = lapply(Var_Condition,function(Var_Condition) {
            if (Var_Condition == "Zone90") { order <- c("Monimolimnion","Mixolimnion","Shared")}
            if (Var_Condition == "Cycle90") { order <- c("Day","Night","Shared")}
            if (Var_Condition == "Fraction90") { order <- c("Large","Small","Shared")}
            ggplot(get(paste0("data_seq",Var_Periods,"_tax")), aes(y = `Dim 2`, x = `Dim 1`, color = factor(get(Var_Condition),level=order))) + geom_point(size = 2) +
              geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
              geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
              stat_ellipse(aes(linetype = factor(get(Var_Condition),level=order), alpha = factor(get(Var_Condition),level=order)),geom = "polygon",type = "norm") +
              scale_linetype_manual(values = c(2,2,0)) +
              theme(axis.title = element_text(face="bold", size=12), 
                    axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
                    title = element_text(face="bold", size=14),
                    legend.title = element_text(face="bold"),
                    legend.position = "right",
                    legend.text = element_text(size=10)) +
              scale_alpha_manual(values = c(0.1,0.1,0)) +
              scale_color_manual(values = c("#F8766D","#00A5FF","lightgrey")) +
              labs(x=Dim1seq,y=Dim2seq,color = Var_Condition, alpha = Var_Condition, linetype = Var_Condition)
          }
          )
          title <- paste0("AFC-Distribution/AFC-Sequence-",Var_Periods,"-90.svg")
          svglite(title,width = 12.00,height = 8.00)
          print(plot_grid(plot_grid(plotlist = myplots_Each_Periods), ncol = 1, nrow = 1, rel_widths = c(3,3),rel_heights = c(3,3)))
          dev.off()
        }
# 
  # ASV ----------------------------------------------------------------
    data_asv_tax <- data_asv
    colnames(tax_tablemix) <- c("ASV_Id","Taxonomy")
    data_asv_tax <- merge(x = data_asv_tax, y =  tax_tablemix, by = "ASV_Id")
    data_asv_tax <- separate(data_asv_tax, Taxonomy,c("Domain","Supergroup","Division","Class","Order","Family","Genus","Species"), sep =";")
    data_asv_tax[data_asv_tax == ""] <- NA
    data_asv_tax[data_asv_tax == " "] <- NA
    data_asv_tax[data_asv_tax == "NA"] <- NA
    data_asv_tax[is.na(data_asv_tax) == TRUE] <- "Not Affiliated"
    sortDomain <- c("Eukaryota","Bacteria","Not Affiliated")
    data_asv_tax <- data_asv_tax %>% filter(Domain %in% sortDomain)
    for (i in row.names(data_asv_tax)) {
      if (data_asv_tax[i,"Domain"] == "Bacteria") { 
        data_asv_tax[i,"Supergroup"] <- "Bacteria_X"
        data_asv_tax[i,"Division"] <- "Bacteria_XX"
        data_asv_tax[i,"Class"] <- "Bacteria_XXX"
        data_asv_tax[i,"Order"] <- "Bacteria_XXXX"
        data_asv_tax[i,"Family"] <- "Bacteria_XXXXX"
        data_asv_tax[i,"Genus"] <- "Bacteria_XXXXXX"
        data_asv_tax[i,"Species"] <- "Bacteria_XXXXXXX"
      }
    }
#
    # Sort interest RangInterest ------------------------------------------------------
      data_asv_tax$RangInterest <- data_asv_tax$Supergroup
#
    # Special Zone ------------------------------------------------------------
      data_asv_tax_specialZ <- data_asv_tax
      data_asv_tax_specialZ$Period<-"Shared"
      for (i in row.names(data_asv_tax_specialZ)) { 
        if (data_asv_tax_specialZ[i,"Total04"]>0 && data_asv_tax_specialZ[i,"Total06"]==0 && data_asv_tax_specialZ[i,"Total09"]==0 && data_asv_tax_specialZ[i,"Total11"]==0) { data_asv_tax_specialZ[i,"Period"] <- "Avril"}
        if (data_asv_tax_specialZ[i,"Total06"]>0 && data_asv_tax_specialZ[i,"Total04"]==0 && data_asv_tax_specialZ[i,"Total09"]==0 && data_asv_tax_specialZ[i,"Total11"]==0) { data_asv_tax_specialZ[i,"Period"] <- "Juin"}
        if (data_asv_tax_specialZ[i,"Total09"]>0 && data_asv_tax_specialZ[i,"Total06"]==0 && data_asv_tax_specialZ[i,"Total04"]==0 && data_asv_tax_specialZ[i,"Total11"]==0) { data_asv_tax_specialZ[i,"Period"] <- "Septembre"}
        if (data_asv_tax_specialZ[i,"Total11"]>0 && data_asv_tax_specialZ[i,"Total06"]==0 && data_asv_tax_specialZ[i,"Total09"]==0 && data_asv_tax_specialZ[i,"Total04"]==0) { data_asv_tax_specialZ[i,"Period"] <- "Novembre"}
      }
      data_asv_tax_specialZ$Condition<-NULL
      for (i in row.names(data_asv_tax_specialZ)) { if (is.na(data_asv_tax_specialZ[i,"Period"])==FALSE ) {data_asv_tax_specialZ[i,"Condition"] <- paste(data_asv_tax_specialZ[i,"Zone90"],data_asv_tax_specialZ[i,"Period"],sep=" & ") }}
      data_asv_tax_specialZ <- data_asv_tax_specialZ %>% filter(is.na(Condition)==FALSE)
    ##Figure Zone
      multiZoneasv <- ggplot(data_asv_tax_specialZ, aes(y = `Dim 2`, x = `Dim 1`, color = Period)) + geom_point(size = 2) +
        geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
        geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
        stat_ellipse(aes(alpha = Period,linetype =Period),geom = "polygon",type = "norm") +
        scale_linetype_manual(values = c(2,2,2,2,0)) +
        theme(axis.title = element_text(face="bold", size=12), 
              axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
              title = element_text(face="bold", size=14),
              legend.title = element_text(face="bold"),
              legend.position = "right",
              legend.text = element_text(size=10)) +
        scale_alpha_manual(values = c(0.1,0.1,0.1,0.1,0)) +
        scale_color_manual(values = c("#F8766D","#80B69B","#7F76A7","#00A5FF","lightgrey")) +
        labs(x=Dim1asv,y=Dim2asv,color = "Period", alpha = "Period", linetype = "Period")
      print(multiZoneasv)
      ggsave("Betadisp/PeriodASV.svg", device = "svg", width = 5, height = 5)


    # All_2018 ----------------------------------------------------------
      Var_Condition <- c("Zone","Cycle","Fraction")
      myplots100_All_2018 = lapply(Var_Condition,function(Var_Condition) {
        if (Var_Condition == "Zone") { order <- c("Monimolimnion","Mixolimnion","Shared")}
        if (Var_Condition == "Cycle") { order <- c("Day","Night","Shared")}
        if (Var_Condition == "Fraction") { order <- c("Large","Small","Shared")}
        ggplot(data_asv_tax, aes(y = `Dim 2`, x = `Dim 1`, color = factor(get(Var_Condition),level=order))) + geom_point(size = 2) +
          geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
          geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
          stat_ellipse(aes(linetype = factor(get(Var_Condition),level=order), alpha = factor(get(Var_Condition),level=order)),geom = "polygon",type = "norm") +
          scale_linetype_manual(values = c(2,2,0)) +
          theme(axis.title = element_text(face="bold", size=12), 
                axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
                title = element_text(face="bold", size=14),
                legend.title = element_text(face="bold"),
                legend.position = "right",
                legend.text = element_text(size=10)) +
          scale_alpha_manual(values = c(0.1,0.1,0)) +
          scale_color_manual(values = c("#F8766D","#00A5FF","lightgrey")) +
          labs(x=Dim1asv,y=Dim2asv,color = Var_Condition, alpha = Var_Condition, linetype = Var_Condition)
      }
      )
      svglite("AFC-Distribution/AFC-ASV.svg",width = 12.00,height = 8.00)
      print(plot_grid(plot_grid(plotlist = myplots100_All_2018), ncol = 1, nrow = 1))
      dev.off()
#
    # Each Periods ------------------------------------------------------------
      for (Var_Periods in c("04","06","09","11")){ print(Var_Periods)
        assign(paste0("data_asv",Var_Periods,"_tax"),merge(x=get(paste0("data_asv",Var_Periods)),y=tax_tablemix,by = "ASV_Id"))
        Var_Condition <- c("Zone","Cycle","Fraction")
        myplots_Each_Periods = lapply(Var_Condition,function(Var_Condition) {
          if (Var_Condition == "Zone") { order <- c("Monimolimnion","Mixolimnion","Shared")}
          if (Var_Condition == "Cycle") { order <- c("Day","Night","Shared")}
          if (Var_Condition == "Fraction") { order <- c("Large","Small","Shared")}
          ggplot(get(paste0("data_asv",Var_Periods,"_tax")), aes(y = `Dim 2`, x = `Dim 1`, color = factor(get(Var_Condition),level=order))) + geom_point(size = 2) +
            geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
            geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
            stat_ellipse(aes(linetype = factor(get(Var_Condition),level=order), alpha = factor(get(Var_Condition),level=order)),geom = "polygon",type = "norm") +
            scale_linetype_manual(values = c(2,2,0)) +
            theme(axis.title = element_text(face="bold", size=12), 
                  axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
                  title = element_text(face="bold", size=14),
                  legend.title = element_text(face="bold"),
                  legend.position = "right",
                  legend.text = element_text(size=10)) +
            scale_alpha_manual(values = c(0.1,0.1,0)) +
            scale_color_manual(values = c("#F8766D","#00A5FF","lightgrey")) +
            labs(x=Dim1asv,y=Dim2asv,color = Var_Condition, alpha = Var_Condition, linetype = Var_Condition)
        }
        )
        title <- paste0("AFC-Distribution/AFC-ASV-",Var_Periods,"-100.svg")
        svglite(title,width = 12.00,height = 8.00)
        print(plot_grid(plot_grid(plotlist = myplots_Each_Periods), ncol = 1, nrow = 1, rel_widths = c(3,3),rel_heights = c(3,3)))
        dev.off()
      }
      #
# Venn diagramm -----------------------------------------------------------
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
## Cycle
  CycleVenn <- venn.diagram(x = list(data_asv_tax %>% filter(Day == 1) %>% select(ASV_Id) %>% unlist(), data_asv_tax %>% filter(Night == 2) %>% select(ASV_Id) %>% unlist()), 
                            category.names = c("Day" , "Night"),
                            filename = NULL,
                            compression = "lzw",
                            lwd = 1,
                            main = "ASVs",
                            col=c("#440154ff", '#21908dff'),
                            fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
                            cex = 0.5,
                            main.fontfamily = "sans",
                            main.fontface = "bold",
                            fontfamily = "sans",
                            cat.cex = 0.6,
                            cat.default.pos = "outer",
                            cat.pos = c(-35,35),
                            cat.fontfamily = "sans",
                            cat.fontface = "bold",
                            cat.col = c("#440154ff", '#21908dff'))
  ggsave(CycleVenn, file="Venn/Cycle.svg", device = "svg", width = 3, height = 3)
## Zone
  ZoneVenn <- venn.diagram(x = list(data_asv_tax %>% filter(Mixolimnion == 1) %>% select(ASV_Id) %>% unlist(), data_asv_tax %>% filter(Monimolimnion == 2) %>% select(ASV_Id) %>% unlist()),
                          category.names = c("Mixolimnion" , "Monimolimnion"),
                          filename = NULL,
                          compression = "lzw",
                          lwd = 1,
                          main = "ASVs",
                          col=c("#440154ff", '#21908dff'),
                          fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
                          cex = 0.5,
                          main.fontfamily = "sans",
                          main.fontface = "bold",
                          fontfamily = "sans",
                          cat.cex = 0.6,
                          cat.default.pos = "outer",
                          cat.pos = c(-35,35),
                          cat.fontfamily = "sans",
                          cat.fontface = "bold",
                          cat.col = c("#440154ff", '#21908dff'))
  ggsave(ZoneVenn, file="Venn/Zone.svg", device = "svg", width = 3, height = 3)
## Fraction
  FractionVenn <- venn.diagram(x = list(data_asv_tax %>% filter(Small == 1) %>% select(ASV_Id) %>% unlist(), data_asv_tax %>% filter(Large == 2) %>% select(ASV_Id) %>% unlist()),
                           category.names = c("Small" , "Large"),
                           filename = NULL,
                           compression = "lzw",
                           lwd = 1,
                           main = "ASVs",
                           col=c("#440154ff", '#21908dff'),
                           fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
                           cex = 0.5,
                           main.fontfamily = "sans",
                           main.fontface = "bold",
                           fontfamily = "sans",
                           cat.cex = 0.6,
                           cat.default.pos = "outer",
                           cat.pos = c(-35,35),
                           cat.fontfamily = "sans",
                           cat.fontface = "bold",
                           cat.col = c("#440154ff", '#21908dff'))
  ggsave(FractionVenn, file="Venn/Fraction.svg", device = "svg", width = 3, height = 3)
## Period
  PeriodVenn <- venn.diagram(x = list(data_asv_tax %>% filter(Total04 == 1) %>% select(ASV_Id) %>% unlist(), data_asv_tax %>% filter(Total06 == 1) %>% select(ASV_Id) %>% unlist(),data_asv_tax %>% filter(Total09 == 1) %>% select(ASV_Id) %>% unlist(),data_asv_tax %>% filter(Total11 == 1) %>% select(ASV_Id) %>% unlist()),
                               category.names = c("Avril","Juin","Septembre","Novembre"),
                               filename = NULL,
                               compression = "lzw",
                               lwd = 1,
                               main = "ASVs",
                               col=c("#440154ff", "#638fb4ff", "#0062b4ff",'#21908dff'),
                               fill = c(alpha("#440154ff",0.3), alpha('#638fb4ff',0.3), alpha('#0062b4ff',0.3),alpha('#21908dff',0.3)),
                               cex = 0.5,
                               main.fontfamily = "sans",
                               main.fontface = "bold",
                               fontfamily = "sans",
                               cat.cex = 0.6,
                               cat.default.pos = "outer",
                               #cat.pos = c(-35,35),
                               cat.fontfamily = "sans",
                               cat.fontface = "bold",
                               cat.col = c("#440154ff", "#638fb4ff", "#0062b4ff",'#21908dff'))
  ggsave(PeriodVenn, file="Venn/Period.svg", device = "svg", width = 3, height = 3)
#
# HIST  --------------------------------------------------------------------
  # Color for Tax ------------------------------------------------------------
  color_tax <- data_seq_tax %>% dplyr::select(Supergroup,Division) %>% distinct(Supergroup,Division)
  color_tax[,"palet"] <- ""
  color_tax[,"Label"] <- color_tax[,"Division"]
  for (i in row.names(color_tax)){
    if (color_tax[i,"Label"] == "Not Affiliated" && color_tax[i,"Supergroup"]!="Not Affiliated") { color_tax[i,"Label"] <- paste("Unaffiliated",color_tax[i,"Supergroup"])}
    if (grepl(x = color_tax[i,"Division"], pattern = "_X") == TRUE) { color_tax[i,"Label"] <- paste("Unaffiliated",color_tax[i,"Supergroup"])}
  }
  color_tax
  color_tax <- color_tax %>% arrange(Supergroup,Label)
  
  for (i in 1:length(unique(color_tax$Supergroup))) {
    B_courant <- unique(color_tax$Supergroup)[i]
    B_countLabel <- nrow(unique(color_tax %>% filter(Supergroup == B_courant) %>% dplyr::select('Label')))
    #palet_liste <- c("red_material","pink_material","purple_material","deep_purple_material","indigo_material","blue_material","light_blue_material","cyan_material","teal_material","green_material","light_green_material","lime_material","yellow_material","amber_material","orange_material","deep_orange_material","brown_material","grey_material","blue_grey_material")
    palet_liste <- c("red_material","deep_purple_material","blue_material","light_blue_material","light_green_material","yellow_material","orange_material","grey_material","pink_material","cyan_material","green_material","lime_material","indigo_material","purple_material","amber_material","blue_grey_material")
    palet_name <- paste0("ggsci::",palet_liste[i])
    print(palet_name)
    palet_courante <- c(paletteer_d(palet_name,direction = -1,n=10),"#fae6eb")[seq(1,floor(B_countLabel+(B_countLabel/2))+2,2)]
    print(B_courant)
    print(B_countLabel)
    print(palet_courante)
    for (j in 1:B_countLabel) {
      Label_courant <- unique(color_tax %>% filter(Supergroup==B_courant) %>% dplyr::select(Label))[[1]][j]
      print(Label_courant)
      color_tax[,"palet"][color_tax[,"Label"]==Label_courant] <- palet_courante[j]
    }
  }
  color_tax <- color_tax %>% arrange(Supergroup,Label)
  color_tax_table <- color_tax %>% dplyr::select(palet,Label) %>% distinct()
  color_tax_table_order <- color_tax_table$Label
  paletspe <- color_tax_table$palet
  
  
#
  # Zone X Periods ----------------------------------------------------------
    for (type in c("Sequence","ASV")) {
      if (type == "ASV") { input_data_table <- asv_mat_rare }
      if (type == "Sequence") { input_data_table <- seq_mat_rare }
      for (zonexs in c("Mixolimnion","Monimolimnion")) {
        for (periodxs in c("04","06","09","11")) {
          tableI <- paste(zonexs,periodxs,sep="X")
          assign(tableI,samples_df %>% filter(Zone == zonexs) %>% filter(Period == periodxs) %>% filter(Replicat == 1))
          assign(tableI,get(tableI)$Amplicon)
          assign(paste(tableI,"seq",sep="_"),as.data.frame(rowSums(input_data_table %>% select(all_of(get(tableI))))))
          tableIt <- get(paste(tableI,"seq",sep="_"))
          colnames(tableIt) <- paste(zonexs,periodxs,sep="X")
          tableIt$ASV_Id <- row.names(tableIt)
          assign(paste(tableI,"seq",sep="_"),tableIt)
        }
      }
  ## Merging
    ZoneX04 <- merge(x = MixolimnionX04_seq,y = MonimolimnionX04_seq, by = "ASV_Id")
    ZoneX06 <- merge(x = MixolimnionX06_seq,y = MonimolimnionX06_seq, by = "ASV_Id")
    ZoneX09 <- merge(x = MixolimnionX09_seq,y = MonimolimnionX09_seq, by = "ASV_Id")
    ZoneX11 <- merge(x = MixolimnionX11_seq,y = MonimolimnionX11_seq, by = "ASV_Id")
    ZoneXi <- merge(x = ZoneX04,y = ZoneX06, by = "ASV_Id")
    ZoneXj <- merge(x = ZoneXi,y = ZoneX09, by = "ASV_Id")
    ZoneXPeriods <- merge(x = ZoneXj,y = ZoneX11, by = "ASV_Id")
    if (type == "ASV") { ZoneXPeriods[,2:9][ZoneXPeriods[,2:9]>1]<-1}
    ZoneXPeriods_tax <- merge(data_seq_tax %>% select(ASV_Id,RangInterest,Division), ZoneXPeriods, by = "ASV_Id")
    ZoneXPeriods_tax <- ZoneXPeriods_tax %>% select(-ASV_Id)
    for (i in rownames(ZoneXPeriods_tax)) { if (ZoneXPeriods_tax[i,"Division"]=="Not Affiliated" && ZoneXPeriods_tax[i,"RangInterest"]!="Not Affiliated") { ZoneXPeriods_tax[i,"Division"] <- paste("Unaffiliated",ZoneXPeriods_tax[i,"RangInterest"],sep=" ")}
      if (ZoneXPeriods_tax[i,"Division"]=="Alveolata_X" ) { ZoneXPeriods_tax[i,"Division"] <- paste("Unaffiliated",ZoneXPeriods_tax[i,"RangInterest"],sep=" ")}}
    ZoneXPeriods_tax <- ZoneXPeriods_tax %>% group_by(RangInterest,Division) %>% summarise_all(sum)
    Orderspe <- ZoneXPeriods_tax$Division
    ZoneXPeriods_tax_melt <-melt(ZoneXPeriods_tax,id.vars = c("RangInterest","Division"))
    ZoneXPeriods_tax_melt$Sum <- 0
    for (i in row.names(ZoneXPeriods_tax_melt)) { var <- ZoneXPeriods_tax_melt[i,"variable"]
      ZoneXPeriods_tax_melt[i,"Sum"] <- sum(ZoneXPeriods_tax_melt %>% filter(variable == var) %>% select(value))
      ZoneXPeriods_tax_melt[i,"Proportion"] <- ZoneXPeriods_tax_melt[i,"value"]*100/ZoneXPeriods_tax_melt[i,"Sum"]}
      ZoneXPeriods_tax_melt <- ZoneXPeriods_tax_melt %>% mutate()
      ZoneXPeriods_tax_melt <- separate(ZoneXPeriods_tax_melt,"variable",c("Zone","Periods"),sep="X")
      assign(paste("ZoneXPeriods_tax_melt",type,sep = "_"),ZoneXPeriods_tax_melt)
    }
  ## Hist
  ### Sequence
    svglite("Hist-Taxomy/ZoneXPeriods-Sequence-Total-bar.svg",width = 7.00,height = 9.00)
    Total_ZoneXPeriods_Sequence_figa <- ggplot(ZoneXPeriods_tax_melt_Sequence, mapping = aes(y= Proportion, x = Periods,fill = factor(Division,level=Orderspe), group = factor(Division,level=Orderspe)), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity",color="black") +
      facet_grid(.~Zone,scales="free") + scale_fill_manual(values = paletspe) + guides(fill=guide_legend(ncol=1)) +
      labs(x="Periods",y="Sequences %",fill="Divisions") + theme_unique_art()
    print(Total_ZoneXPeriods_Sequence_figa)
    dev.off()
    write.table(ZoneXPeriods_tax_melt_Sequence,file="dataTables/ZoneXPeriods-Sequence-Total-Table.csv",row.names = F, quote = F, sep ="\t")
  ### ASV
    svglite("Hist-Taxomy/ZoneXPeriods-ASV-Total-bar.svg",width = 7.00,height = 9.00)
    Total_ZoneXPeriods_ASV_figa <- ggplot(ZoneXPeriods_tax_melt_ASV, mapping = aes(y= Proportion, x = Periods, fill = factor(Division,level=Orderspe), group = factor(Division,level=Orderspe)), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity",color="black") +
      facet_grid(.~Zone,scales="free") + scale_fill_manual(values = paletspe) + guides(fill=guide_legend(ncol=1)) +
      labs(x="Periods",y="ASVs %",fill="Divisions") + geom_label(aes(y = 106,label = paste(Sum," ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") + theme_unique_art()
    print(Total_ZoneXPeriods_ASV_figa)
    dev.off()
    write.table(ZoneXPeriods_tax_melt_ASV,file="dataTables/ZoneXPeriods-ASV-Total-Table.csv",row.names = F, quote = F, sep ="\t")
  ## Area
  ### Sequence
    svglite("Hist-Taxomy/ZoneXPeriods-Sequence-Total-area.svg",width = 7.00,height = 9.00)
    Total_ZoneXPeriods_Sequence_figb <- ggplot(ZoneXPeriods_tax_melt_Sequence, mapping = aes(y= Proportion, x = Periods,fill = factor(Division,level=Orderspe), group = factor(Division,level=Orderspe)), Rowv = NA, col = colMain, scale = "column") + geom_area(stat="identity",color="black",alpha=0.95) +
      facet_grid(.~Zone,scales="free") + scale_fill_manual(values = paletspe) + guides(fill=guide_legend(ncol=1)) +
      labs(x="Periods",y="Sequences %",fill="Divisions") + theme_unique_art()
    print(Total_ZoneXPeriods_Sequence_figb)
    dev.off()
  ### ASV
    svglite("Hist-Taxomy/ZoneXPeriods-ASV-Total-area.svg",width = 7.00,height = 9.00)
    Total_ZoneXPeriods_ASV_figb <- ggplot(ZoneXPeriods_tax_melt_ASV, mapping = aes(y= Proportion, x = Periods, fill = factor(Division,level=Orderspe), group = factor(Division,level=Orderspe)), Rowv = NA, col = colMain, scale = "column") + geom_area(stat="identity",color="black",alpha=0.95) +
      facet_grid(.~Zone,scales="free") + scale_fill_manual(values = paletspe) + guides(fill=guide_legend(ncol=1)) +
      labs(x="Periods",y="ASVs %",fill="Divisions") + theme_unique_art()
    print(Total_ZoneXPeriods_ASV_figb)
    dev.off()
  ## Coplot
    legendSequenceArea <- get_legend(Total_ZoneXPeriods_ASV_figb)
    Total_ZoneXPeriods_ASV_figbj <- Total_ZoneXPeriods_ASV_figb + theme(legend.position = "none", axis.text.x = element_blank(),axis.title.x = element_blank()) + scale_x_discrete(position="top")
    Total_ZoneXPeriods_Sequence_figbj <- Total_ZoneXPeriods_Sequence_figb + theme(legend.position = "none",strip.text = element_blank())
  ### plot
    xmi_plot <- plot_grid(Total_ZoneXPeriods_ASV_figbj,Total_ZoneXPeriods_Sequence_figbj, ncol = 1, nrow = 2, rel_widths = c(3),rel_heights = c(3,3))
    xmi_plot <- plot_grid(xmi_plot,legendSequenceArea, ncol = 2, nrow = 1, rel_widths = c(5,2),rel_heights = c(3))
    svglite("Hist-Taxomy/ZoneXPeriods-ASV&Seq-Total-area.svg",width = 9.00,height = 9.00)
    print(xmi_plot)
    dev.off()
#
  # ASV Total ---------------------------------------------------------------------
    # Cycle -------------------------------------------------------------------
    ## Day
      totalDayASV <- data_asv_tax %>% select(ASV_Id,TotalDay,RangInterest)
      row.names(totalDayASV)<-totalDayASV$ASV_Id ; totalDayASV <- totalDayASV %>% select(-ASV_Id)
      totalDayASV <- totalDayASV %>% group_by(RangInterest) %>% summarise_all(sum)
      totalDayASV$TotalDay <- totalDayASV$TotalDay*100/sum(totalDayASV$TotalDay)
      totalDayASV <- totalDayASV %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(TotalDay) - 0.5*TotalDay)
      totalDayASV$label <- paste(round(totalDayASV$TotalDay,1), "%", sep = "")
      for (i in rownames(totalDayASV)) {
        if (totalDayASV[i,"label"] == "0%") { totalDayASV[i,"label"] <- NA}}
      for (i in rownames(totalDayASV)) {
        if (is.na(totalDayASV[i,"label"]) == FALSE) { totalDayASV[i,"label"] <- paste(totalDayASV[i,"RangInterest"]," : ",totalDayASV[i,"label"], sep = "")}}
      totalDayASV$Cycle<- rep("Day", each = nrow(totalDayASV))
    ## Night
      totalNightASV <- data_asv_tax %>% select(ASV_Id,TotalNight,RangInterest)
      row.names(totalNightASV)<-totalNightASV$ASV_Id ; totalNightASV <- totalNightASV %>% select(-ASV_Id)
      totalNightASV <- totalNightASV %>% group_by(RangInterest) %>% summarise_all(sum)
      totalNightASV$TotalNight <- totalNightASV$TotalNight*100/sum(totalNightASV$TotalNight)
      totalNightASV <- totalNightASV %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(TotalNight) - 0.5*TotalNight)
      totalNightASV$label <- paste(round(totalNightASV$TotalNight,1), "%", sep = "")
      for (i in rownames(totalNightASV)) {
        if (totalNightASV[i,"label"] == "0%") { totalNightASV[i,"label"] <- NA}}
      for (i in rownames(totalNightASV)) {
        if (is.na(totalNightASV[i,"label"]) == FALSE) { totalNightASV[i,"label"] <- paste(totalNightASV[i,"RangInterest"]," : ",totalNightASV[i,"label"], sep = "")}}
      totalNightASV$Cycle<- rep("Night", each = nrow(totalNightASV))
    ## Cycle
      colnames(totalDayASV)[2]  <- "value"
      totalDayASV$Sum <- rep(0, each = nrow(totalDayASV))
      totalDayASV$Count <- rep(0, each = nrow(totalDayASV))
      for (i in totalDayASV$RangInterest) { totalDayASV$Count[which(totalDayASV$RangInterest == i)] <- sum(data_asv_tax  %>% filter(RangInterest == i) %>% select(TotalDay))}
      totalDayASV$Sum <- sum(totalDayASV$Count)
      colnames(totalNightASV)[2]  <- "value"
      totalNightASV$Sum <- rep(0, each = nrow(totalNightASV))
      totalNightASV$Count <- rep(0, each = nrow(totalNightASV))
      for (i in totalNightASV$RangInterest) { totalNightASV$Count[which(totalNightASV$RangInterest == i)] <- sum(data_asv_tax  %>% filter(RangInterest == i) %>% select(TotalNight))}
      totalNightASV$Sum <- sum(totalNightASV$Count)
      totalCycleASV <- rbind(totalDayASV,totalNightASV)
    ## Figure
      ay <- ggplot(totalCycleASV, mapping = aes(y= value, x = Cycle, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + scale_fill_manual(values = palette) + 
        geom_label(aes(y = 106,label = paste(Sum," ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
      legendASV <- get_legend(ay)
      ay <- ay + labs(x="Cycles",y="ASVs (%)") + theme(legend.position = "none")
      print(ay)
#
    # Zone -------------------------------------------------------------------
    ## Mixolimnion
      totalMixolimnionASV <- data_asv_tax %>% select(ASV_Id,TotalMixolimnion,RangInterest)
      row.names(totalMixolimnionASV)<-totalMixolimnionASV$ASV_Id ; totalMixolimnionASV <- totalMixolimnionASV %>% select(-ASV_Id)
      totalMixolimnionASV <- totalMixolimnionASV %>% group_by(RangInterest) %>% summarise_all(sum)
      totalMixolimnionASV$TotalMixolimnion <- totalMixolimnionASV$TotalMixolimnion*100/sum(totalMixolimnionASV$TotalMixolimnion)
      totalMixolimnionASV <- totalMixolimnionASV %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(TotalMixolimnion) - 0.5*TotalMixolimnion)
      totalMixolimnionASV$label <- paste(round(totalMixolimnionASV$TotalMixolimnion,1), "%", sep = "")
      for (i in rownames(totalMixolimnionASV)) {
        if (totalMixolimnionASV[i,"label"] == "0%") { totalMixolimnionASV[i,"label"] <- NA}}
      for (i in rownames(totalMixolimnionASV)) {
        if (is.na(totalMixolimnionASV[i,"label"]) == FALSE) { totalMixolimnionASV[i,"label"] <- paste(totalMixolimnionASV[i,"RangInterest"]," : ",totalMixolimnionASV[i,"label"], sep = "")}}
      totalMixolimnionASV$Zone<- rep("Mixolimnion", each = nrow(totalMixolimnionASV))
    ## Monimolimnion
      totalMonimolimnionASV <- data_asv_tax %>% select(ASV_Id,TotalMonimolimnion,RangInterest)
      row.names(totalMonimolimnionASV)<-totalMonimolimnionASV$ASV_Id ; totalMonimolimnionASV <- totalMonimolimnionASV %>% select(-ASV_Id)
      totalMonimolimnionASV <- totalMonimolimnionASV %>% group_by(RangInterest) %>% summarise_all(sum)
      totalMonimolimnionASV$TotalMonimolimnion <- totalMonimolimnionASV$TotalMonimolimnion*100/sum(totalMonimolimnionASV$TotalMonimolimnion)
      totalMonimolimnionASV <- totalMonimolimnionASV %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(TotalMonimolimnion) - 0.5*TotalMonimolimnion)
      totalMonimolimnionASV$label <- paste(round(totalMonimolimnionASV$TotalMonimolimnion,1), "%", sep = "")
      for (i in rownames(totalMonimolimnionASV)) {
        if (totalMonimolimnionASV[i,"label"] == "0%") { totalMonimolimnionASV[i,"label"] <- NA}}
      for (i in rownames(totalMonimolimnionASV)) {
        if (is.na(totalMonimolimnionASV[i,"label"]) == FALSE) { totalMonimolimnionASV[i,"label"] <- paste(totalMonimolimnionASV[i,"RangInterest"]," : ",totalMonimolimnionASV[i,"label"], sep = "")}}
      totalMonimolimnionASV$Zone<- rep("Monimolimnion", each = nrow(totalMonimolimnionASV))
    ## Zone
      colnames(totalMixolimnionASV)[2]  <- "value"
      totalMixolimnionASV$Sum <- rep(0, each = nrow(totalMixolimnionASV))
      totalMixolimnionASV$Count <- rep(0, each = nrow(totalMixolimnionASV))
      for (i in totalMixolimnionASV$RangInterest) { totalMixolimnionASV$Count[which(totalMixolimnionASV$RangInterest == i)] <- sum(data_asv_tax  %>% filter(RangInterest == i) %>% select(TotalMixolimnion))}
      totalMixolimnionASV$Sum <- sum(totalMixolimnionASV$Count)
      colnames(totalMonimolimnionASV)[2]  <- "value"
      totalMonimolimnionASV$Sum <- rep(0, each = nrow(totalMonimolimnionASV))
      totalMonimolimnionASV$Count <- rep(0, each = nrow(totalMonimolimnionASV))
      for (i in totalMonimolimnionASV$RangInterest) { totalMonimolimnionASV$Count[which(totalMonimolimnionASV$RangInterest == i)] <- sum(data_asv_tax  %>% filter(RangInterest == i) %>% select(TotalMonimolimnion))}
      totalMonimolimnionASV$Sum <- sum(totalMonimolimnionASV$Count)
      totalZoneASV <- rbind(totalMixolimnionASV,totalMonimolimnionASV)
    ## plot
      by <- ggplot(totalZoneASV, mapping = aes(y= value, x = Zone, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + scale_fill_manual(values = palette) + 
        geom_label(aes(y = 106,label = paste(Sum," ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
      by <- by + labs(x="Zones",y="ASVs (%)") + theme(legend.position = "none")
      print(by)
#
    # Fraction -------------------------------------------------------------------
    ## Small
      totalSmallASV <- data_asv_tax %>% select(ASV_Id,TotalSmall,RangInterest)
      row.names(totalSmallASV)<-totalSmallASV$ASV_Id ; totalSmallASV <- totalSmallASV %>% select(-ASV_Id)
      totalSmallASV <- totalSmallASV %>% group_by(RangInterest) %>% summarise_all(sum)
      totalSmallASV$TotalSmall <- totalSmallASV$TotalSmall*100/sum(totalSmallASV$TotalSmall)
      totalSmallASV <- totalSmallASV %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(TotalSmall) - 0.5*TotalSmall)
      totalSmallASV$label <- paste(round(totalSmallASV$TotalSmall,1), "%", sep = "")
      for (i in rownames(totalSmallASV)) {
        if (totalSmallASV[i,"label"] == "0%") { totalSmallASV[i,"label"] <- NA}}
      for (i in rownames(totalSmallASV)) {
        if (is.na(totalSmallASV[i,"label"]) == FALSE) { totalSmallASV[i,"label"] <- paste(totalSmallASV[i,"RangInterest"]," : ",totalSmallASV[i,"label"], sep = "")}}
      totalSmallASV$Fraction<- rep("Small", each = nrow(totalSmallASV))
    ## Large
      totalLargeASV <- data_asv_tax %>% select(ASV_Id,TotalLarge,RangInterest)
      row.names(totalLargeASV)<-totalLargeASV$ASV_Id ; totalLargeASV <- totalLargeASV %>% select(-ASV_Id)
      totalLargeASV <- totalLargeASV %>% group_by(RangInterest) %>% summarise_all(sum)
      totalLargeASV$TotalLarge <- totalLargeASV$TotalLarge*100/sum(totalLargeASV$TotalLarge)
      totalLargeASV <- totalLargeASV %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(TotalLarge) - 0.5*TotalLarge)
      totalLargeASV$label <- paste(round(totalLargeASV$TotalLarge,1), "%", sep = "")
      for (i in rownames(totalLargeASV)) {
        if (totalLargeASV[i,"label"] == "0%") { totalLargeASV[i,"label"] <- NA}}
      for (i in rownames(totalLargeASV)) {
        if (is.na(totalLargeASV[i,"label"]) == FALSE) { totalLargeASV[i,"label"] <- paste(totalLargeASV[i,"RangInterest"]," : ",totalLargeASV[i,"label"], sep = "")}}
      totalLargeASV$Fraction<- rep("Large", each = nrow(totalLargeASV))
    ## Fraction
      colnames(totalSmallASV)[2]  <- "value"
      totalSmallASV$Sum <- rep(0, each = nrow(totalSmallASV))
      totalSmallASV$Count <- rep(0, each = nrow(totalSmallASV))
      for (i in totalSmallASV$RangInterest) { totalSmallASV$Count[which(totalSmallASV$RangInterest == i)] <- sum(data_asv_tax  %>% filter(RangInterest == i) %>% select(TotalSmall))}
      totalSmallASV$Sum <- sum(totalSmallASV$Count)
      colnames(totalLargeASV)[2]  <- "value"
      totalLargeASV$Sum <- rep(0, each = nrow(totalLargeASV))
      totalLargeASV$Count <- rep(0, each = nrow(totalLargeASV))
      for (i in totalLargeASV$RangInterest) { totalLargeASV$Count[which(totalLargeASV$RangInterest == i)] <- sum(data_asv_tax  %>% filter(RangInterest == i) %>% select(TotalLarge))}
      totalLargeASV$Sum <- sum(totalLargeASV$Count)
      totalFractionASV <- rbind(totalSmallASV,totalLargeASV)
    ## Plot
      cy <- ggplot(totalFractionASV, mapping = aes(y= value, x = Fraction, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + scale_fill_manual(values = palette)+ 
        geom_label(aes(y = 106,label = paste(Sum," ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
      cy <- cy + labs(x="Fractions",y="ASVs (%)") + theme(legend.position = "none")
      print(cy)
#
    # Period -------------------------------------------------------------------
    ## Avril
      total04ASV <- data_asv_tax %>% select(ASV_Id,Total04,RangInterest)
      row.names(total04ASV)<-total04ASV$ASV_Id ; total04ASV <- total04ASV %>% select(-ASV_Id)
      total04ASV <- total04ASV %>% group_by(RangInterest) %>% summarise_all(sum)
      total04ASV$Total04 <- total04ASV$Total04*100/sum(total04ASV$Total04)
      total04ASV <- total04ASV %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(Total04) - 0.5*Total04)
      total04ASV$label <- paste(round(total04ASV$Total04,1), "%", sep = "")
      for (i in rownames(total04ASV)) {
        if (total04ASV[i,"label"] == "0%") { total04ASV[i,"label"] <- NA}}
      for (i in rownames(total04ASV)) {
        if (is.na(total04ASV[i,"label"]) == FALSE) { total04ASV[i,"label"] <- paste(total04ASV[i,"RangInterest"]," : ",total04ASV[i,"label"], sep = "")}}
      total04ASV$Period<- rep("04", each = nrow(total04ASV))
    ## Juin
      total06ASV <- data_asv_tax %>% select(ASV_Id,Total06,RangInterest)
      row.names(total06ASV)<-total06ASV$ASV_Id ; total06ASV <- total06ASV %>% select(-ASV_Id)
      total06ASV <- total06ASV %>% group_by(RangInterest) %>% summarise_all(sum)
      total06ASV$Total06 <- total06ASV$Total06*100/sum(total06ASV$Total06)
      total06ASV <- total06ASV %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(Total06) - 0.5*Total06)
      total06ASV$label <- paste(round(total06ASV$Total06,1), "%", sep = "")
      for (i in rownames(total06ASV)) {
        if (total06ASV[i,"label"] == "0%") { total06ASV[i,"label"] <- NA}}
      for (i in rownames(total06ASV)) {
        if (is.na(total06ASV[i,"label"]) == FALSE) { total06ASV[i,"label"] <- paste(total06ASV[i,"RangInterest"]," : ",total06ASV[i,"label"], sep = "")}}
      total06ASV$Period<- rep("06", each = nrow(total06ASV))
    ## Septembre
      total09ASV <- data_asv_tax %>% select(ASV_Id,Total09,RangInterest)
      row.names(total09ASV)<-total09ASV$ASV_Id ; total09ASV <- total09ASV %>% select(-ASV_Id)
      total09ASV <- total09ASV %>% group_by(RangInterest) %>% summarise_all(sum)
      total09ASV$Total09 <- total09ASV$Total09*100/sum(total09ASV$Total09)
      total09ASV <- total09ASV %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(Total09) - 0.5*Total09)
      total09ASV$label <- paste(round(total09ASV$Total09,1), "%", sep = "")
      for (i in rownames(total09ASV)) {
        if (total09ASV[i,"label"] == "0%") { total09ASV[i,"label"] <- NA}}
      for (i in rownames(total09ASV)) {
        if (is.na(total09ASV[i,"label"]) == FALSE) { total09ASV[i,"label"] <- paste(total09ASV[i,"RangInterest"]," : ",total09ASV[i,"label"], sep = "")}}
      total09ASV$Period<- rep("09", each = nrow(total09ASV))
    ## Octobre
      total11ASV <- data_asv_tax %>% select(ASV_Id,Total11,RangInterest)
      row.names(total11ASV)<-total11ASV$ASV_Id ; total11ASV <- total11ASV %>% select(-ASV_Id)
      total11ASV <- total11ASV %>% group_by(RangInterest) %>% summarise_all(sum)
      total11ASV$Total11 <- total11ASV$Total11*100/sum(total11ASV$Total11)
      total11ASV <- total11ASV %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(Total11) - 0.5*Total11)
      total11ASV$label <- paste(round(total11ASV$Total11,1), "%", sep = "")
      for (i in rownames(total11ASV)) {
        if (total11ASV[i,"label"] == "0%") { total11ASV[i,"label"] <- NA}}
      for (i in rownames(total11ASV)) {
        if (is.na(total11ASV[i,"label"]) == FALSE) { total11ASV[i,"label"] <- paste(total11ASV[i,"RangInterest"]," : ",total11ASV[i,"label"], sep = "")}}
      total11ASV$Period<- rep("11", each = nrow(total11ASV))
    ## Period
      colnames(total04ASV)[2]  <- "value"
      total04ASV$Sum <- rep(0, each = nrow(total04ASV))
      total04ASV$Count <- rep(0, each = nrow(total04ASV))
      for (i in total04ASV$RangInterest) { total04ASV$Count[which(total04ASV$RangInterest == i)] <- sum(data_asv_tax  %>% filter(RangInterest == i) %>% select(Total04))}
      total04ASV$Sum <- sum(total04ASV$Count)
      colnames(total06ASV)[2]  <- "value"
      total06ASV$Sum <- rep(0, each = nrow(total06ASV))
      total06ASV$Count <- rep(0, each = nrow(total06ASV))
      for (i in total06ASV$RangInterest) { total06ASV$Count[which(total06ASV$RangInterest == i)] <- sum(data_asv_tax  %>% filter(RangInterest == i) %>% select(Total06))}
      total06ASV$Sum <- sum(total06ASV$Count)
      colnames(total09ASV)[2]  <- "value"
      total09ASV$Sum <- rep(0, each = nrow(total09ASV))
      total09ASV$Count <- rep(0, each = nrow(total09ASV))
      for (i in total09ASV$RangInterest) { total09ASV$Count[which(total09ASV$RangInterest == i)] <- sum(data_asv_tax  %>% filter(RangInterest == i) %>% select(Total09))}
      total09ASV$Sum <- sum(total09ASV$Count)
      colnames(total11ASV)[2]  <- "value"
      total11ASV$Sum <- rep(0, each = nrow(total11ASV))
      total11ASV$Count <- rep(0, each = nrow(total11ASV))
      for (i in total11ASV$RangInterest) { total11ASV$Count[which(total11ASV$RangInterest == i)] <- sum(data_asv_tax  %>% filter(RangInterest == i) %>% select(Total11))}
      total11ASV$Sum <- sum(total11ASV$Count)
      totalPeriodASV <- rbind(total04ASV,total06ASV,total09ASV,total11ASV)
    ## Plot
      dy <- ggplot(totalPeriodASV, mapping = aes(y= value, x = Period, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + scale_fill_manual(values = palette)+ 
        geom_label(aes(y = 106,label = paste(Sum," ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
      dy <- dy + theme(legend.position = "none") + labs(x="Periods",y="ASVs (%)")
      print(dy)
#
    # Coplot -------------------------------------------------------------------
      svglite("Hist-Taxomy/ASV-Total.svg",width = 12.00,height = 6.00)
      b_plot <- plot_grid(ay,by,cy,dy,legendASV, ncol = 5, nrow = 1, rel_widths = c(3,3,3,5,3),rel_heights = c(3))
      print(b_plot)
      dev.off()
  # ASV Only ---------------------------------------------------------------------
    # Cycles -------------------------------------------------------------------
    ## Day
      onlyDayASV <- data_asv_tax %>% filter(Cycle == "Day") %>% select(ASV_Id,TotalDay,RangInterest)
      row.names(onlyDayASV) <- onlyDayASV$ASV_Id ; onlyDayASV <- onlyDayASV %>% select(-ASV_Id)
      onlyDayASV <- onlyDayASV %>% group_by(RangInterest) %>% summarise_all(sum)
      onlyDayASV$TotalDay <- onlyDayASV$TotalDay*100/sum(onlyDayASV$TotalDay)
      onlyDayASV <- onlyDayASV %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(TotalDay) - 0.5*TotalDay)
      onlyDayASV$label <- paste(round(onlyDayASV$TotalDay,1), "%", sep = "")
      for (i in rownames(onlyDayASV)) {
        if (onlyDayASV[i,"label"] == "0%") { onlyDayASV[i,"label"] <- NA}}
      for (i in rownames(onlyDayASV)) {
        if (is.na(onlyDayASV[i,"label"]) == FALSE) { onlyDayASV[i,"label"] <- paste(onlyDayASV[i,"RangInterest"]," : ",onlyDayASV[i,"label"], sep = "")}}
      onlyDayASV$Cycle<- rep("Day", each = nrow(onlyDayASV))
    ## Night
      onlyNightASV <- data_asv_tax %>% filter(Cycle == "Night") %>% select(ASV_Id,TotalNight,RangInterest)
      row.names(onlyNightASV)<-onlyNightASV$ASV_Id ; onlyNightASV <- onlyNightASV %>% select(-ASV_Id)
      onlyNightASV <- onlyNightASV %>% group_by(RangInterest) %>% summarise_all(sum)
      onlyNightASV$TotalNight <- onlyNightASV$TotalNight*100/sum(onlyNightASV$TotalNight)
      onlyNightASV <- onlyNightASV %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(TotalNight) - 0.5*TotalNight)
      onlyNightASV$label <- paste(round(onlyNightASV$TotalNight,1), "%", sep = "")
      for (i in rownames(onlyNightASV)) {
        if (onlyNightASV[i,"label"] == "0%") { onlyNightASV[i,"label"] <- NA}}
      for (i in rownames(onlyNightASV)) {
        if (is.na(onlyNightASV[i,"label"]) == FALSE) { onlyNightASV[i,"label"] <- paste(onlyNightASV[i,"RangInterest"]," : ",onlyNightASV[i,"label"], sep = "")}}
      onlyNightASV$Cycle<- rep("Night", each = nrow(onlyNightASV))
    ## Cycle
      colnames(onlyDayASV)[2]  <- "value"
      onlyDayASV$Sum <- rep(0, each = nrow(onlyDayASV))
      onlyDayASV$Count <- rep(0, each = nrow(onlyDayASV))
      for (i in onlyDayASV$RangInterest) { onlyDayASV$Count[which(onlyDayASV$RangInterest == i)] <- nrow(data_asv_tax %>% filter(Cycle == "Day") %>% filter(RangInterest == i))}
      onlyDayASV$Sum <- sum(onlyDayASV$Count)
      colnames(onlyNightASV)[2]  <- "value"
      onlyNightASV$Sum <- rep(0, each = nrow(onlyNightASV))
      onlyNightASV$Count <- rep(0, each = nrow(onlyNightASV))
      for (i in onlyNightASV$RangInterest) { onlyNightASV$Count[which(onlyNightASV$RangInterest == i)] <- nrow(data_asv_tax %>% filter(Cycle == "Night") %>% filter(RangInterest == i))}
      onlyNightASV$Sum <- sum(onlyNightASV$Count)
      onlyCycleASV <- rbind(onlyDayASV,onlyNightASV)
    ## plot
      az <- ggplot(onlyCycleASV, mapping = aes(y= value, x = Cycle, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
        geom_label(aes(y = 106,label = paste(Sum," ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") + scale_fill_manual(values = palette)
      az <- az + labs(x="Cycles",y="ASVs (%)")
      print(az)
#
    # Zone -------------------------------------------------------------------
    ## Mixolimnion
      onlyMixolimnionASV <- data_asv_tax %>% filter(Zone == "Mixolimnion") %>% select(ASV_Id,TotalMixolimnion,RangInterest)
      row.names(onlyMixolimnionASV)<-onlyMixolimnionASV$ASV_Id ; onlyMixolimnionASV <- onlyMixolimnionASV %>% select(-ASV_Id)
      onlyMixolimnionASV <- onlyMixolimnionASV %>% group_by(RangInterest) %>% summarise_all(sum)
      onlyMixolimnionASV$TotalMixolimnion <- onlyMixolimnionASV$TotalMixolimnion*100/sum(onlyMixolimnionASV$TotalMixolimnion)
      onlyMixolimnionASV <- onlyMixolimnionASV %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(TotalMixolimnion) - 0.5*TotalMixolimnion)
      onlyMixolimnionASV$label <- paste(round(onlyMixolimnionASV$TotalMixolimnion,1), "%", sep = "")
      for (i in rownames(onlyMixolimnionASV)) {
        if (onlyMixolimnionASV[i,"label"] == "0%") { onlyMixolimnionASV[i,"label"] <- NA}}
      for (i in rownames(onlyMixolimnionASV)) {
        if (is.na(onlyMixolimnionASV[i,"label"]) == FALSE) { onlyMixolimnionASV[i,"label"] <- paste(onlyMixolimnionASV[i,"RangInterest"]," : ",onlyMixolimnionASV[i,"label"], sep = "")}}
      onlyMixolimnionASV$Zone<- rep("Mixolimnion", each = nrow(onlyMixolimnionASV))
    ## Monimolimnion
      onlyMonimolimnionASV <- data_asv_tax %>% filter(Zone == "Monimolimnion") %>% select(ASV_Id,TotalMonimolimnion,RangInterest)
      row.names(onlyMonimolimnionASV)<-onlyMonimolimnionASV$ASV_Id ; onlyMonimolimnionASV <- onlyMonimolimnionASV %>% select(-ASV_Id)
      onlyMonimolimnionASV <- onlyMonimolimnionASV %>% group_by(RangInterest) %>% summarise_all(sum)
      onlyMonimolimnionASV$TotalMonimolimnion <- onlyMonimolimnionASV$TotalMonimolimnion*100/sum(onlyMonimolimnionASV$TotalMonimolimnion)
      onlyMonimolimnionASV <- onlyMonimolimnionASV %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(TotalMonimolimnion) - 0.5*TotalMonimolimnion)
      onlyMonimolimnionASV$label <- paste(round(onlyMonimolimnionASV$TotalMonimolimnion,1), "%", sep = "")
      for (i in rownames(onlyMonimolimnionASV)) {
        if (onlyMonimolimnionASV[i,"label"] == "0%") { onlyMonimolimnionASV[i,"label"] <- NA}}
      for (i in rownames(onlyMonimolimnionASV)) {
        if (is.na(onlyMonimolimnionASV[i,"label"]) == FALSE) { onlyMonimolimnionASV[i,"label"] <- paste(onlyMonimolimnionASV[i,"RangInterest"]," : ",onlyMonimolimnionASV[i,"label"], sep = "")}}
      onlyMonimolimnionASV$Zone<- rep("Monimolimnion", each = nrow(onlyMonimolimnionASV))
    ## Zone
      colnames(onlyMixolimnionASV)[2]  <- "value"
      onlyMixolimnionASV$Sum <- rep(0, each = nrow(onlyMixolimnionASV))
      onlyMixolimnionASV$Count <- rep(0, each = nrow(onlyMixolimnionASV))
      for (i in onlyMixolimnionASV$RangInterest) { onlyMixolimnionASV$Count[which(onlyMixolimnionASV$RangInterest == i)] <- nrow(data_asv_tax %>% filter(Zone == "Mixolimnion") %>% filter(RangInterest == i))}
      onlyMixolimnionASV$Sum <- sum(onlyMixolimnionASV$Count)
      colnames(onlyMonimolimnionASV)[2]  <- "value"
      onlyMonimolimnionASV$Sum <- rep(0, each = nrow(onlyMonimolimnionASV))
      onlyMonimolimnionASV$Count <- rep(0, each = nrow(onlyMonimolimnionASV))
      for (i in onlyMonimolimnionASV$RangInterest) { onlyMonimolimnionASV$Count[which(onlyMonimolimnionASV$RangInterest == i)] <- nrow(data_asv_tax %>% filter(Zone == "Monimolimnion") %>% filter(RangInterest == i))}
      onlyMonimolimnionASV$Sum <- sum(onlyMonimolimnionASV$Count)
      onlyZoneASV <- rbind(onlyMixolimnionASV,onlyMonimolimnionASV)
    ## Plot
      bz <- ggplot(onlyZoneASV, mapping = aes(y= value, x = Zone, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
        geom_label(aes(y = 106,label = paste(Sum," ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") + scale_fill_manual(values = palette)
      bz <- bz + labs(x="Zones",y="ASVs (%)")
      print(bz)
#
    # Fraction -------------------------------------------------------------------
    ## Small
      onlySmallASV <- data_asv_tax %>% filter(Fraction == "Small") %>% select(ASV_Id,TotalSmall,RangInterest)
      row.names(onlySmallASV)<-onlySmallASV$ASV_Id ; onlySmallASV <- onlySmallASV %>% select(-ASV_Id)
      onlySmallASV <- onlySmallASV %>% group_by(RangInterest) %>% summarise_all(sum)
      onlySmallASV$TotalSmall <- onlySmallASV$TotalSmall*100/sum(onlySmallASV$TotalSmall)
      onlySmallASV <- onlySmallASV %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(TotalSmall) - 0.5*TotalSmall)
      onlySmallASV$label <- paste(round(onlySmallASV$TotalSmall,1), "%", sep = "")
      for (i in rownames(onlySmallASV)) {
        if (onlySmallASV[i,"label"] == "0%") { onlySmallASV[i,"label"] <- NA}}
      for (i in rownames(onlySmallASV)) {
        if (is.na(onlySmallASV[i,"label"]) == FALSE) { onlySmallASV[i,"label"] <- paste(onlySmallASV[i,"RangInterest"]," : ",onlySmallASV[i,"label"], sep = "")}}
      onlySmallASV$Fraction<- rep("Small", each = nrow(onlySmallASV))
    ## Large
      onlyLargeASV <- data_asv_tax %>% filter(Fraction == "Large") %>% select(ASV_Id,TotalLarge,RangInterest)
      row.names(onlyLargeASV)<-onlyLargeASV$ASV_Id
      onlyLargeASV <- onlyLargeASV %>% select(-ASV_Id)
      onlyLargeASV <- onlyLargeASV %>% group_by(RangInterest) %>% summarise_all(sum)
      onlyLargeASV$TotalLarge <- onlyLargeASV$TotalLarge*100/sum(onlyLargeASV$TotalLarge)
      onlyLargeASV <- onlyLargeASV %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(TotalLarge) - 0.5*TotalLarge)
      onlyLargeASV$label <- paste(round(onlyLargeASV$TotalLarge,1), "%", sep = "")
      for (i in rownames(onlyLargeASV)) {
        if (onlyLargeASV[i,"label"] == "0%") { onlyLargeASV[i,"label"] <- NA}}
      for (i in rownames(onlyLargeASV)) {
        if (is.na(onlyLargeASV[i,"label"]) == FALSE) { onlyLargeASV[i,"label"] <- paste(onlyLargeASV[i,"RangInterest"]," : ",onlyLargeASV[i,"label"], sep = "")}}
      onlyLargeASV$Fraction<- rep("Large", each = nrow(onlyLargeASV))
    ## Fraction
      colnames(onlySmallASV)[2]  <- "value"
      onlySmallASV$Sum <- rep(0, each = nrow(onlySmallASV))
      onlySmallASV$Count <- rep(0, each = nrow(onlySmallASV))
      for (i in onlySmallASV$RangInterest) { onlySmallASV$Count[which(onlySmallASV$RangInterest == i)] <- nrow(data_asv_tax %>% filter(Fraction == "Small") %>% filter(RangInterest == i))}
      onlySmallASV$Sum <- sum(onlySmallASV$Count)
      colnames(onlyLargeASV)[2]  <- "value"
      onlyLargeASV$Sum <- rep(0, each = nrow(onlyLargeASV))
      onlyLargeASV$Count <- rep(0, each = nrow(onlyLargeASV))
      for (i in onlyLargeASV$RangInterest) { onlyLargeASV$Count[which(onlyLargeASV$RangInterest == i)] <- nrow(data_asv_tax %>% filter(Fraction == "Large") %>% filter(RangInterest == i))}
      onlyLargeASV$Sum <- sum(onlyLargeASV$Count)
      onlyFractionASV <- rbind(onlySmallASV,onlyLargeASV)
    ## plot
      cz <- ggplot(onlyFractionASV, mapping = aes(y= value, x = Fraction, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
        geom_label(aes(y = 106,label = paste(Sum," ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
        scale_fill_manual(values = palette)
      cz <- cz + labs(x="Fractions",y="ASVs (%)")
      print(cz)
#
    # Coplot -------------------------------------------------------------------
      svglite("Hist-Taxomy/ASV-only.svg",width = 12.00,height = 6.00)
      b_plot <- plot_grid(az,bz,cz,  ncol = 3, nrow = 1, rel_widths = c(3,3,3),rel_heights = c(3))
      print(b_plot)
      dev.off()  
  # Séquence Total ---------------------------------------------------------------------
    # Cycle -------------------------------------------------------------------
    ## Day
      totalDaySequence <- data_seq_tax %>% select(ASV_Id,TotalDay,RangInterest)
      row.names(totalDaySequence)<-totalDaySequence$ASV_Id ; totalDaySequence <- totalDaySequence %>% select(-ASV_Id)
      totalDaySequence <- totalDaySequence %>% group_by(RangInterest) %>% summarise_all(sum)
      totalDaySequence$TotalDay <- totalDaySequence$TotalDay*100/sum(totalDaySequence$TotalDay)
      totalDaySequence <- totalDaySequence %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(TotalDay) - 0.5*TotalDay)
      totalDaySequence$label <- paste(round(totalDaySequence$TotalDay,1), "%", sep = "")
      for (i in rownames(totalDaySequence)) {
        if (totalDaySequence[i,"label"] == "0%") { totalDaySequence[i,"label"] <- NA}}
      for (i in rownames(totalDaySequence)) {
        if (is.na(totalDaySequence[i,"label"]) == FALSE) { totalDaySequence[i,"label"] <- paste(totalDaySequence[i,"RangInterest"]," : ",totalDaySequence[i,"label"], sep = "")}}
      totalDaySequence$Cycle<- rep("Day", each = nrow(totalDaySequence))
    ## Night
      totalNightSequence <- data_seq_tax %>% select(ASV_Id,TotalNight,RangInterest)
      row.names(totalNightSequence)<-totalNightSequence$ASV_Id ; totalNightSequence <- totalNightSequence %>% select(-ASV_Id)
      totalNightSequence <- totalNightSequence %>% group_by(RangInterest) %>% summarise_all(sum)
      totalNightSequence$TotalNight <- totalNightSequence$TotalNight*100/sum(totalNightSequence$TotalNight)
      totalNightSequence <- totalNightSequence %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(TotalNight) - 0.5*TotalNight)
      totalNightSequence$label <- paste(round(totalNightSequence$TotalNight,1), "%", sep = "")
      for (i in rownames(totalNightSequence)) {
        if (totalNightSequence[i,"label"] == "0%") { totalNightSequence[i,"label"] <- NA}}
      for (i in rownames(totalNightSequence)) {
        if (is.na(totalNightSequence[i,"label"]) == FALSE) { totalNightSequence[i,"label"] <- paste(totalNightSequence[i,"RangInterest"]," : ",totalNightSequence[i,"label"], sep = "")}}
      totalNightSequence$Cycle<- rep("Night", each = nrow(totalNightSequence))
    ## Cycle
      colnames(totalDaySequence)[2]  <- "value"
      totalDaySequence$Sum <- rep(0, each = nrow(totalDaySequence))
      totalDaySequence$Count <- rep(0, each = nrow(totalDaySequence))
      for (i in totalDaySequence$RangInterest) { totalDaySequence$Count[which(totalDaySequence$RangInterest == i)] <- sum(data_seq_tax  %>% filter(RangInterest == i) %>% select(TotalDay))}
      totalDaySequence$Sum <- sum(totalDaySequence$Count)
      colnames(totalNightSequence)[2]  <- "value"
      totalNightSequence$Sum <- rep(0, each = nrow(totalNightSequence))
      totalNightSequence$Count <- rep(0, each = nrow(totalNightSequence))
      for (i in totalNightSequence$RangInterest) { totalNightSequence$Count[which(totalNightSequence$RangInterest == i)] <- sum(data_seq_tax  %>% filter(RangInterest == i) %>% select(TotalNight))}
      totalNightSequence$Sum <- sum(totalNightSequence$Count)
      totalCycleSequence <- rbind(totalDaySequence,totalNightSequence)
      totalCycleSequence$percent <- paste("(",round(totalCycleSequence$Sum*100/(colSums(statRarefy %>% select(`apRarefy-Sequence`))/2),1)," %)",sep ="")
    ## Plot
      iy <- ggplot(totalCycleSequence, mapping = aes(y= value, x = Cycle, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
        scale_fill_manual(values = palette) + 
        geom_label(aes(y = 108,label = paste(Sum,"séquences",percent,sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
      legendSequence <- get_legend(iy)
      iy <- iy + labs(x="Cycles",y="Sequences (%)") + theme(legend.position = "none")
      print(iy)
#
    # Zone -------------------------------------------------------------------
    ## Mixolimnion
      totalMixolimnionSequence <- data_seq_tax %>% select(ASV_Id,TotalMixolimnion,RangInterest)
      row.names(totalMixolimnionSequence)<-totalMixolimnionSequence$ASV_Id ; totalMixolimnionSequence <- totalMixolimnionSequence %>% select(-ASV_Id)
      totalMixolimnionSequence <- totalMixolimnionSequence %>% group_by(RangInterest) %>% summarise_all(sum)
      totalMixolimnionSequence$TotalMixolimnion <- totalMixolimnionSequence$TotalMixolimnion*100/sum(totalMixolimnionSequence$TotalMixolimnion)
      totalMixolimnionSequence <- totalMixolimnionSequence %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(TotalMixolimnion) - 0.5*TotalMixolimnion)
      totalMixolimnionSequence$label <- paste(round(totalMixolimnionSequence$TotalMixolimnion,1), "%", sep = "")
      for (i in rownames(totalMixolimnionSequence)) {
        if (totalMixolimnionSequence[i,"label"] == "0%") { totalMixolimnionSequence[i,"label"] <- NA}}
      for (i in rownames(totalMixolimnionSequence)) {
        if (is.na(totalMixolimnionSequence[i,"label"]) == FALSE) { totalMixolimnionSequence[i,"label"] <- paste(totalMixolimnionSequence[i,"RangInterest"]," : ",totalMixolimnionSequence[i,"label"], sep = "")}}
      totalMixolimnionSequence$Zone<- rep("Mixolimnion", each = nrow(totalMixolimnionSequence))
    ## Monimolimnion
      totalMonimolimnionSequence <- data_seq_tax %>% select(ASV_Id,TotalMonimolimnion,RangInterest)
      row.names(totalMonimolimnionSequence)<-totalMonimolimnionSequence$ASV_Id ; totalMonimolimnionSequence <- totalMonimolimnionSequence %>% select(-ASV_Id)
      totalMonimolimnionSequence <- totalMonimolimnionSequence %>% group_by(RangInterest) %>% summarise_all(sum)
      totalMonimolimnionSequence$TotalMonimolimnion <- totalMonimolimnionSequence$TotalMonimolimnion*100/sum(totalMonimolimnionSequence$TotalMonimolimnion)
      totalMonimolimnionSequence <- totalMonimolimnionSequence %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(TotalMonimolimnion) - 0.5*TotalMonimolimnion)
      totalMonimolimnionSequence$label <- paste(round(totalMonimolimnionSequence$TotalMonimolimnion,1), "%", sep = "")
      for (i in rownames(totalMonimolimnionSequence)) {
        if (totalMonimolimnionSequence[i,"label"] == "0%") { totalMonimolimnionSequence[i,"label"] <- NA}}
      for (i in rownames(totalMonimolimnionSequence)) {
        if (is.na(totalMonimolimnionSequence[i,"label"]) == FALSE) { totalMonimolimnionSequence[i,"label"] <- paste(totalMonimolimnionSequence[i,"RangInterest"]," : ",totalMonimolimnionSequence[i,"label"], sep = "")}}
      totalMonimolimnionSequence$Zone<- rep("Monimolimnion", each = nrow(totalMonimolimnionSequence))
    ## Zone
      colnames(totalMixolimnionSequence)[2]  <- "value"
      totalMixolimnionSequence$Sum <- rep(0, each = nrow(totalMixolimnionSequence))
      totalMixolimnionSequence$Count <- rep(0, each = nrow(totalMixolimnionSequence))
      for (i in totalMixolimnionSequence$RangInterest) { totalMixolimnionSequence$Count[which(totalMixolimnionSequence$RangInterest == i)] <- sum(data_seq_tax  %>% filter(RangInterest == i) %>% select(TotalMixolimnion))}
      totalMixolimnionSequence$Sum <- sum(totalMixolimnionSequence$Count)
      colnames(totalMonimolimnionSequence)[2]  <- "value"
      totalMonimolimnionSequence$Sum <- rep(0, each = nrow(totalMonimolimnionSequence))
      totalMonimolimnionSequence$Count <- rep(0, each = nrow(totalMonimolimnionSequence))
      for (i in totalMonimolimnionSequence$RangInterest) { totalMonimolimnionSequence$Count[which(totalMonimolimnionSequence$RangInterest == i)] <- sum(data_seq_tax  %>% filter(RangInterest == i) %>% select(TotalMonimolimnion))}
      totalMonimolimnionSequence$Sum <- sum(totalMonimolimnionSequence$Count)
      totalZoneSequence <- rbind(totalMixolimnionSequence,totalMonimolimnionSequence)
      totalZoneSequence$percent <- paste("(",round(totalZoneSequence$Sum*100/(colSums(statRarefy %>% select(`apRarefy-Sequence`))/2),1)," %)",sep ="")
    ## Plot
      jy <- ggplot(totalZoneSequence, mapping = aes(y= value, x = Zone, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
        scale_fill_manual(values = palette) + 
        geom_label(aes(y = 108,label = paste(Sum,"séquences",percent,sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
      jy <- jy + labs(x="Zones",y="Sequences (%)") + theme(legend.position = "none")
      print(jy)
#
    # Fraction -------------------------------------------------------------------
    ## Small
      totalSmallSequence <- data_seq_tax %>% select(ASV_Id,TotalSmall,RangInterest)
      row.names(totalSmallSequence)<-totalSmallSequence$ASV_Id ; totalSmallSequence <- totalSmallSequence %>% select(-ASV_Id)
      totalSmallSequence <- totalSmallSequence %>% group_by(RangInterest) %>% summarise_all(sum)
      totalSmallSequence$TotalSmall <- totalSmallSequence$TotalSmall*100/sum(totalSmallSequence$TotalSmall)
      totalSmallSequence <- totalSmallSequence %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(TotalSmall) - 0.5*TotalSmall)
      totalSmallSequence$label <- paste(round(totalSmallSequence$TotalSmall,1), "%", sep = "")
      for (i in rownames(totalSmallSequence)) {
        if (totalSmallSequence[i,"label"] == "0%") { totalSmallSequence[i,"label"] <- NA}}
      for (i in rownames(totalSmallSequence)) {
        if (is.na(totalSmallSequence[i,"label"]) == FALSE) { totalSmallSequence[i,"label"] <- paste(totalSmallSequence[i,"RangInterest"]," : ",totalSmallSequence[i,"label"], sep = "")}}
      totalSmallSequence$Fraction<- rep("Small", each = nrow(totalSmallSequence))
    ## Large
      totalLargeSequence <- data_seq_tax %>% select(ASV_Id,TotalLarge,RangInterest)
      row.names(totalLargeSequence)<-totalLargeSequence$ASV_Id ; totalLargeSequence <- totalLargeSequence %>% select(-ASV_Id)
      totalLargeSequence <- totalLargeSequence %>% group_by(RangInterest) %>% summarise_all(sum)
      totalLargeSequence$TotalLarge <- totalLargeSequence$TotalLarge*100/sum(totalLargeSequence$TotalLarge)
      totalLargeSequence <- totalLargeSequence %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(TotalLarge) - 0.5*TotalLarge)
      totalLargeSequence$label <- paste(round(totalLargeSequence$TotalLarge,1), "%", sep = "")
      for (i in rownames(totalLargeSequence)) {
        if (totalLargeSequence[i,"label"] == "0%") { totalLargeSequence[i,"label"] <- NA}}
      for (i in rownames(totalLargeSequence)) {
        if (is.na(totalLargeSequence[i,"label"]) == FALSE) { totalLargeSequence[i,"label"] <- paste(totalLargeSequence[i,"RangInterest"]," : ",totalLargeSequence[i,"label"], sep = "")}}
      totalLargeSequence$Fraction<- rep("Large", each = nrow(totalLargeSequence))
    ## Fraction
      colnames(totalSmallSequence)[2]  <- "value"
      totalSmallSequence$Sum <- rep(0, each = nrow(totalSmallSequence))
      totalSmallSequence$Count <- rep(0, each = nrow(totalSmallSequence))
      for (i in totalSmallSequence$RangInterest) { totalSmallSequence$Count[which(totalSmallSequence$RangInterest == i)] <- sum(data_seq_tax  %>% filter(RangInterest == i) %>% select(TotalSmall))}
      totalSmallSequence$Sum <- sum(totalSmallSequence$Count)
      colnames(totalLargeSequence)[2]  <- "value"
      totalLargeSequence$Sum <- rep(0, each = nrow(totalLargeSequence))
      totalLargeSequence$Count <- rep(0, each = nrow(totalLargeSequence))
      for (i in totalLargeSequence$RangInterest) { totalLargeSequence$Count[which(totalLargeSequence$RangInterest == i)] <- sum(data_seq_tax  %>% filter(RangInterest == i) %>% select(TotalLarge))}
      totalLargeSequence$Sum <- sum(totalLargeSequence$Count)
      totalFractionSequence <- rbind(totalSmallSequence,totalLargeSequence)
      totalFractionSequence$percent <- paste("(",round(totalFractionSequence$Sum*100/(colSums(statRarefy %>% select(`apRarefy-Sequence`))/2),1)," %)",sep ="")
    ## Plot
      ky <- ggplot(totalFractionSequence, mapping = aes(y= value, x = Fraction, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
        scale_fill_manual(values = palette) + 
        geom_label(aes(y = 108,label = paste(Sum,"séquences",percent,sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
      ky <- ky + labs(x="Fractions",y="Sequences (%)") + theme(legend.position = "none")
      print(ky)
#
    # Period -------------------------------------------------------------------
    ## Avril
      total04Sequence <- data_seq_tax %>% select(ASV_Id,Total04,RangInterest)
      row.names(total04Sequence)<-total04Sequence$ASV_Id ; total04Sequence <- total04Sequence %>% select(-ASV_Id)
      total04Sequence <- total04Sequence %>% group_by(RangInterest) %>% summarise_all(sum)
      total04Sequence$Total04 <- total04Sequence$Total04*100/sum(total04Sequence$Total04)
      total04Sequence <- total04Sequence %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(Total04) - 0.5*Total04)
      total04Sequence$label <- paste(round(total04Sequence$Total04,1), "%", sep = "")
      for (i in rownames(total04Sequence)) {
        if (total04Sequence[i,"label"] == "0%") { total04Sequence[i,"label"] <- NA}}
      for (i in rownames(total04Sequence)) {
        if (is.na(total04Sequence[i,"label"]) == FALSE) { total04Sequence[i,"label"] <- paste(total04Sequence[i,"RangInterest"]," : ",total04Sequence[i,"label"], sep = "")}}
      total04Sequence$Period<- rep("04", each = nrow(total04Sequence))
    ## Juin
      total06Sequence <- data_seq_tax %>% select(ASV_Id,Total06,RangInterest)
      row.names(total06Sequence)<-total06Sequence$ASV_Id ; total06Sequence <- total06Sequence %>% select(-ASV_Id)
      total06Sequence <- total06Sequence %>% group_by(RangInterest) %>% summarise_all(sum)
      total06Sequence$Total06 <- total06Sequence$Total06*100/sum(total06Sequence$Total06)
      total06Sequence <- total06Sequence %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(Total06) - 0.5*Total06)
      total06Sequence$label <- paste(round(total06Sequence$Total06,1), "%", sep = "")
      for (i in rownames(total06Sequence)) {
        if (total06Sequence[i,"label"] == "0%") { total06Sequence[i,"label"] <- NA}}
      for (i in rownames(total06Sequence)) {
        if (is.na(total06Sequence[i,"label"]) == FALSE) { total06Sequence[i,"label"] <- paste(total06Sequence[i,"RangInterest"]," : ",total06Sequence[i,"label"], sep = "")}}
      total06Sequence$Period<- rep("06", each = nrow(total06Sequence))
    ## Septembre
      total09Sequence <- data_seq_tax %>% select(ASV_Id,Total09,RangInterest)
      row.names(total09Sequence)<-total09Sequence$ASV_Id ; total09Sequence <- total09Sequence %>% select(-ASV_Id)
      total09Sequence <- total09Sequence %>% group_by(RangInterest) %>% summarise_all(sum)
      total09Sequence$Total09 <- total09Sequence$Total09*100/sum(total09Sequence$Total09)
      total09Sequence <- total09Sequence %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(Total09) - 0.5*Total09)
      total09Sequence$label <- paste(round(total09Sequence$Total09,1), "%", sep = "")
      for (i in rownames(total09Sequence)) {
        if (total09Sequence[i,"label"] == "0%") { total09Sequence[i,"label"] <- NA}}
      for (i in rownames(total09Sequence)) {
        if (is.na(total09Sequence[i,"label"]) == FALSE) { total09Sequence[i,"label"] <- paste(total09Sequence[i,"RangInterest"]," : ",total09Sequence[i,"label"], sep = "")}}
      total09Sequence$Period<- rep("09", each = nrow(total09Sequence))
    ## Octobre
      total11Sequence <- data_seq_tax %>% select(ASV_Id,Total11,RangInterest)
      row.names(total11Sequence)<-total11Sequence$ASV_Id ; total11Sequence <- total11Sequence %>% select(-ASV_Id)
      total11Sequence <- total11Sequence %>% group_by(RangInterest) %>% summarise_all(sum)
      total11Sequence$Total11 <- total11Sequence$Total11*100/sum(total11Sequence$Total11)
      total11Sequence <- total11Sequence %>%
        arrange(desc(RangInterest)) %>%
        mutate(lab.ypos = cumsum(Total11) - 0.5*Total11)
      total11Sequence$label <- paste(round(total11Sequence$Total11,1), "%", sep = "")
      for (i in rownames(total11Sequence)) {
        if (total11Sequence[i,"label"] == "0%") { total11Sequence[i,"label"] <- NA}}
      for (i in rownames(total11Sequence)) {
        if (is.na(total11Sequence[i,"label"]) == FALSE) { total11Sequence[i,"label"] <- paste(total11Sequence[i,"RangInterest"]," : ",total11Sequence[i,"label"], sep = "")}}
      total11Sequence$Period<- rep("11", each = nrow(total11Sequence))
    ## Periods
      colnames(total04Sequence)[2]  <- "value"
      total04Sequence$Sum <- rep(0, each = nrow(total04Sequence))
      total04Sequence$Count <- rep(0, each = nrow(total04Sequence))
      for (i in total04Sequence$RangInterest) { total04Sequence$Count[which(total04Sequence$RangInterest == i)] <- sum(data_seq_tax  %>% filter(RangInterest == i) %>% select(Total04))}
      total04Sequence$Sum <- sum(total04Sequence$Count)
      colnames(total06Sequence)[2]  <- "value"
      total06Sequence$Sum <- rep(0, each = nrow(total06Sequence))
      total06Sequence$Count <- rep(0, each = nrow(total06Sequence))
      for (i in total06Sequence$RangInterest) { total06Sequence$Count[which(total06Sequence$RangInterest == i)] <- sum(data_seq_tax  %>% filter(RangInterest == i) %>% select(Total06))}
      total06Sequence$Sum <- sum(total06Sequence$Count)
      colnames(total09Sequence)[2]  <- "value"
      total09Sequence$Sum <- rep(0, each = nrow(total09Sequence))
      total09Sequence$Count <- rep(0, each = nrow(total09Sequence))
      for (i in total09Sequence$RangInterest) { total09Sequence$Count[which(total09Sequence$RangInterest == i)] <- sum(data_seq_tax  %>% filter(RangInterest == i) %>% select(Total09))}
      total09Sequence$Sum <- sum(total09Sequence$Count)
      colnames(total11Sequence)[2]  <- "value"
      total11Sequence$Sum <- rep(0, each = nrow(total11Sequence))
      total11Sequence$Count <- rep(0, each = nrow(total11Sequence))
      for (i in total11Sequence$RangInterest) { total11Sequence$Count[which(total11Sequence$RangInterest == i)] <- sum(data_seq_tax  %>% filter(RangInterest == i) %>% select(Total11))}
      total11Sequence$Sum <- sum(total11Sequence$Count)
      totalPeriodSequence <- rbind(total04Sequence,total06Sequence,total09Sequence,total11Sequence)
      totalPeriodSequence$percent <- paste("(",round(totalPeriodSequence$Sum*100/(colSums(statRarefy %>% select(`apRarefy-Sequence`))/4),1)," %)",sep ="")
    ## Plot
      ly <- ggplot(totalPeriodSequence, mapping = aes(y= value, x = Period, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
        scale_fill_manual(values = palette) + 
        geom_label(aes(y = 108,label = paste(Sum,"séquences",percent,sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
      ly <- ly + labs(x="Periods",y="Sequences (%)") + theme(legend.position = "none")
      print(ly)
#
    # Coplot -------------------------------------------------------------------
      svglite("Hist-Taxomy/Sequence-Total.svg",width = 12.00,height = 7.00)
      b_plot <- plot_grid(iy,jy,ky,ly, legendSequence, ncol = 5, nrow = 1, rel_widths = c(3,3,3,5,2.5),rel_heights = c(3))
      print(b_plot)
      dev.off()
  # Séquence Only ---------------------------------------------------------------------
    # 100% - 0% ---------------------------------------------------------------
      # Cycle -------------------------------------------------------------------
      ## Day
        onlyDaySequence <- data_seq_tax %>% filter(Cycle == "Day") %>% select(ASV_Id,TotalDay,RangInterest)
        row.names(onlyDaySequence)<-onlyDaySequence$ASV_Id ; onlyDaySequence <- onlyDaySequence %>% select(-ASV_Id)
        onlyDaySequence <- onlyDaySequence %>% group_by(RangInterest) %>% summarise_all(sum)
        onlyDaySequence$TotalDay <- onlyDaySequence$TotalDay*100/sum(onlyDaySequence$TotalDay)
        onlyDaySequence <- onlyDaySequence %>%
          arrange(desc(RangInterest)) %>%
          mutate(lab.ypos = cumsum(TotalDay) - 0.5*TotalDay)
        onlyDaySequence$label <- paste(round(onlyDaySequence$TotalDay,1), "%", sep = "")
        for (i in rownames(onlyDaySequence)) {
          if (onlyDaySequence[i,"label"] == "0%") { onlyDaySequence[i,"label"] <- NA}}
        for (i in rownames(onlyDaySequence)) {
          if (is.na(onlyDaySequence[i,"label"]) == FALSE) { onlyDaySequence[i,"label"] <- paste(onlyDaySequence[i,"RangInterest"]," : ",onlyDaySequence[i,"label"], sep = "")}}
        onlyDaySequence$Cycle<- rep("Day", each = nrow(onlyDaySequence))
      ## Night
        onlyNightSequence <- data_seq_tax %>% filter(Cycle == "Night") %>% select(ASV_Id,TotalNight,RangInterest)
        row.names(onlyNightSequence)<-onlyNightSequence$ASV_Id ; onlyNightSequence <- onlyNightSequence %>% select(-ASV_Id)
        onlyNightSequence <- onlyNightSequence %>% group_by(RangInterest) %>% summarise_all(sum)
        onlyNightSequence$TotalNight <- onlyNightSequence$TotalNight*100/sum(onlyNightSequence$TotalNight)
        onlyNightSequence <- onlyNightSequence %>%
          arrange(desc(RangInterest)) %>%
          mutate(lab.ypos = cumsum(TotalNight) - 0.5*TotalNight)
        onlyNightSequence$label <- paste(round(onlyNightSequence$TotalNight,1), "%", sep = "")
        for (i in rownames(onlyNightSequence)) {
          if (onlyNightSequence[i,"label"] == "0%") { onlyNightSequence[i,"label"] <- NA}}
        for (i in rownames(onlyNightSequence)) {
          if (is.na(onlyNightSequence[i,"label"]) == FALSE) { onlyNightSequence[i,"label"] <- paste(onlyNightSequence[i,"RangInterest"]," : ",onlyNightSequence[i,"label"], sep = "")}}
        onlyNightSequence$Cycle<- rep("Night", each = nrow(onlyNightSequence))
      ## Cycle
        colnames(onlyDaySequence)[2]  <- "value"
        onlyDaySequence$Sum <- rep(0, each = nrow(onlyDaySequence))
        onlyDaySequence$Count <- rep(0, each = nrow(onlyDaySequence))
        for (i in onlyDaySequence$RangInterest) { onlyDaySequence$Count[which(onlyDaySequence$RangInterest == i)] <- sum(data_seq_tax %>% filter(Cycle == "Day") %>% filter(RangInterest == i) %>% select(TotalDay))}
        onlyDaySequence$Sum <- sum(onlyDaySequence$Count)
        colnames(onlyNightSequence)[2]  <- "value"
        onlyNightSequence$Sum <- rep(0, each = nrow(onlyNightSequence))
        onlyNightSequence$Count <- rep(0, each = nrow(onlyNightSequence))
        for (i in onlyNightSequence$RangInterest) { onlyNightSequence$Count[which(onlyNightSequence$RangInterest == i)] <- sum(data_seq_tax %>% filter(Cycle == "Night") %>% filter(RangInterest == i) %>% select(TotalNight))}
        onlyNightSequence$Sum <- sum(onlyNightSequence$Count)
        onlyCycleSequence <- rbind(onlyDaySequence,onlyNightSequence)
      ## Plot
        iz <- ggplot(onlyCycleSequence, mapping = aes(y= value, x = Cycle, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
          geom_label(aes(y = 106,label = paste(Sum,"séquences",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
          scale_fill_manual(values = palette)
        iz <- iz + labs(x="Cycles",y="Sequences (%)")
        print(iz)
#
      # Zone -------------------------------------------------------------------
      ## Mixolimnion
        onlyMixolimnionSequence <- data_seq_tax %>% filter(Zone == "Mixolimnion") %>% select(ASV_Id,TotalMixolimnion,RangInterest)
        row.names(onlyMixolimnionSequence)<-onlyMixolimnionSequence$ASV_Id ; onlyMixolimnionSequence <- onlyMixolimnionSequence %>% select(-ASV_Id)
        onlyMixolimnionSequence <- onlyMixolimnionSequence %>% group_by(RangInterest) %>% summarise_all(sum)
        onlyMixolimnionSequence$TotalMixolimnion <- onlyMixolimnionSequence$TotalMixolimnion*100/sum(onlyMixolimnionSequence$TotalMixolimnion)
        onlyMixolimnionSequence <- onlyMixolimnionSequence %>%
          arrange(desc(RangInterest)) %>%
          mutate(lab.ypos = cumsum(TotalMixolimnion) - 0.5*TotalMixolimnion)
        onlyMixolimnionSequence$label <- paste(round(onlyMixolimnionSequence$TotalMixolimnion,1), "%", sep = "")
        for (i in rownames(onlyMixolimnionSequence)) {
          if (onlyMixolimnionSequence[i,"label"] == "0%") { onlyMixolimnionSequence[i,"label"] <- NA}}
        for (i in rownames(onlyMixolimnionSequence)) {
          if (is.na(onlyMixolimnionSequence[i,"label"]) == FALSE) { onlyMixolimnionSequence[i,"label"] <- paste(onlyMixolimnionSequence[i,"RangInterest"]," : ",onlyMixolimnionSequence[i,"label"], sep = "")}}
        onlyMixolimnionSequence$Zone<- rep("Mixolimnion", each = nrow(onlyMixolimnionSequence))
      ## Monimolimnion
        onlyMonimolimnionSequence <- data_seq_tax %>% filter(Zone == "Monimolimnion") %>% select(ASV_Id,TotalMonimolimnion,RangInterest)
        row.names(onlyMonimolimnionSequence)<-onlyMonimolimnionSequence$ASV_Id ; onlyMonimolimnionSequence <- onlyMonimolimnionSequence %>% select(-ASV_Id)
        onlyMonimolimnionSequence <- onlyMonimolimnionSequence %>% group_by(RangInterest) %>% summarise_all(sum)
        onlyMonimolimnionSequence$TotalMonimolimnion <- onlyMonimolimnionSequence$TotalMonimolimnion*100/sum(onlyMonimolimnionSequence$TotalMonimolimnion)
        onlyMonimolimnionSequence <- onlyMonimolimnionSequence %>%
          arrange(desc(RangInterest)) %>%
          mutate(lab.ypos = cumsum(TotalMonimolimnion) - 0.5*TotalMonimolimnion)
        onlyMonimolimnionSequence$label <- paste(round(onlyMonimolimnionSequence$TotalMonimolimnion,1), "%", sep = "")
        for (i in rownames(onlyMonimolimnionSequence)) {
          if (onlyMonimolimnionSequence[i,"label"] == "0%") { onlyMonimolimnionSequence[i,"label"] <- NA}}
        for (i in rownames(onlyMonimolimnionSequence)) {
          if (is.na(onlyMonimolimnionSequence[i,"label"]) == FALSE) { onlyMonimolimnionSequence[i,"label"] <- paste(onlyMonimolimnionSequence[i,"RangInterest"]," : ",onlyMonimolimnionSequence[i,"label"], sep = "")}}
        onlyMonimolimnionSequence$Zone<- rep("Monimolimnion", each = nrow(onlyMonimolimnionSequence))
      ## Zone
        colnames(onlyMixolimnionSequence)[2]  <- "value"
        onlyMixolimnionSequence$Sum <- rep(0, each = nrow(onlyMixolimnionSequence))
        onlyMixolimnionSequence$Count <- rep(0, each = nrow(onlyMixolimnionSequence))
        for (i in onlyMixolimnionSequence$RangInterest) { onlyMixolimnionSequence$Count[which(onlyMixolimnionSequence$RangInterest == i)] <- sum(data_seq_tax %>% filter(Zone == "Mixolimnion") %>% filter(RangInterest == i) %>% select(TotalMixolimnion))}
        onlyMixolimnionSequence$Sum <- sum(onlyMixolimnionSequence$Count)
        colnames(onlyMonimolimnionSequence)[2]  <- "value"
        onlyMonimolimnionSequence$Sum <- rep(0, each = nrow(onlyMonimolimnionSequence))
        onlyMonimolimnionSequence$Count <- rep(0, each = nrow(onlyMonimolimnionSequence))
        for (i in onlyMonimolimnionSequence$RangInterest) { onlyMonimolimnionSequence$Count[which(onlyMonimolimnionSequence$RangInterest == i)] <- sum(data_seq_tax %>% filter(Zone == "Monimolimnion") %>% filter(RangInterest == i) %>% select(TotalMonimolimnion))}
        onlyMonimolimnionSequence$Sum <- sum(onlyMonimolimnionSequence$Count)
        onlyZoneSequence <- rbind(onlyMixolimnionSequence,onlyMonimolimnionSequence)
      ## Plot
        jz <- ggplot(onlyZoneSequence, mapping = aes(y= value, x = Zone, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
          geom_label(aes(y = 106,label = paste(Sum,"séquences",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
          scale_fill_manual(values = palette)
        jz <- jz + labs(x="Zones",y="Sequences (%)") #+ theme(legend.position = "none") 
        print(jz)
#
      # Fraction -------------------------------------------------------------------
      ## Small
        onlySmallSequence <- data_seq_tax %>% filter(Fraction == "Small") %>% select(ASV_Id,TotalSmall,RangInterest)
        row.names(onlySmallSequence)<-onlySmallSequence$ASV_Id ; onlySmallSequence <- onlySmallSequence %>% select(-ASV_Id)
        onlySmallSequence <- onlySmallSequence %>% group_by(RangInterest) %>% summarise_all(sum)
        onlySmallSequence$TotalSmall <- onlySmallSequence$TotalSmall*100/sum(onlySmallSequence$TotalSmall)
        onlySmallSequence <- onlySmallSequence %>%
          arrange(desc(RangInterest)) %>%
          mutate(lab.ypos = cumsum(TotalSmall) - 0.5*TotalSmall)
        onlySmallSequence$label <- paste(round(onlySmallSequence$TotalSmall,1), "%", sep = "")
        for (i in rownames(onlySmallSequence)) {
          if (onlySmallSequence[i,"label"] == "0%") { onlySmallSequence[i,"label"] <- NA}}
        for (i in rownames(onlySmallSequence)) {
          if (is.na(onlySmallSequence[i,"label"]) == FALSE) { onlySmallSequence[i,"label"] <- paste(onlySmallSequence[i,"RangInterest"]," : ",onlySmallSequence[i,"label"], sep = "")}}
        onlySmallSequence$Fraction<- rep("Small", each = nrow(onlySmallSequence))
      ## Large
        onlyLargeSequence <- data_seq_tax %>% filter(Fraction == "Large") %>% select(ASV_Id,TotalLarge,RangInterest)
        row.names(onlyLargeSequence)<-onlyLargeSequence$ASV_Id ; onlyLargeSequence <- onlyLargeSequence %>% select(-ASV_Id)
        onlyLargeSequence <- onlyLargeSequence %>% group_by(RangInterest) %>% summarise_all(sum)
        onlyLargeSequence$TotalLarge <- onlyLargeSequence$TotalLarge*100/sum(onlyLargeSequence$TotalLarge)
        onlyLargeSequence <- onlyLargeSequence %>%
          arrange(desc(RangInterest)) %>%
          mutate(lab.ypos = cumsum(TotalLarge) - 0.5*TotalLarge)
        onlyLargeSequence$label <- paste(round(onlyLargeSequence$TotalLarge,1), "%", sep = "")
        for (i in rownames(onlyLargeSequence)) {
          if (onlyLargeSequence[i,"label"] == "0%") { onlyLargeSequence[i,"label"] <- NA}}
        for (i in rownames(onlyLargeSequence)) {
          if (is.na(onlyLargeSequence[i,"label"]) == FALSE) { onlyLargeSequence[i,"label"] <- paste(onlyLargeSequence[i,"RangInterest"]," : ",onlyLargeSequence[i,"label"], sep = "")}}
        onlyLargeSequence$Fraction<- rep("Large", each = nrow(onlyLargeSequence))
      ## Fraction
        colnames(onlySmallSequence)[2]  <- "value"
        onlySmallSequence$Sum <- rep(0, each = nrow(onlySmallSequence))
        onlySmallSequence$Count <- rep(0, each = nrow(onlySmallSequence))
        for (i in onlySmallSequence$RangInterest) { onlySmallSequence$Count[which(onlySmallSequence$RangInterest == i)] <- sum(data_seq_tax %>% filter(Fraction == "Small") %>% filter(RangInterest == i) %>% select(TotalSmall))}
        onlySmallSequence$Sum <- sum(onlySmallSequence$Count)
        colnames(onlyLargeSequence)[2]  <- "value"
        onlyLargeSequence$Sum <- rep(0, each = nrow(onlyLargeSequence))
        onlyLargeSequence$Count <- rep(0, each = nrow(onlyLargeSequence))
        for (i in onlyLargeSequence$RangInterest) { onlyLargeSequence$Count[which(onlyLargeSequence$RangInterest == i)] <- sum(data_seq_tax %>% filter(Fraction == "Large") %>% filter(RangInterest == i) %>% select(TotalLarge))}
        onlyLargeSequence$Sum <- sum(onlyLargeSequence$Count)
        onlyFractionSequence <- rbind(onlySmallSequence,onlyLargeSequence)
      ## Plot
        kz <- ggplot(onlyFractionSequence, mapping = aes(y= value, x = Fraction, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
          geom_label(aes(y = 106,label = paste(Sum,"séquences",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
          scale_fill_manual(values = palette)
        kz <- kz + labs(x="Fractions",y="Sequences (%)") #+ theme(legend.position = "none")
        print(kz)
#
      # Coplot -------------------------------------------------------------------
        svglite("Hist-Taxomy/Sequence-only-100.svg",width = 12.00,height = 6.00)
        b_plot <- plot_grid(iz,jz,kz,  ncol = 3, nrow = 1, rel_widths = c(3,3,3),rel_heights = c(3))
        print(b_plot)
        dev.off()  
#
    # 90% - 10% ---------------------------------------------------------------
      # Cycle -------------------------------------------------------------------
      ## Day
        onlyDaySequence90 <- data_seq_tax %>% filter(Cycle90 == "Day") %>% select(ASV_Id,TotalDay,RangInterest)
        row.names(onlyDaySequence90)<-onlyDaySequence90$ASV_Id ; onlyDaySequence90 <- onlyDaySequence90 %>% select(-ASV_Id)
        onlyDaySequence90 <- onlyDaySequence90 %>% group_by(RangInterest) %>% summarise_all(sum)
        onlyDaySequence90$TotalDay <- onlyDaySequence90$TotalDay*100/sum(onlyDaySequence90$TotalDay)
        onlyDaySequence90 <- onlyDaySequence90 %>%
          arrange(desc(RangInterest)) %>%
          mutate(lab.ypos = cumsum(TotalDay) - 0.5*TotalDay)
        onlyDaySequence90$label <- paste(round(onlyDaySequence90$TotalDay,1), "%", sep = "")
        for (i in rownames(onlyDaySequence90)) {
          if (onlyDaySequence90[i,"label"] == "0%") { onlyDaySequence90[i,"label"] <- NA}}
        for (i in rownames(onlyDaySequence90)) {
          if (is.na(onlyDaySequence90[i,"label"]) == FALSE) { onlyDaySequence90[i,"label"] <- paste(onlyDaySequence90[i,"RangInterest"]," : ",onlyDaySequence90[i,"label"], sep = "")}}
        onlyDaySequence90$Cycle90<- rep("Day", each = nrow(onlyDaySequence90))
      ## Night
        onlyNightSequence90 <- data_seq_tax %>% filter(Cycle90 == "Night") %>% select(ASV_Id,TotalNight,RangInterest)
        row.names(onlyNightSequence90)<-onlyNightSequence90$ASV_Id ; onlyNightSequence90 <- onlyNightSequence90 %>% select(-ASV_Id)
        onlyNightSequence90 <- onlyNightSequence90 %>% group_by(RangInterest) %>% summarise_all(sum)
        onlyNightSequence90$TotalNight <- onlyNightSequence90$TotalNight*100/sum(onlyNightSequence90$TotalNight)
        onlyNightSequence90 <- onlyNightSequence90 %>%
          arrange(desc(RangInterest)) %>%
          mutate(lab.ypos = cumsum(TotalNight) - 0.5*TotalNight)
        onlyNightSequence90$label <- paste(round(onlyNightSequence90$TotalNight,1), "%", sep = "")
        for (i in rownames(onlyNightSequence90)) {
          if (onlyNightSequence90[i,"label"] == "0%") { onlyNightSequence90[i,"label"] <- NA}}
        for (i in rownames(onlyNightSequence90)) {
          if (is.na(onlyNightSequence90[i,"label"]) == FALSE) { onlyNightSequence90[i,"label"] <- paste(onlyNightSequence90[i,"RangInterest"]," : ",onlyNightSequence90[i,"label"], sep = "")}}
        onlyNightSequence90$Cycle90<- rep("Night", each = nrow(onlyNightSequence90))
      ## Cycle
        colnames(onlyDaySequence90)[2]  <- "value"
        onlyDaySequence90$Sum <- rep(0, each = nrow(onlyDaySequence90))
        onlyDaySequence90$Count <- rep(0, each = nrow(onlyDaySequence90))
        for (i in onlyDaySequence90$RangInterest) { onlyDaySequence90$Count[which(onlyDaySequence90$RangInterest == i)] <- sum(data_seq_tax %>% filter(Cycle90 == "Day") %>% filter(RangInterest == i) %>% select(TotalDay))}
        onlyDaySequence90$Sum <- sum(onlyDaySequence90$Count)
        colnames(onlyNightSequence90)[2]  <- "value"
        onlyNightSequence90$Sum <- rep(0, each = nrow(onlyNightSequence90))
        onlyNightSequence90$Count <- rep(0, each = nrow(onlyNightSequence90))
        for (i in onlyNightSequence90$RangInterest) { onlyNightSequence90$Count[which(onlyNightSequence90$RangInterest == i)] <- sum(data_seq_tax %>% filter(Cycle90 == "Night") %>% filter(RangInterest == i) %>% select(TotalNight))}
        onlyNightSequence90$Sum <- sum(onlyNightSequence90$Count)
        onlyCycleSequence90 <- rbind(onlyDaySequence90,onlyNightSequence90)
      ## Plot
        iz90 <- ggplot(onlyCycleSequence90, mapping = aes(y= value, x = Cycle90, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
          geom_label(aes(y = 106,label = paste(Sum,"séquences",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
          scale_fill_manual(values = palette)
        iz90 <- iz90 + labs(x="Cycles",y="Sequences (%)")
        print(iz90)
#
      # Zone -------------------------------------------------------------------
      ## Mixolimnion
        onlyMixolimnionSequence90 <- data_seq_tax %>% filter(Zone90 == "Mixolimnion") %>% select(ASV_Id,TotalMixolimnion,RangInterest)
        row.names(onlyMixolimnionSequence90)<-onlyMixolimnionSequence90$ASV_Id ; onlyMixolimnionSequence90 <- onlyMixolimnionSequence90 %>% select(-ASV_Id)
        onlyMixolimnionSequence90 <- onlyMixolimnionSequence90 %>% group_by(RangInterest) %>% summarise_all(sum)
        onlyMixolimnionSequence90$TotalMixolimnion <- onlyMixolimnionSequence90$TotalMixolimnion*100/sum(onlyMixolimnionSequence90$TotalMixolimnion)
        onlyMixolimnionSequence90 <- onlyMixolimnionSequence90 %>%
          arrange(desc(RangInterest)) %>%
          mutate(lab.ypos = cumsum(TotalMixolimnion) - 0.5*TotalMixolimnion)
        onlyMixolimnionSequence90$label <- paste(round(onlyMixolimnionSequence90$TotalMixolimnion,1), "%", sep = "")
        for (i in rownames(onlyMixolimnionSequence90)) {
          if (onlyMixolimnionSequence90[i,"label"] == "0%") { onlyMixolimnionSequence90[i,"label"] <- NA}}
        for (i in rownames(onlyMixolimnionSequence90)) {
          if (is.na(onlyMixolimnionSequence90[i,"label"]) == FALSE) { onlyMixolimnionSequence90[i,"label"] <- paste(onlyMixolimnionSequence90[i,"RangInterest"]," : ",onlyMixolimnionSequence90[i,"label"], sep = "")}}
        onlyMixolimnionSequence90$Zone90<- rep("Mixolimnion", each = nrow(onlyMixolimnionSequence90))
      ## Monimolimnion
        onlyMonimolimnionSequence90 <- data_seq_tax %>% filter(Zone90 == "Monimolimnion") %>% select(ASV_Id,TotalMonimolimnion,RangInterest)
        row.names(onlyMonimolimnionSequence90)<-onlyMonimolimnionSequence90$ASV_Id ; onlyMonimolimnionSequence90 <- onlyMonimolimnionSequence90 %>% select(-ASV_Id)
        onlyMonimolimnionSequence90 <- onlyMonimolimnionSequence90 %>% group_by(RangInterest) %>% summarise_all(sum)
        onlyMonimolimnionSequence90$TotalMonimolimnion <- onlyMonimolimnionSequence90$TotalMonimolimnion*100/sum(onlyMonimolimnionSequence90$TotalMonimolimnion)
        onlyMonimolimnionSequence90 <- onlyMonimolimnionSequence90 %>%
          arrange(desc(RangInterest)) %>%
          mutate(lab.ypos = cumsum(TotalMonimolimnion) - 0.5*TotalMonimolimnion)
        onlyMonimolimnionSequence90$label <- paste(round(onlyMonimolimnionSequence90$TotalMonimolimnion,1), "%", sep = "")
        for (i in rownames(onlyMonimolimnionSequence90)) {
          if (onlyMonimolimnionSequence90[i,"label"] == "0%") { onlyMonimolimnionSequence90[i,"label"] <- NA}}
        for (i in rownames(onlyMonimolimnionSequence90)) {
          if (is.na(onlyMonimolimnionSequence90[i,"label"]) == FALSE) { onlyMonimolimnionSequence90[i,"label"] <- paste(onlyMonimolimnionSequence90[i,"RangInterest"]," : ",onlyMonimolimnionSequence90[i,"label"], sep = "")}}
        onlyMonimolimnionSequence90$Zone90<- rep("Monimolimnion", each = nrow(onlyMonimolimnionSequence90))
      ## Zone
        colnames(onlyMixolimnionSequence90)[2]  <- "value"
        onlyMixolimnionSequence90$Sum <- rep(0, each = nrow(onlyMixolimnionSequence90))
        onlyMixolimnionSequence90$Count <- rep(0, each = nrow(onlyMixolimnionSequence90))
        for (i in onlyMixolimnionSequence90$RangInterest) { onlyMixolimnionSequence90$Count[which(onlyMixolimnionSequence90$RangInterest == i)] <- sum(data_seq_tax %>% filter(Zone90 == "Mixolimnion") %>% filter(RangInterest == i) %>% select(TotalMixolimnion))}
        onlyMixolimnionSequence90$Sum <- sum(onlyMixolimnionSequence90$Count)
        colnames(onlyMonimolimnionSequence90)[2]  <- "value"
        onlyMonimolimnionSequence90$Sum <- rep(0, each = nrow(onlyMonimolimnionSequence90))
        onlyMonimolimnionSequence90$Count <- rep(0, each = nrow(onlyMonimolimnionSequence90))
        for (i in onlyMonimolimnionSequence90$RangInterest) { onlyMonimolimnionSequence90$Count[which(onlyMonimolimnionSequence90$RangInterest == i)] <- sum(data_seq_tax %>% filter(Zone90 == "Monimolimnion") %>% filter(RangInterest == i) %>% select(TotalMonimolimnion))}
        onlyMonimolimnionSequence90$Sum <- sum(onlyMonimolimnionSequence90$Count)
        onlyZoneSequence90 <- rbind(onlyMixolimnionSequence90,onlyMonimolimnionSequence90)
      ## Plot
        jz90 <- ggplot(onlyZoneSequence90, mapping = aes(y= value, x = Zone90, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
          geom_label(aes(y = 106,label = paste(Sum,"séquences",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
          scale_fill_manual(values = palette)
        jz90 <- jz90 + labs(x="Zones",y="Séquences (%)")
        print(jz90)
#
      # Fraction -------------------------------------------------------------------
      ## Small
        onlySmallSequence90 <- data_seq_tax %>% filter(Fraction90 == "Small") %>% select(ASV_Id,TotalSmall,RangInterest)
        row.names(onlySmallSequence90)<-onlySmallSequence90$ASV_Id ; onlySmallSequence90 <- onlySmallSequence90 %>% select(-ASV_Id)
        onlySmallSequence90 <- onlySmallSequence90 %>% group_by(RangInterest) %>% summarise_all(sum)
        onlySmallSequence90$TotalSmall <- onlySmallSequence90$TotalSmall*100/sum(onlySmallSequence90$TotalSmall)
        onlySmallSequence90 <- onlySmallSequence90 %>%
          arrange(desc(RangInterest)) %>%
          mutate(lab.ypos = cumsum(TotalSmall) - 0.5*TotalSmall)
        onlySmallSequence90$label <- paste(round(onlySmallSequence90$TotalSmall,1), "%", sep = "")
        for (i in rownames(onlySmallSequence90)) {
          if (onlySmallSequence90[i,"label"] == "0%") { onlySmallSequence90[i,"label"] <- NA}}
        for (i in rownames(onlySmallSequence90)) {
          if (is.na(onlySmallSequence90[i,"label"]) == FALSE) { onlySmallSequence90[i,"label"] <- paste(onlySmallSequence90[i,"RangInterest"]," : ",onlySmallSequence90[i,"label"], sep = "")}}
        onlySmallSequence90$Fraction90<- rep("Small", each = nrow(onlySmallSequence90))
      ## Large
        onlyLargeSequence90 <- data_seq_tax %>% filter(Fraction90 == "Large") %>% select(ASV_Id,TotalLarge,RangInterest)
        row.names(onlyLargeSequence90)<-onlyLargeSequence90$ASV_Id ; onlyLargeSequence90 <- onlyLargeSequence90 %>% select(-ASV_Id)
        onlyLargeSequence90 <- onlyLargeSequence90 %>% group_by(RangInterest) %>% summarise_all(sum)
        onlyLargeSequence90$TotalLarge <- onlyLargeSequence90$TotalLarge*100/sum(onlyLargeSequence90$TotalLarge)
        onlyLargeSequence90 <- onlyLargeSequence90 %>%
          arrange(desc(RangInterest)) %>%
          mutate(lab.ypos = cumsum(TotalLarge) - 0.5*TotalLarge)
        onlyLargeSequence90$label <- paste(round(onlyLargeSequence90$TotalLarge,1), "%", sep = "")
        for (i in rownames(onlyLargeSequence90)) {
          if (onlyLargeSequence90[i,"label"] == "0%") { onlyLargeSequence90[i,"label"] <- NA}}
        for (i in rownames(onlyLargeSequence90)) {
          if (is.na(onlyLargeSequence90[i,"label"]) == FALSE) { onlyLargeSequence90[i,"label"] <- paste(onlyLargeSequence90[i,"RangInterest"]," : ",onlyLargeSequence90[i,"label"], sep = "")}}
        onlyLargeSequence90$Fraction90<- rep("Large", each = nrow(onlyLargeSequence90))
      ## Fraction
        colnames(onlySmallSequence90)[2]  <- "value"
        onlySmallSequence90$Sum <- rep(0, each = nrow(onlySmallSequence90))
        onlySmallSequence90$Count <- rep(0, each = nrow(onlySmallSequence90))
        for (i in onlySmallSequence90$RangInterest) { onlySmallSequence90$Count[which(onlySmallSequence90$RangInterest == i)] <- sum(data_seq_tax %>% filter(Fraction90 == "Small") %>% filter(RangInterest == i) %>% select(TotalSmall))}
        onlySmallSequence90$Sum <- sum(onlySmallSequence90$Count)
        colnames(onlyLargeSequence90)[2]  <- "value"
        onlyLargeSequence90$Sum <- rep(0, each = nrow(onlyLargeSequence90))
        onlyLargeSequence90$Count <- rep(0, each = nrow(onlyLargeSequence90))
        for (i in onlyLargeSequence90$RangInterest) { onlyLargeSequence90$Count[which(onlyLargeSequence90$RangInterest == i)] <- sum(data_seq_tax %>% filter(Fraction90 == "Large") %>% filter(RangInterest == i) %>% select(TotalLarge))}
        onlyLargeSequence90$Sum <- sum(onlyLargeSequence90$Count)
        onlyFractionSequence90 <- rbind(onlySmallSequence90,onlyLargeSequence90)
      ## Plot
        kz90 <- ggplot(onlyFractionSequence90, mapping = aes(y= value, x = Fraction90, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
          geom_label(aes(y = 106,label = paste(Sum,"séquences",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
          scale_fill_manual(values = palette)
        kz90 <- kz90 + labs(x="Fractions",y="Séquences (%)")
        print(kz90)
#
      # Coplot -------------------------------------------------------------------
        svglite("Hist-Taxomy/Sequence-only-90.svg",width = 12.00,height = 6.00)
        b_plot <- plot_grid(iz90,jz90,kz90,  ncol = 3, nrow = 1, rel_widths = c(3,3,3),rel_heights = c(3))
        print(b_plot)
        dev.off()  
#
# Small vs Large ----------------------------------------------------------
  data_seq_tax_Small <- data_seq_tax %>% filter(Small==1)
  data_seq_tax_Large <- data_seq_tax %>% filter(Large==2)
  data_asv_tax_Small <- data_asv_tax %>% filter(Small==1)
  data_asv_tax_Large <- data_asv_tax %>% filter(Large==2)
  for (type in c("Sequence","ASV")) {
    if (type == "ASV") { input_data_table <- asv_mat_rare }
    if (type == "Sequence") { input_data_table <- seq_mat_rare }
    ##
    MixolimnionXSmall <- samples_df %>% filter(Zone == "Mixolimnion") %>% filter(Fraction == "Small") %>% filter(Replicat == 1)
    MixolimnionXSmall<-MixolimnionXSmall$Amplicon
    MixolimnionXSmall_seq <- as.data.frame(rowSums(input_data_table %>% select(all_of(MixolimnionXSmall))))
    colnames(MixolimnionXSmall_seq) <- "MixolimnionXSmall"
    MixolimnionXSmall_seq$ASV_Id <- row.names(MixolimnionXSmall_seq)
    ##
    MixolimnionXLarge <- samples_df %>% filter(Zone == "Mixolimnion") %>% filter(Fraction == "Large") %>% filter(Replicat == 1)
    MixolimnionXLarge<-MixolimnionXLarge$Amplicon
    MixolimnionXLarge_seq <- as.data.frame(rowSums(input_data_table %>% select(all_of(MixolimnionXLarge))))
    colnames(MixolimnionXLarge_seq) <- "MixolimnionXLarge"
    MixolimnionXLarge_seq$ASV_Id <- row.names(MixolimnionXLarge_seq)
    ##
    MonimolimnionXSmall <- samples_df %>% filter(Zone == "Monimolimnion") %>% filter(Fraction == "Small") %>% filter(Replicat == 1)
    MonimolimnionXSmall<-MonimolimnionXSmall$Amplicon
    MonimolimnionXSmall_seq <- as.data.frame(rowSums(input_data_table %>% select(all_of(MonimolimnionXSmall))))
    colnames(MonimolimnionXSmall_seq) <- "MonimolimnionXSmall"
    MonimolimnionXSmall_seq$ASV_Id <- row.names(MonimolimnionXSmall_seq)
    ##
    MonimolimnionXLarge <- samples_df %>% filter(Zone == "Monimolimnion") %>% filter(Fraction == "Large") %>% filter(Replicat == 1)
    MonimolimnionXLarge<-MonimolimnionXLarge$Amplicon
    MonimolimnionXLarge_seq <- as.data.frame(rowSums(input_data_table %>% select(all_of(MonimolimnionXLarge))))
    colnames(MonimolimnionXLarge_seq) <- "MonimolimnionXLarge"
    MonimolimnionXLarge_seq$ASV_Id <- row.names(MonimolimnionXLarge_seq)
    ##
    MixolimnionXFraction_seq <- merge(x = MixolimnionXSmall_seq,y = MixolimnionXLarge_seq, by = "ASV_Id")
    AnxiqueXFraction_seq <- merge(x = MonimolimnionXSmall_seq,y = MonimolimnionXLarge_seq, by = "ASV_Id")
    ZoneXFraction_seq <- merge(x = MixolimnionXFraction_seq,y = AnxiqueXFraction_seq, by = "ASV_Id")
    if (type == "ASV") { ZoneXFraction_seq[,2:5][ZoneXFraction_seq[,2:5]>1]<-1}
    ZoneXFraction_seq_tax <- merge(data_seq_tax %>% select(ASV_Id,RangInterest,Division), ZoneXFraction_seq, by = "ASV_Id")
    ZoneXFraction_seq_tax <- ZoneXFraction_seq_tax %>% select(-ASV_Id)
    for (i in rownames(ZoneXFraction_seq_tax)) { if (ZoneXFraction_seq_tax[i,"Division"]=="Not Affiliated" && ZoneXFraction_seq_tax[i,"RangInterest"]!="Not Affiliated") { ZoneXFraction_seq_tax[i,"Division"] <- paste("Unaffiliated",ZoneXFraction_seq_tax[i,"RangInterest"],sep=" ")}
      if (ZoneXFraction_seq_tax[i,"Division"]=="Alveolata_X" ) { ZoneXFraction_seq_tax[i,"Division"] <- paste("Unaffiliated",ZoneXFraction_seq_tax[i,"RangInterest"],sep=" ")}}
    ZoneXFraction_seq_tax <- ZoneXFraction_seq_tax %>% group_by(RangInterest,Division) %>% summarise_all(sum)
    Orderspe <- ZoneXFraction_seq_tax$Division
    ZoneXFraction_seq_tax_melt <-melt(ZoneXFraction_seq_tax,id.vars = c("RangInterest","Division"))
    ZoneXFraction_seq_tax_melt$Sum <- 0
    for (i in row.names(ZoneXFraction_seq_tax_melt)) { var <- ZoneXFraction_seq_tax_melt[i,"variable"]
    ZoneXFraction_seq_tax_melt[i,"Sum"] <- sum(ZoneXFraction_seq_tax_melt %>% filter(variable == var) %>% select(value))
    ZoneXFraction_seq_tax_melt[i,"Proportion"] <- ZoneXFraction_seq_tax_melt[i,"value"]*100/ZoneXFraction_seq_tax_melt[i,"Sum"]}
    ZoneXFraction_seq_tax_melt <- ZoneXFraction_seq_tax_melt %>% mutate()
    ZoneXFraction_seq_tax_melt <- separate(ZoneXFraction_seq_tax_melt,"variable",c("Zone","Fraction"),sep="X")
    assign(paste("ZoneXFraction_seq_tax_melt",type,sep = "_"),ZoneXFraction_seq_tax_melt)
  }
## Plot Sequence
  svglite("Hist-Taxomy/ZoneXFraction-Sequence-Total.svg",width = 7.00,height = 9.00)
  Total_ZoneXFraction_Sequence_fig <- ggplot(ZoneXFraction_seq_tax_melt_Sequence, mapping = aes(y= Proportion, x = Zone,fill = factor(Division,level=Orderspe), group = factor(Division,level=Orderspe)), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity",color="black") +
    facet_grid(.~Fraction,scales="free") + scale_fill_manual(values = paletspe) + guides(fill=guide_legend(ncol=1)) +
    labs(x="Zones",y="Sequences %",fill="Divisions") + theme_unique_art()
  print(Total_ZoneXFraction_Sequence_fig)
  dev.off()
  write.table(ZoneXFraction_seq_tax_melt_Sequence,file="dataTables/ZoneXFraction-Sequence-Total-Table.csv",row.names = F, quote = F, sep ="\t")
## Plot ASV
  svglite("Hist-Taxomy/ZoneXFraction-ASV-Total.svg",width = 7.00,height = 9.00)
  Total_ZoneXFraction_ASV_fig <- ggplot(ZoneXFraction_seq_tax_melt_ASV, mapping = aes(y= Proportion, x = Zone, fill = factor(Division,level=Orderspe), group = factor(Division,level=Orderspe)), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity",color="black") +
    facet_grid(.~Fraction,scales="free") + scale_fill_manual(values = paletspe) + guides(fill=guide_legend(ncol=1)) +
    labs(x="Zones",y="ASVs %",fill="Divisions") + geom_label(aes(y = 106,label = paste(Sum," ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B") + theme_unique_art()
  print(Total_ZoneXFraction_ASV_fig)
  dev.off()
  write.table(ZoneXFraction_seq_tax_melt_ASV,file="dataTables/ZoneXFraction-ASV-Total-Table.csv",row.names = F, quote = F, sep ="\t")
#    
  # Balance plot ------------------------------------------------------------
    # Sequence ---------------------------------------------------------------------
    ## Mixolimnion
      balance_totalZoneXFractionSequence_Mixolimnion <- ZoneXFraction_seq_tax_melt_Sequence %>% filter(Zone == "Mixolimnion") %>% select(RangInterest,Division,value,Zone,Fraction,Sum) %>% group_by(RangInterest,Division,Zone,Fraction,Sum) %>% summarise_all(sum)
    ## Monimolimnion
      balance_totalZoneXFractionSequence_Monimolimnion <- ZoneXFraction_seq_tax_melt_Sequence %>% filter(Zone == "Monimolimnion") %>% select(RangInterest,Division,value,Zone,Fraction,Sum) %>% group_by(RangInterest,Division,Zone,Fraction) %>% summarise_all(sum)
    ## Bind
      balance_totalZoneXFractionSequence<-balance_totalZoneXFractionSequence_Mixolimnion %>% select(-Sum)
      for (i in row.names(balance_totalZoneXFractionSequence)) { balance_totalZoneXFractionSequence[i,"value"] <- foldchange(balance_totalZoneXFractionSequence_Mixolimnion[i,"value"],balance_totalZoneXFractionSequence_Monimolimnion[i,"value"])
      balance_totalZoneXFractionSequence[i,"value"] <- foldchange2logratio(balance_totalZoneXFractionSequence[i,"value"][[1]])
      if (balance_totalZoneXFractionSequence_Mixolimnion[i,"value"] < balance_totalZoneXFractionSequence_Monimolimnion[i,"value"]) {balance_totalZoneXFractionSequence[i,"Zone"] <- "Monimolimnion"}
      if (balance_totalZoneXFractionSequence_Mixolimnion[i,"value"] > balance_totalZoneXFractionSequence_Monimolimnion[i,"value"]) {balance_totalZoneXFractionSequence[i,"Zone"] <- "Mixolimnion"}
      balance_totalZoneXFractionSequence[i,"label"] <- paste(max(balance_totalZoneXFractionSequence_Mixolimnion[i,"value"],balance_totalZoneXFractionSequence_Monimolimnion[i,"value"]),"VS",min(balance_totalZoneXFractionSequence_Mixolimnion[i,"value"],balance_totalZoneXFractionSequence_Monimolimnion[i,"value"]),sep=" ")
      balance_totalZoneXFractionSequence[i,"Abundance"] <- max(balance_totalZoneXFractionSequence_Mixolimnion[i,"value"],balance_totalZoneXFractionSequence_Monimolimnion[i,"value"])}
    ## plot
      balance_palet <- paletteer_d("ggthemes::Nuriel_Stone",n=2, direction=-1)
      svglite("Hist-Taxomy/Sequence-Balance-ZoneXFraction.svg",width = 14.00,height = 8.00)
      balxi <- ggplot(balance_totalZoneXFractionSequence, aes(x = factor(Division,level=rev(Orderspe)), y = value)) +
        geom_bar(aes(fill = Zone,alpha=log10(Abundance)), stat = 'identity',color="black") +  facet_grid(.~Fraction) +
        coord_flip() +  # horizontal bars
        geom_text(aes(label = label ,y = 0,vjust = ifelse(value >= 0, 0.5, 0.5), hjust = ifelse(value>= 0, 1.2, -0.2)),color = "black",size = 3,show.legend = FALSE, fill = "white") +
        #geom_label(aes(label = round(Abundance),y = 0, vjust = ifelse(value >= 0, 0.5, 0.5),hjust = ifelse(value>= 0, 1.2, -0.2)),color = "black",size = 3,show.legend = FALSE, fill = "white") +
        labs(x="",y="Sequences abundance log2-ratio",fill="Zones") +
        theme(axis.text.y = element_blank()) + theme_unique_art() + scale_fill_manual(values=balance_palet) + guides(alpha = "none")
      print(balxi)
      dev.off()
    ## Coplot
      blank <- grid.rect(gp=gpar(col="white"))
    ## Seq
      balxim <- balxi + theme(axis.text.y = element_blank(),legend.position = c(0.08, 0.92),legend.background = element_rect(fill = "white", color = "black"))
      balxim <- plot_grid(blank,balxim, ncol = 1, nrow = 2, rel_widths = c(2),rel_heights = c(1,35))
      ymibal_plot <- plot_grid(Total_ZoneXFraction_Sequence_fig,balxim, ncol = 2, nrow = 1, labels=c("C","D"),rel_widths = c(2,3.3),rel_heights = c(3))
      svglite("Hist-Taxomy/Sequence-ALL-ZoneXFraction.svg",width = 18.00,height = 8.60)
      print(ymibal_plot)
      dev.off()
#
    # ASV ---------------------------------------------------------------------
    ## Mixolimnion
      balance_totalZoneXFractionASV_Mixolimnion <- ZoneXFraction_seq_tax_melt_ASV %>% filter(Zone == "Mixolimnion") %>% select(RangInterest,Division,value,Zone,Fraction,Sum) %>% group_by(RangInterest,Division,Zone,Fraction,Sum) %>% summarise_all(sum)
    ## Monimolimnion
      balance_totalZoneXFractionASV_Monimolimnion <- ZoneXFraction_seq_tax_melt_ASV %>% filter(Zone == "Monimolimnion") %>% select(RangInterest,Division,value,Zone,Fraction,Sum) %>% group_by(RangInterest,Division,Zone,Fraction) %>% summarise_all(sum)
    ## Bind
      balance_totalZoneXFractionASV<-balance_totalZoneXFractionASV_Mixolimnion %>% select(-Sum)
      for (i in row.names(balance_totalZoneXFractionASV)) { balance_totalZoneXFractionASV[i,"value"] <- foldchange(balance_totalZoneXFractionASV_Mixolimnion[i,"value"],balance_totalZoneXFractionASV_Monimolimnion[i,"value"])
      balance_totalZoneXFractionASV[i,"value"] <- foldchange2logratio(balance_totalZoneXFractionASV[i,"value"][[1]])
      if (balance_totalZoneXFractionASV_Mixolimnion[i,"value"] < balance_totalZoneXFractionASV_Monimolimnion[i,"value"]) {balance_totalZoneXFractionASV[i,"Zone"] <- "Monimolimnion"}
      if (balance_totalZoneXFractionASV_Mixolimnion[i,"value"] > balance_totalZoneXFractionASV_Monimolimnion[i,"value"]) {balance_totalZoneXFractionASV[i,"Zone"] <- "Mixolimnion"}
      balance_totalZoneXFractionASV[i,"label"] <- paste(max(balance_totalZoneXFractionASV_Mixolimnion[i,"value"],balance_totalZoneXFractionASV_Monimolimnion[i,"value"]),"VS",min(balance_totalZoneXFractionASV_Mixolimnion[i,"value"],balance_totalZoneXFractionASV_Monimolimnion[i,"value"]),sep=" ")
      balance_totalZoneXFractionASV[i,"Abundance"] <- max(balance_totalZoneXFractionASV_Mixolimnion[i,"value"],balance_totalZoneXFractionASV_Monimolimnion[i,"value"])}
    ## plot
      svglite("Hist-Taxomy/ASV-Balance-ZoneXFraction.svg",width = 14.00,height = 8.00)
      balxj <- ggplot(balance_totalZoneXFractionASV, aes(x = factor(Division,level=rev(Orderspe)), y = value)) +
        geom_bar(aes(fill = Zone,alpha=log10(Abundance)), stat = 'identity',color="black") +  facet_grid(.~Fraction) +
        coord_flip() +  # horizontal bars
        geom_text(aes(label = label ,y = 0,vjust = ifelse(value >= 0, 0.5, 0.5), hjust = ifelse(value>= 0, 1.2, -0.2)),color = "black",size = 3,show.legend = FALSE, fill = "white") +
        #geom_label(aes(label = round(Abundance),y = 0, vjust = ifelse(value >= 0, 0.5, 0.5),hjust = ifelse(value>= 0, 1.2, -0.2)),color = "black",size = 3,show.legend = FALSE, fill = "white") +
        labs(x="",y="ASVs log2-ratio",fill="Zones") +
        theme(axis.text.y = element_blank()) + theme_unique_art() + scale_fill_manual(values=balance_palet) + guides(alpha = "none")
      print(balxj)
      dev.off()
    ## Coplot
      blank <- grid.rect(gp=gpar(col="white"))
    ## ASV
      balxjm <- balxj + theme(axis.text.y = element_blank(),legend.position = c(0.08, 0.92),legend.background = element_rect(fill = "white", color = "black"))
      balxjm <- plot_grid(blank,balxjm, ncol = 1, nrow = 2, rel_widths = c(2),rel_heights = c(1,35))
      ymjbal_plot <- plot_grid(Total_ZoneXFraction_ASV_fig,balxjm, ncol = 2, nrow = 1, labels=c("A","B"),rel_widths = c(2,3.3),rel_heights = c(3))
      svglite("Hist-Taxomy/ASV-ALL-ZoneXFraction.svg",width = 18.00,height = 8.60)
      print(ymjbal_plot)
      dev.off()
#
    # cowplot -----------------------------------------------------------------
      ymjibal_plot <- plot_grid(ymjbal_plot,ymibal_plot, ncol = 1, nrow = 2,rel_widths = c(2),rel_heights = c(3,3))
      svglite("Hist-Taxomy/ALL-ALL-ZoneXFraction.svg",width = 18.00,height = 17.2)
      print(ymjibal_plot)
      dev.off()
      
# Polar ------------------------------------------------------
  # Séquence ----------------------------------------------------------------
  ## Supergroup
  ### Total column
    Polar_seq <- data_seq_tax
    Polar_seq$TotalAnnée <- 0
    for (i in rownames(Polar_seq)) { Polar_seq[i,"TotalAnnée"] <- Polar_seq[i,"Total04"] + Polar_seq[i,"Total06"] + Polar_seq[i,"Total09"] + Polar_seq[i,"Total11"]}
    Polar_seq <- Polar_seq %>% select(TotalAnnée,RangInterest)
  ### Total Figure
    Polar_seq <- Polar_seq %>% group_by(RangInterest) %>% summarise_all(sum)
    for (i in rownames(Polar_seq)) { Polar_seq[i,"Total"] <- (Polar_seq[i,"TotalAnnée"] * 100) / sum(Polar_seq$TotalAnnée)}
  ### Label
    Polar_seq <- Polar_seq %>%
      arrange(desc(RangInterest)) %>%
      mutate(lab.ypos = cumsum(Total) - 0.5*Total)
    Polar_seq$label <- paste(round(Polar_seq$Total,1), "%", sep = "")
    for (i in rownames(Polar_seq)) {
      if (Polar_seq[i,"label"] == "0%") { Polar_seq[i,"label"] <- NA}}
    for (i in rownames(Polar_seq)) {
      if (is.na(Polar_seq[i,"label"]) == FALSE) { Polar_seq[i,"label"] <- paste(Polar_seq[i,"RangInterest"]," : ",Polar_seq[i,"label"], sep = "")}}
    Polar_seq$color <- rev(c("#c62728cc","#6b1b9acc","#273493cc","#0277bdcc","#558b30cc","#f9a825cc","#ef6c00cc","#25252590","#880e4fcc","#206064cc"))
  ### Figure
    svglite("Composition/Polar-Total-Supergroup-seq.svg",width = 4.50,height = 4.50)
    ax <- ggplot(Polar_seq, mapping = aes(y= Total, x = 2, fill = RangInterest), Rowv = NA, col = colMain, scale = "column") +
      geom_bar(stat="identity", color = "white", width = 1,fill=Polar_seq$color) + coord_polar("y") + 
      geom_label_repel(aes(y = lab.ypos,label = label),color = "white",size = 3,segment.color = "black",show.legend = FALSE, nudge_x = 0.5,fill=Polar_seq$color) +
      scale_y_continuous(limits=c(0,sum(Polar_seq %>% select(Total))))+
      xlim(1,3) +
      theme_unique_darkbis() + 
      #facet_wrap( ~ variable , nrow = 4) +
      labs(y = "Sequences",x="") #+ scale_fill_manual(values = rev(palette))
    print(ax)
    dev.off()
    write.table(Polar_seq,file="Composition/Table-Total-Supergroup-seq.csv",quote = F,row.names = F,sep="\t")
  ## Division
  ### Total column
    Polar_seq_div <- data_seq_tax
    Polar_seq_div$TotalAnnée <- 0
    for (i in rownames(Polar_seq_div)) { Polar_seq_div[i,"TotalAnnée"] <- Polar_seq_div[i,"Total04"] + Polar_seq_div[i,"Total06"] + Polar_seq_div[i,"Total09"] + Polar_seq_div[i,"Total11"]}
    Polar_seq_div <- Polar_seq_div %>% select(TotalAnnée,RangInterest,Division)
  ### Total Figure
    Polar_seq_div <- Polar_seq_div %>% group_by(RangInterest,Division) %>% summarise_all(sum)
    for (i in rownames(Polar_seq_div)) { if (Polar_seq_div[i,"Division"]=="Not Affiliated" && Polar_seq_div[i,"RangInterest"]!="Not Affiliated" ) { Polar_seq_div[i,"Division"] <- paste("Unaffiliated",Polar_seq_div[i,"RangInterest"],sep=" ")}}
    for (i in rownames(Polar_seq_div)) { if (length(grep("_X",Polar_seq_div[i,"Division"]))>0) { Polar_seq_div[i,"Division"] <- paste("Unaffiliated",Polar_seq_div[i,"RangInterest"],sep=" ")}}
    Polar_seq_div <- Polar_seq_div %>% group_by(RangInterest,Division) %>% summarise_all(sum)
    Polar_seq_div <- as.data.frame(Polar_seq_div) %>% select(-RangInterest)
    for (i in rownames(Polar_seq_div)) { Polar_seq_div[i,"Total"] <- (Polar_seq_div[i,"TotalAnnée"] * 100) / sum(Polar_seq_div$TotalAnnée)}
  ### Label
    Polar_seq_div <- Polar_seq_div %>%
      arrange(-row_number()) %>%
      mutate(lab.ypos = cumsum(Total) - 0.5*Total)
    Polar_seq_div$label <- paste(round(Polar_seq_div$Total,1), "%", sep = "")
    for (i in rownames(Polar_seq_div)) {
      if (Polar_seq_div[i,"label"] == "0%") { Polar_seq_div[i,"label"] <- NA}}
    for (i in rownames(Polar_seq_div)) {
      if (is.na(Polar_seq_div[i,"label"]) == FALSE) { Polar_seq_div[i,"label"] <- paste(Polar_seq_div[i,"Division"]," : ",Polar_seq_div[i,"label"], sep = "")}}
    Polar_seq_div$color <- rev(paletspe)
    ### Figure
      svglite("Composition/Polar-Total-Division-seq.svg",width = 5.50,height = 5.50)
      axd <- ggplot(Polar_seq_div, mapping = aes(y= Total, x = 2, fill = factor(Division,level=Polar_seq_div$Division)), Rowv = NA, col = colMain, scale = "column") +
        geom_bar(stat="identity", color = "white", width = 1,fill=Polar_seq_div$color) + coord_polar("y") + 
        geom_label_repel(aes(y = lab.ypos,label = label),color = "white",size = 3,segment.color = "black",show.legend = FALSE, nudge_x = 0.5,fill=Polar_seq_div$color) +
        scale_y_continuous(limits=c(0,sum(Polar_seq_div %>% select(Total)))) +
        xlim(0.5,2.5) +
        theme_unique_darkbis() + 
        #facet_wrap( ~ variable , nrow = 4) +
        labs(y = "Séquences",x="") #+ scale_fill_paletteer_d("ggsci::default_igv")
      print(axd)
      dev.off()
      write.table(Polar_seq_div,file="Composition/Table-Total-Division-seq.csv",quote = F,row.names = F,sep="\t")
#
  #Try test
      axCop <- ggplot(Rowv = NA, col = colMain, scale = "column") +
        geom_bar(data = Polar_seq,  mapping = aes(y= Total, x = 2, fill = RangInterest), stat="identity", color = "white", width = 1,fill=Polar_seq$color) + 
        #xlim(1,3) +
        theme_unique_darkbis() +
        labs(y = "Sequences",x="") +
        geom_bar(data = Polar_seq_div, mapping = aes(y= Total, x = 3, fill = factor(Division,level=Polar_seq_div$Division)), stat="identity", color = "white", width = 1,fill=Polar_seq_div$color) + coord_polar("y") + 
        #geom_label_repel(data = Polar_seq_div, aes(x=3,y = lab.ypos,label = label),color = "white",size = 3,segment.color = "black",show.legend = FALSE, nudge_y = 0.5,fill=Polar_seq_div$color) +
        geom_text_repel(data = Polar_seq_div, aes(x=3,y = lab.ypos,label = label),fontface= "bold",color = "black",size = 3,segment.color = "black",show.legend = FALSE, nudge_x = 0.8, force_pull = 0) +
        scale_y_continuous(limits=c(0,sum(Polar_seq_div %>% select(Total)))) +
        xlim(1,3.8)
        print(axCop)
#
  # ASV ----------------------------------------------------------------
  ## Supergroup
  ### Total column
    Polar_asv <- data_asv_tax
    Polar_asv$TotalAnnée <- 0
    for (i in rownames(Polar_asv)) { if (Polar_asv[i,"Total04"] + Polar_asv[i,"Total06"] + Polar_asv[i,"Total09"] + Polar_asv[i,"Total11"] > 0) {Polar_asv[i,"TotalAnnée"] <- 1}}
    Polar_asv <- Polar_asv %>% select(TotalAnnée,RangInterest)
  ### Total Figure
    Polar_asv <- Polar_asv %>% group_by(RangInterest) %>% summarise_all(sum)
    for (i in rownames(Polar_asv)) { Polar_asv[i,"Total"] <- (Polar_asv[i,"TotalAnnée"] * 100) / sum(Polar_asv$TotalAnnée)}
  ### Label
    Polar_asv <- Polar_asv %>%
      arrange(desc(RangInterest)) %>%
      mutate(lab.ypos = cumsum(Total) - 0.5*Total)
    Polar_asv$label <- paste(round(Polar_asv$Total,1), "%", sep = "")
    for (i in rownames(Polar_asv)) {
      if (Polar_asv[i,"label"] == "0%") { Polar_asv[i,"label"] <- NA}}
    for (i in rownames(Polar_asv)) {
      if (is.na(Polar_asv[i,"label"]) == FALSE) { Polar_asv[i,"label"] <- paste(Polar_asv[i,"RangInterest"]," : ",Polar_asv[i,"label"], sep = "")}}
    Polar_asv$color <- rev(c("#c62728cc","#6b1b9acc","#273493cc","#0277bdcc","#558b30cc","#f9a825cc","#ef6c00cc","#25252590","#880e4fcc","#206064cc"))
  ### Figure
    svglite("Composition/Polar-Total-Supergroup-asv.svg",width = 4.50,height = 4.50)
    bx <- ggplot(Polar_asv, mapping = aes(y= Total, x = 2, fill = RangInterest), Rowv = NA, col = colMain, scale = "column") +
      geom_bar(stat="identity", color = "white", width = 1,fill=Polar_asv$color) + coord_polar("y") + 
      geom_label_repel(aes(y = lab.ypos,label = label),color = "white",size = 3,segment.color = "black",show.legend = FALSE, nudge_x = 0.5,fill=Polar_asv$color) +
      scale_y_continuous(limits=c(0,sum(Polar_asv %>% select(Total))))+
      xlim(1,3) +
      theme_unique_darkbis() + 
      #facet_wrap( ~ variable , nrow = 4) +
      labs(y = "ASVs",x="") #+ scale_fill_manual(values = rev(palette))
    print(bx)
    dev.off()
    write.table(Polar_asv,file="Composition/Table-Total-Supergroup-asv.csv",quote = F,row.names = F,sep="\t")
  ## Division
  ### Total column
    Polar_asv_div <- data_asv_tax
    Polar_asv_div$TotalAnnée <- 0
    for (i in rownames(Polar_asv_div)) { if (Polar_asv_div[i,"Total04"] + Polar_asv_div[i,"Total06"] + Polar_asv_div[i,"Total09"] + Polar_asv_div[i,"Total11"] > 0) {Polar_asv_div[i,"TotalAnnée"] <- 1}}
    Polar_asv_div <- Polar_asv_div %>% select(TotalAnnée,RangInterest,Division)
  ### Total Figure
    Polar_asv_div <- Polar_asv_div %>% group_by(RangInterest,Division) %>% summarise_all(sum)
    for (i in rownames(Polar_asv_div)) { if (Polar_asv_div[i,"Division"]=="Not Affiliated" && Polar_asv_div[i,"RangInterest"]!="Not Affiliated") { Polar_asv_div[i,"Division"] <- paste("Unaffiliated",Polar_asv_div[i,"RangInterest"],sep=" ")}}
    for (i in rownames(Polar_asv_div)) { if (length(grep("_X",Polar_asv_div[i,"Division"]))>0) { Polar_asv_div[i,"Division"] <- paste("Unaffiliated",Polar_asv_div[i,"RangInterest"],sep=" ")}}
    Polar_asv_div <- Polar_asv_div %>% group_by(RangInterest,Division) %>% summarise_all(sum)
    Polar_asv_div <- as.data.frame(Polar_asv_div) %>% select(-RangInterest)
    for (i in rownames(Polar_asv_div)) { Polar_asv_div[i,"Total"] <- (Polar_asv_div[i,"TotalAnnée"] * 100) / sum(Polar_asv_div$TotalAnnée)}
  ### Label
    Polar_asv_div <- Polar_asv_div %>%
      arrange(-row_number()) %>%
      mutate(lab.ypos = cumsum(Total) - 0.5*Total)
    Polar_asv_div$label <- paste(round(Polar_asv_div$Total,1), "%", sep = "")
    for (i in rownames(Polar_asv_div)) {
      if (Polar_asv_div[i,"label"] == "0%") { Polar_asv_div[i,"label"] <- NA}}
    for (i in rownames(Polar_asv_div)) {
      if (is.na(Polar_asv_div[i,"label"]) == FALSE) { Polar_asv_div[i,"label"] <- paste(Polar_asv_div[i,"Division"]," : ",Polar_asv_div[i,"label"], sep = "")}}
    Polar_asv_div$color <- rev(paletspe)
  ### Figure
    svglite("Composition/Polar-Total-Division-asv.svg",width = 5.50,height = 5.50)
    bxd <- ggplot(Polar_asv, mapping = aes(y= Total, x = 2, fill = factor(Division,level=Polar_asv$Division)), Rowv = NA, col = colMain, scale = "column") +
      geom_bar(stat="identity", color = "white", width = 1,fill=Polar_asv$color) + coord_polar("y") + 
      geom_label_repel(aes(y = lab.ypos,label = label),color = "white",size = 3,segment.color = "black",show.legend = FALSE, nudge_x = 0.5,fill=Polar_asv$color) +
      scale_y_continuous(limits=c(0,sum(Polar_asv %>% select(Total)))) +
      xlim(0.5,2.5) +
      theme_unique_darkbis() + 
      #facet_wrap( ~ variable , nrow = 4) +
      labs(y = "ASVs",x="") #+ scale_fill_paletteer_d("ggsci::default_igv")
    print(bxd)
    dev.off()
    write.table(Polar_asv,file="Composition/Table-Total-Division-asv.csv",quote = F,row.names = F,sep="\t")
  #Try test
    bxCop <- ggplot(Rowv = NA, col = colMain, scale = "column") +
      geom_bar(data = Polar_asv,  mapping = aes(y= Total, x = 2, fill = RangInterest), stat="identity", color = "white", width = 1,fill=Polar_asv$color) + 
      #xlim(1,3) +
      theme_unique_darkbis() +
      labs(y = "ASVs",x="") +
      geom_bar(data = Polar_asv_div, mapping = aes(y= Total, x = 3, fill = factor(Division,level=Polar_asv_div$Division)), stat="identity", color = "white", width = 1,fill=Polar_asv_div$color) + coord_polar("y") + 
      #geom_label_repel(data = Polar_asv_div, aes(x=3,y = lab.ypos,label = label),color = "white",size = 3,segment.color = "black",show.legend = FALSE, nudge_y = 0.5,fill=Polar_asv_div$color) +
      geom_text_repel(data = Polar_asv_div, aes(x=3,y = lab.ypos,label = label),fontface= "bold",color = "black",size = 3,segment.color = "black",show.legend = FALSE, nudge_x = 0.8, force_pull = 0) +
      scale_y_continuous(limits=c(0,sum(Polar_asv_div %>% select(Total)))) +
      xlim(1,3.8)
    print(bxCop)
#
  # Coplot ------------------------------------------------------------------
    svglite("Composition/Table-Total-Full.svg",width = 6.50,height = 13.00)
    b_plot <- plot_grid(bxCop,axCop, labels = c("A","B"), ncol = 1, nrow = 2, rel_widths = c(3),rel_heights = c(3,3))
    print(b_plot)
    dev.off()

# Table ASVs majoritaires -------------------------------------------------
  # Table Only ---------------------------------------------------------------------
    for (i in c("Cycle","Fraction","Zone","Cycle90","Fraction90","Zone90")) {
      mdt <- unique(data_seq_tax[,i])
      mdt <- mdt[!mdt %in% "Shared"]
      cdt <- str_split(i,"90")[[1]][1]
      for (f in mdt) {
        assign(paste0("only",f,"Table"),data_seq_tax %>% filter(get(i) == f) %>% select(ASV_Id,paste0("Total",f),RangInterest,i))
        assign(paste0("only",f,"Table"), get(paste0("only",f,"Table")) %>% filter(get(paste0("Total",f)) > 0.01*sum(get(paste0("only",f,"Table"))[,paste0("Total",f)])))
        dt <- get(paste0("only",f,"Table"))
        colnames(dt)[2] <- "value"
        assign(paste0("only",f,"Table"),dt)
      }
      assign(paste0("only",i,"Table"),rbind(get(paste0("only",mdt[1],"Table")),get(paste0("only",mdt[2],"Table"))))
      assign(paste0("only",i,"Table"),merge(get(paste0("only",i,"Table")),tax_tablemix, by = "ASV_Id"))
      write.table(get(paste0("only",i,"Table")), file = paste0("dataTables/Only_",i,".txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
      if (i %in% c("Cycle90","Fraction90","Zone90")) { 
        write.table(get(paste0("only",cdt,"Sequence90")), file = paste0("dataTables/Only_",i,"_Sequence.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
      }
      else {
      write.table(get(paste0("only",i,"ASV")), file = paste0("dataTables/Only_",i,"_ASV.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
      write.table(get(paste0("only",i,"Sequence")), file = paste0("dataTables/Only_",i,"_Sequence.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
      }
    }
  # Table Total --------------------------------------------------------------------- 
  
    for (i in c("Cycle","Fraction","Zone")) {
      mdt <- unique(data_seq_tax[,i])
      mdt <- mdt[!mdt %in% "Shared"]
      for (f in mdt) {
        assign(paste0("total",f,"Table"),data_seq_tax %>% select(ASV_Id,paste0("Total",f),RangInterest,i))
        dt <- get(paste0("total",f,"Table"))
        dt[,i] <- f
        assign(paste0("total",f,"Table"),dt)
        assign(paste0("total",f,"Table"), get(paste0("total",f,"Table")) %>% filter(get(paste0("Total",f)) > 0.05*sum(get(paste0("total",f,"Table"))[,paste0("Total",f)])))
        dt <- get(paste0("total",f,"Table"))
        colnames(dt)[2] <- "value"
        assign(paste0("total",f,"Table"),dt)
      }
      assign(paste0("total",i,"Table"),rbind(get(paste0("total",mdt[1],"Table")),get(paste0("total",mdt[2],"Table"))))
      assign(paste0("total",i,"Table"),merge(get(paste0("total",i,"Table")),tax_tablemix, by = "ASV_Id"))
      write.table(get(paste0("total",i,"Table")), file = paste0("dataTables/Total_",i,".txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
      write.table(get(paste0("total",i,"ASV")), file = paste0("dataTables/Total_",i,"_ASV.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
      write.table(get(paste0("total",i,"Sequence")), file = paste0("dataTables/Total_",i,"_Sequence.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    }
  ## Period
    total04Table <- data_seq_tax %>% select(ASV_Id,Total04,RangInterest)
    total04Table$Periods <- rep(04,each = nrow(total04Table))
    total04Table <- total04Table %>% filter(Total04 > 0.05*sum(total04Table$Total04))
    colnames(total04Table)[2]  <- "value"
    ##
    total06Table <- data_seq_tax %>% select(ASV_Id,Total06,RangInterest)
    total06Table$Periods <- rep(06,each = nrow(total06Table))
    total06Table <- total06Table %>% filter(Total06 > 0.05*sum(total06Table$Total06))
    colnames(total06Table)[2]  <- "value"
    ##
    total09Table <- data_seq_tax %>% select(ASV_Id,Total09,RangInterest)
    total09Table$Periods <- rep(09,each = nrow(total09Table))
    total09Table <- total09Table %>% filter(Total09 > 0.05*sum(total09Table$Total09))
    colnames(total09Table)[2]  <- "value"
    ##
    total11Table <- data_seq_tax %>% select(ASV_Id,Total11,RangInterest)
    total11Table$Periods <- rep(11,each = nrow(total11Table))
    total11Table <- total11Table %>% filter(Total11 > 0.05*sum(total11Table$Total11))
    colnames(total11Table)[2]  <- "value"
    ##
    totalPeriodTable <- rbind(total04Table,total06Table,total09Table,total11Table)
    totalPeriodTable <- merge(totalPeriodTable,tax_tablemix, by = "ASV_Id")
    ##
    write.table(totalPeriodTable, file = "dataTables/Total_Periods.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(totalPeriodASV, file = "dataTables/Total_Periods_ASV.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(totalPeriodSequence, file = "dataTables/Total_Periods_Sequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
#
# Hist & Table Mix ----------------------------------------------------------------
  data_seq_tax$RangInterestMix <- data_seq_tax$Division 
  data_asv_tax$RangInterestMix <- data_asv_tax$Division 
  # Hist Seq Total Mix ---------------------------------------------------------------------
    # Cycle -------------------------------------------------------------------
    ## Day
      totalDaySequence <- data_seq_tax %>% select(ASV_Id,TotalDay,RangInterest,RangInterestMix)
      row.names(totalDaySequence)<-totalDaySequence$ASV_Id ; totalDaySequence <- totalDaySequence %>% select(-ASV_Id)
      totalDaySequence <- totalDaySequence %>% group_by(RangInterest,RangInterestMix) %>% summarise_all(sum)
      for (i in rownames(totalDaySequence)) { if (totalDaySequence[i,"RangInterestMix"]=="Not Affiliated" && totalDaySequence[i,"RangInterest"]!="Not Affiliated") { totalDaySequence[i,"RangInterestMix"] <- paste("Unaffiliated",totalDaySequence[i,"RangInterest"],sep=" ")}
        if (totalDaySequence[i,"RangInterestMix"]=="Alveolata_X" ) { totalDaySequence[i,"RangInterestMix"] <- paste("Unaffiliated",totalDaySequence[i,"RangInterest"],sep=" ")}}
      totalDaySequence<-as.data.frame(totalDaySequence) %>% select(-RangInterest)
      totalDaySequence$Count <- totalDaySequence$TotalDay*100/sum(totalDaySequence$TotalDay)
      totalDaySequence <- totalDaySequence %>%
        arrange(desc(RangInterestMix)) %>%
        mutate(lab.ypos = cumsum(Count) - 0.5*Count)
      totalDaySequence$label <- paste(round(totalDaySequence$Count,1), "%", sep = "")
      for (i in rownames(totalDaySequence)) {
        if (totalDaySequence[i,"label"] == "0%") { totalDaySequence[i,"label"] <- NA}}
      for (i in rownames(totalDaySequence)) {
        if (is.na(totalDaySequence[i,"label"]) == FALSE) { totalDaySequence[i,"label"] <- paste(totalDaySequence[i,"RangInterestMix"]," : ",totalDaySequence[i,"label"], sep = "")}}
      totalDaySequence$Cycle<- rep("Day", each = nrow(totalDaySequence))
    ## Night
      totalNightSequence <- data_seq_tax %>% select(ASV_Id,TotalNight,RangInterest,RangInterestMix)
      row.names(totalNightSequence)<-totalNightSequence$ASV_Id ; totalNightSequence <- totalNightSequence %>% select(-ASV_Id)
      totalNightSequence <- totalNightSequence %>% group_by(RangInterest,RangInterestMix) %>% summarise_all(sum)
      for (i in rownames(totalNightSequence)) { if (totalNightSequence[i,"RangInterestMix"]=="Not Affiliated" && totalNightSequence[i,"RangInterest"]!="Not Affiliated") { totalNightSequence[i,"RangInterestMix"] <- paste("Unaffiliated",totalNightSequence[i,"RangInterest"],sep=" ")}
        if (totalNightSequence[i,"RangInterestMix"]=="Alveolata_X" ) { totalNightSequence[i,"RangInterestMix"] <- paste("Unaffiliated",totalNightSequence[i,"RangInterest"],sep=" ")}}
      totalNightSequence<-as.data.frame(totalNightSequence) %>% select(-RangInterest)
      totalNightSequence$Count <- totalNightSequence$TotalNight*100/sum(totalNightSequence$TotalNight)
      totalNightSequence <- totalNightSequence %>%
        arrange(desc(RangInterestMix)) %>%
        mutate(lab.ypos = cumsum(Count) - 0.5*Count)
      totalNightSequence$label <- paste(round(totalNightSequence$Count,1), "%", sep = "")
      for (i in rownames(totalNightSequence)) {
        if (totalNightSequence[i,"label"] == "0%") { totalNightSequence[i,"label"] <- NA}}
      for (i in rownames(totalNightSequence)) {
        if (is.na(totalNightSequence[i,"label"]) == FALSE) { totalNightSequence[i,"label"] <- paste(totalNightSequence[i,"RangInterestMix"]," : ",totalNightSequence[i,"label"], sep = "")}}
      totalNightSequence$Cycle<- rep("Night", each = nrow(totalNightSequence))
    ## Cycle
      colnames(totalDaySequence)[2]  <- "value"
      totalDaySequence$Sum <- rep(0, each = nrow(totalDaySequence))
      totalDaySequence$Sum <- sum(totalDaySequence$value)
      colnames(totalNightSequence)[2]  <- "value"
      totalNightSequence$Sum <- rep(0, each = nrow(totalNightSequence))
      totalNightSequence$Sum <- sum(totalNightSequence$value)
      totalCycleSequence <- rbind(totalDaySequence,totalNightSequence)
      totalCycleSequence <- totalCycleSequence %>% mutate(percent=paste("(",round(Sum*100/(colSums(statRarefy %>% select(`apRarefy-Sequence`))),1)," %)",sep =""))
    ## Plot
      iy <- ggplot(totalCycleSequence, mapping = aes(y= Count, x = Cycle, fill = factor(RangInterestMix,level=Orderspe), group = factor(RangInterestMix,level=Orderspe)), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity",color="black") +
        scale_fill_manual(values = paletspe) + theme_unique_art()  + guides(fill=guide_legend(ncol=1)) +
        geom_label(aes(y = 106,label = paste(Sum,"sequences",percent,sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") + labs(x="Cycle",y="Sequences (%)",fill="Division")
      legendSequence <- get_legend(iy)
      iy <- iy + theme(legend.position = "none")
      print(iy)
#
    # Zone -------------------------------------------------------------------
    ## Mixolimnion
      totalMixolimnionSequence <- data_seq_tax %>% select(ASV_Id,TotalMixolimnion,RangInterest,RangInterestMix)
      row.names(totalMixolimnionSequence)<-totalMixolimnionSequence$ASV_Id ; totalMixolimnionSequence <- totalMixolimnionSequence %>% select(-ASV_Id)
      totalMixolimnionSequence <- totalMixolimnionSequence %>% group_by(RangInterest,RangInterestMix) %>% summarise_all(sum)
      for (i in rownames(totalMixolimnionSequence)) { if (totalMixolimnionSequence[i,"RangInterestMix"]=="Not Affiliated" && totalMixolimnionSequence[i,"RangInterest"]!="Not Affiliated") { totalMixolimnionSequence[i,"RangInterestMix"] <- paste("Unaffiliated",totalMixolimnionSequence[i,"RangInterest"],sep=" ")}
        if (totalMixolimnionSequence[i,"RangInterestMix"]=="Alveolata_X" ) { totalMixolimnionSequence[i,"RangInterestMix"] <- paste("Unaffiliated",totalMixolimnionSequence[i,"RangInterest"],sep=" ")}}
      totalMixolimnionSequence<-as.data.frame(totalMixolimnionSequence) %>% select(-RangInterest)
      totalMixolimnionSequence$Count <- totalMixolimnionSequence$TotalMixolimnion*100/sum(totalMixolimnionSequence$TotalMixolimnion)
      totalMixolimnionSequence <- totalMixolimnionSequence %>%
        arrange(desc(RangInterestMix)) %>%
        mutate(lab.ypos = cumsum(Count) - 0.5*Count)
      totalMixolimnionSequence$label <- paste(round(totalMixolimnionSequence$Count,1), "%", sep = "")
      for (i in rownames(totalMixolimnionSequence)) {
        if (totalMixolimnionSequence[i,"label"] == "0%") { totalMixolimnionSequence[i,"label"] <- NA}}
      for (i in rownames(totalMixolimnionSequence)) {
        if (is.na(totalMixolimnionSequence[i,"label"]) == FALSE) { totalMixolimnionSequence[i,"label"] <- paste(totalMixolimnionSequence[i,"RangInterestMix"]," : ",totalMixolimnionSequence[i,"label"], sep = "")}}
      totalMixolimnionSequence$Zone<- rep("Mixolimnion", each = nrow(totalMixolimnionSequence))
    ## Monimolimnion
      totalMonimolimnionSequence <- data_seq_tax %>% select(ASV_Id,TotalMonimolimnion,RangInterest,RangInterestMix)
      row.names(totalMonimolimnionSequence)<-totalMonimolimnionSequence$ASV_Id ; totalMonimolimnionSequence <- totalMonimolimnionSequence %>% select(-ASV_Id)
      totalMonimolimnionSequence <- totalMonimolimnionSequence %>% group_by(RangInterest,RangInterestMix) %>% summarise_all(sum)
      for (i in rownames(totalMonimolimnionSequence)) { if (totalMonimolimnionSequence[i,"RangInterestMix"]=="Not Affiliated" && totalMonimolimnionSequence[i,"RangInterest"]!="Not Affiliated") { totalMonimolimnionSequence[i,"RangInterestMix"] <- paste("Unaffiliated",totalMonimolimnionSequence[i,"RangInterest"],sep=" ")}
        if (totalMonimolimnionSequence[i,"RangInterestMix"]=="Alveolata_X" ) { totalMonimolimnionSequence[i,"RangInterestMix"] <- paste("Unaffiliated",totalMonimolimnionSequence[i,"RangInterest"],sep=" ")}}
      totalMonimolimnionSequence<-as.data.frame(totalMonimolimnionSequence) %>% select(-RangInterest)
      totalMonimolimnionSequence$Count <- totalMonimolimnionSequence$TotalMonimolimnion*100/sum(totalMonimolimnionSequence$TotalMonimolimnion)
      totalMonimolimnionSequence <- totalMonimolimnionSequence %>%
        arrange(desc(RangInterestMix)) %>%
        mutate(lab.ypos = cumsum(Count) - 0.5*Count)
      totalMonimolimnionSequence$label <- paste(round(totalMonimolimnionSequence$Count,1), "%", sep = "")
      for (i in rownames(totalMonimolimnionSequence)) {
        if (totalMonimolimnionSequence[i,"label"] == "0%") { totalMonimolimnionSequence[i,"label"] <- NA}}
      for (i in rownames(totalMonimolimnionSequence)) {
        if (is.na(totalMonimolimnionSequence[i,"label"]) == FALSE) { totalMonimolimnionSequence[i,"label"] <- paste(totalMonimolimnionSequence[i,"RangInterestMix"]," : ",totalMonimolimnionSequence[i,"label"], sep = "")}}
      totalMonimolimnionSequence$Zone<- rep("Monimolimnion", each = nrow(totalMonimolimnionSequence))
    ## Zone
      colnames(totalMixolimnionSequence)[2]  <- "value"
      totalMixolimnionSequence$Sum <- rep(0, each = nrow(totalMixolimnionSequence))
      totalMixolimnionSequence$Sum <- sum(totalMixolimnionSequence$value)
      colnames(totalMonimolimnionSequence)[2]  <- "value"
      totalMonimolimnionSequence$Sum <- rep(0, each = nrow(totalMonimolimnionSequence))
      totalMonimolimnionSequence$Sum <- sum(totalMonimolimnionSequence$value)
      totalZoneSequence <- rbind(totalMixolimnionSequence,totalMonimolimnionSequence)
      totalZoneSequence <- totalZoneSequence %>% mutate(percent=paste("(",round(Sum*100/(colSums(statRarefy %>% select(`apRarefy-Sequence`))),1)," %)",sep =""))
    ## Plot
      jy <- ggplot(totalZoneSequence, mapping = aes(y= Count, x = Zone, fill = factor(RangInterestMix,level=Orderspe), group = factor(RangInterestMix,level=Orderspe)), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity",color="black") +
        scale_fill_manual(values = paletspe) + theme_unique_art() +
        geom_label(aes(y = 106,label = paste(Sum,"sequences",percent,sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
      jy <- jy + labs(x="Zone",y="Sequences (%)") + theme(legend.position = "none") + theme(axis.title.y = element_blank())
      print(jy)
#
    # Fraction -------------------------------------------------------------------
    ## Small
      totalSmallSequence <- data_seq_tax %>% select(ASV_Id,TotalSmall,RangInterest,RangInterestMix)
      row.names(totalSmallSequence)<-totalSmallSequence$ASV_Id ; totalSmallSequence <- totalSmallSequence %>% select(-ASV_Id)
      totalSmallSequence <- totalSmallSequence %>% group_by(RangInterest,RangInterestMix) %>% summarise_all(sum)
      for (i in rownames(totalSmallSequence)) { if (totalSmallSequence[i,"RangInterestMix"]=="Not Affiliated" && totalSmallSequence[i,"RangInterest"]!="Not Affiliated") { totalSmallSequence[i,"RangInterestMix"] <- paste("Unaffiliated",totalSmallSequence[i,"RangInterest"],sep=" ")}
        if (totalSmallSequence[i,"RangInterestMix"]=="Alveolata_X" ) { totalSmallSequence[i,"RangInterestMix"] <- paste("Unaffiliated",totalSmallSequence[i,"RangInterest"],sep=" ")}}
      totalSmallSequence<-as.data.frame(totalSmallSequence) %>% select(-RangInterest)
      totalSmallSequence$Count <- totalSmallSequence$TotalSmall*100/sum(totalSmallSequence$TotalSmall)
      totalSmallSequence <- totalSmallSequence %>%
        arrange(desc(RangInterestMix)) %>%
        mutate(lab.ypos = cumsum(Count) - 0.5*Count)
      totalSmallSequence$label <- paste(round(totalSmallSequence$Count,1), "%", sep = "")
      for (i in rownames(totalSmallSequence)) {
        if (totalSmallSequence[i,"label"] == "0%") { totalSmallSequence[i,"label"] <- NA}}
      for (i in rownames(totalSmallSequence)) {
        if (is.na(totalSmallSequence[i,"label"]) == FALSE) { totalSmallSequence[i,"label"] <- paste(totalSmallSequence[i,"RangInterestMix"]," : ",totalSmallSequence[i,"label"], sep = "")}}
      totalSmallSequence$Fraction<- rep("Small", each = nrow(totalSmallSequence))
    ## Large
      totalLargeSequence <- data_seq_tax %>% select(ASV_Id,TotalLarge,RangInterest,RangInterestMix)
      row.names(totalLargeSequence)<-totalLargeSequence$ASV_Id ; totalLargeSequence <- totalLargeSequence %>% select(-ASV_Id)
      totalLargeSequence <- totalLargeSequence %>% group_by(RangInterest,RangInterestMix) %>% summarise_all(sum)
      for (i in rownames(totalLargeSequence)) { if (totalLargeSequence[i,"RangInterestMix"]=="Not Affiliated" && totalLargeSequence[i,"RangInterest"]!="Not Affiliated") { totalLargeSequence[i,"RangInterestMix"] <- paste("Unaffiliated",totalLargeSequence[i,"RangInterest"],sep=" ")}
        if (totalLargeSequence[i,"RangInterestMix"]=="Alveolata_X" ) { totalLargeSequence[i,"RangInterestMix"] <- paste("Unaffiliated",totalLargeSequence[i,"RangInterest"],sep=" ")}}
      totalLargeSequence<-as.data.frame(totalLargeSequence) %>% select(-RangInterest)
      totalLargeSequence$Count <- totalLargeSequence$TotalLarge*100/sum(totalLargeSequence$TotalLarge)
      totalLargeSequence <- totalLargeSequence %>%
        arrange(desc(RangInterestMix)) %>%
        mutate(lab.ypos = cumsum(Count) - 0.5*Count)
      totalLargeSequence$label <- paste(round(totalLargeSequence$Count,1), "%", sep = "")
      for (i in rownames(totalLargeSequence)) {
        if (totalLargeSequence[i,"label"] == "0%") { totalLargeSequence[i,"label"] <- NA}}
      for (i in rownames(totalLargeSequence)) {
        if (is.na(totalLargeSequence[i,"label"]) == FALSE) { totalLargeSequence[i,"label"] <- paste(totalLargeSequence[i,"RangInterestMix"]," : ",totalLargeSequence[i,"label"], sep = "")}}
      totalLargeSequence$Fraction<- rep("Large", each = nrow(totalLargeSequence))
    ## Fraction
      colnames(totalSmallSequence)[2]  <- "value"
      totalSmallSequence$Sum <- rep(0, each = nrow(totalSmallSequence))
      totalSmallSequence$Sum <- sum(totalSmallSequence$value)
      colnames(totalLargeSequence)[2]  <- "value"
      totalLargeSequence$Sum <- rep(0, each = nrow(totalLargeSequence))
      totalLargeSequence$Sum <- sum(totalLargeSequence$value)
      totalFractionSequence <- rbind(totalSmallSequence,totalLargeSequence)
      totalFractionSequence <- totalFractionSequence %>% mutate(percent=paste("(",round(Sum*100/(colSums(statRarefy %>% select(`apRarefy-Sequence`))),1)," %)",sep =""))
    ## Plot
      ky <- ggplot(totalFractionSequence, mapping = aes(y= Count, x = Fraction, fill = factor(RangInterestMix,level=Orderspe), group = factor(RangInterestMix,level=Orderspe)), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity",color="black") +
        scale_fill_manual(values = paletspe) + theme_unique_art() +
        geom_label(aes(y = 106,label = paste(Sum,"sequences",percent,sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
      ky <- ky + labs(x="Fraction",y="Sequences (%)") + theme(legend.position = "none") + theme(axis.title.y = element_blank())
      print(ky)
#
    # Period -------------------------------------------------------------------
    ## Avril
      total04Sequence <- data_seq_tax %>% select(ASV_Id,Total04,RangInterest,RangInterestMix)
      row.names(total04Sequence)<-total04Sequence$ASV_Id ; total04Sequence <- total04Sequence %>% select(-ASV_Id)
      total04Sequence <- total04Sequence %>% group_by(RangInterest,RangInterestMix) %>% summarise_all(sum)
      for (i in rownames(total04Sequence)) { if (total04Sequence[i,"RangInterestMix"]=="Not Affiliated" && total04Sequence[i,"RangInterest"]!="Not Affiliated") { total04Sequence[i,"RangInterestMix"] <- paste("Unaffiliated",total04Sequence[i,"RangInterest"],sep=" ")}
        if (total04Sequence[i,"RangInterestMix"]=="Alveolata_X" ) { total04Sequence[i,"RangInterestMix"] <- paste("Unaffiliated",total04Sequence[i,"RangInterest"],sep=" ")}}
      total04Sequence<-as.data.frame(total04Sequence) %>% select(-RangInterest)
      total04Sequence$Count <- total04Sequence$Total04*100/sum(total04Sequence$Total04)
      total04Sequence <- total04Sequence %>%
        arrange(desc(RangInterestMix)) %>%
        mutate(lab.ypos = cumsum(Count) - 0.5*Count)
      total04Sequence$label <- paste(round(total04Sequence$Count,1), "%", sep = "")
      for (i in rownames(total04Sequence)) {
        if (total04Sequence[i,"label"] == "0%") { total04Sequence[i,"label"] <- NA}}
      for (i in rownames(total04Sequence)) {
        if (is.na(total04Sequence[i,"label"]) == FALSE) { total04Sequence[i,"label"] <- paste(total04Sequence[i,"RangInterestMix"]," : ",total04Sequence[i,"label"], sep = "")}}
      total04Sequence$Period<- rep("04", each = nrow(total04Sequence))
    ## Juin
      total06Sequence <- data_seq_tax %>% select(ASV_Id,Total06,RangInterest,RangInterestMix)
      row.names(total06Sequence)<-total06Sequence$ASV_Id ; total06Sequence <- total06Sequence %>% select(-ASV_Id)
      total06Sequence <- total06Sequence %>% group_by(RangInterest,RangInterestMix) %>% summarise_all(sum)
      for (i in rownames(total06Sequence)) { if (total06Sequence[i,"RangInterestMix"]=="Not Affiliated" && total06Sequence[i,"RangInterest"]!="Not Affiliated") { total06Sequence[i,"RangInterestMix"] <- paste("Unaffiliated",total06Sequence[i,"RangInterest"],sep=" ")}
        if (total06Sequence[i,"RangInterestMix"]=="Alveolata_X" ) { total06Sequence[i,"RangInterestMix"] <- paste("Unaffiliated",total06Sequence[i,"RangInterest"],sep=" ")}}
      total06Sequence<-as.data.frame(total06Sequence) %>% select(-RangInterest)
      total06Sequence$Count <- total06Sequence$Total06*100/sum(total06Sequence$Total06)
      total06Sequence <- total06Sequence %>%
        arrange(desc(RangInterestMix)) %>%
        mutate(lab.ypos = cumsum(Count) - 0.5*Count)
      total06Sequence$label <- paste(round(total06Sequence$Count,1), "%", sep = "")
      for (i in rownames(total06Sequence)) {
        if (total06Sequence[i,"label"] == "0%") { total06Sequence[i,"label"] <- NA}}
      for (i in rownames(total06Sequence)) {
        if (is.na(total06Sequence[i,"label"]) == FALSE) { total06Sequence[i,"label"] <- paste(total06Sequence[i,"RangInterestMix"]," : ",total06Sequence[i,"label"], sep = "")}}
      total06Sequence$Period<- rep("06", each = nrow(total06Sequence))
    ## Septembre
      total09Sequence <- data_seq_tax %>% select(ASV_Id,Total09,RangInterest,RangInterestMix)
      row.names(total09Sequence)<-total09Sequence$ASV_Id ; total09Sequence <- total09Sequence %>% select(-ASV_Id)
      total09Sequence <- total09Sequence %>% group_by(RangInterest,RangInterestMix) %>% summarise_all(sum)
      for (i in rownames(total09Sequence)) { if (total09Sequence[i,"RangInterestMix"]=="Not Affiliated" && total09Sequence[i,"RangInterest"]!="Not Affiliated") { total09Sequence[i,"RangInterestMix"] <- paste("Unaffiliated",total09Sequence[i,"RangInterest"],sep=" ")}
        if (total09Sequence[i,"RangInterestMix"]=="Alveolata_X" ) { total09Sequence[i,"RangInterestMix"] <- paste("Unaffiliated",total09Sequence[i,"RangInterest"],sep=" ")}}
      total09Sequence<-as.data.frame(total09Sequence) %>% select(-RangInterest)
      total09Sequence$Count <- total09Sequence$Total09*100/sum(total09Sequence$Total09)
      total09Sequence <- total09Sequence %>%
        arrange(desc(RangInterestMix)) %>%
        mutate(lab.ypos = cumsum(Count) - 0.5*Count)
      total09Sequence$label <- paste(round(total09Sequence$Count,1), "%", sep = "")
      for (i in rownames(total09Sequence)) {
        if (total09Sequence[i,"label"] == "0%") { total09Sequence[i,"label"] <- NA}}
      for (i in rownames(total09Sequence)) {
        if (is.na(total09Sequence[i,"label"]) == FALSE) { total09Sequence[i,"label"] <- paste(total09Sequence[i,"RangInterestMix"]," : ",total09Sequence[i,"label"], sep = "")}}
      total09Sequence$Period<- rep("09", each = nrow(total09Sequence))
    ## Octobre
      total11Sequence <- data_seq_tax %>% select(ASV_Id,Total11,RangInterest,RangInterestMix)
      row.names(total11Sequence)<-total11Sequence$ASV_Id ; total11Sequence <- total11Sequence %>% select(-ASV_Id)
      total11Sequence <- total11Sequence %>% group_by(RangInterest,RangInterestMix) %>% summarise_all(sum)
      for (i in rownames(total11Sequence)) { if (total11Sequence[i,"RangInterestMix"]=="Not Affiliated" && total11Sequence[i,"RangInterest"]!="Not Affiliated") { total11Sequence[i,"RangInterestMix"] <- paste("Unaffiliated",total11Sequence[i,"RangInterest"],sep=" ")}
        if (total11Sequence[i,"RangInterestMix"]=="Alveolata_X" ) { total11Sequence[i,"RangInterestMix"] <- paste("Unaffiliated",total11Sequence[i,"RangInterest"],sep=" ")}}
      total11Sequence<-as.data.frame(total11Sequence) %>% select(-RangInterest)
      total11Sequence$Count <- total11Sequence$Total11*100/sum(total11Sequence$Total11)
      total11Sequence <- total11Sequence %>%
        arrange(desc(RangInterestMix)) %>%
        mutate(lab.ypos = cumsum(Count) - 0.5*Count)
      total11Sequence$label <- paste(round(total11Sequence$Count,1), "%", sep = "")
      for (i in rownames(total11Sequence)) {
        if (total11Sequence[i,"label"] == "0%") { total11Sequence[i,"label"] <- NA}}
      for (i in rownames(total11Sequence)) {
        if (is.na(total11Sequence[i,"label"]) == FALSE) { total11Sequence[i,"label"] <- paste(total11Sequence[i,"RangInterestMix"]," : ",total11Sequence[i,"label"], sep = "")}}
      total11Sequence$Period<- rep("11", each = nrow(total11Sequence))
    ## Periods
      colnames(total04Sequence)[2]  <- "value"
      total04Sequence$Sum <- rep(0, each = nrow(total04Sequence))
      total04Sequence$Sum <- sum(total04Sequence$value)
      colnames(total06Sequence)[2]  <- "value"
      total06Sequence$Sum <- rep(0, each = nrow(total06Sequence))
      total06Sequence$Sum <- sum(total06Sequence$value)
      colnames(total09Sequence)[2]  <- "value"
      total09Sequence$Sum <- rep(0, each = nrow(total09Sequence))
      total09Sequence$Sum <- sum(total09Sequence$value)
      colnames(total11Sequence)[2]  <- "value"
      total11Sequence$Sum <- rep(0, each = nrow(total11Sequence))
      total11Sequence$Sum <- sum(total11Sequence$value)
      totalPeriodSequence <- rbind(total04Sequence,total06Sequence,total09Sequence,total11Sequence)
      totalPeriodSequence <- totalPeriodSequence %>% mutate(percent=paste("(",round(Sum*100/(colSums(statRarefy %>% select(`apRarefy-Sequence`))),1)," %)",sep =""))
    ## Plot
      ly <- ggplot(totalPeriodSequence, mapping = aes(y= Count, x = Period, fill = factor(RangInterestMix,level=Orderspe), group = factor(RangInterestMix,level=Orderspe)), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity",color="black") +
        scale_fill_manual(values = paletspe) + theme_unique_art() +
        geom_label(aes(y = 106,label = paste(Sum,"sequences",percent,sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
      ly <- ly + labs(x="Period",y="Sequences (%)") + theme(legend.position = "none") + theme(axis.title.y = element_blank())
      print(ly)
#
    # Coplot -------------------------------------------------------------------
      svglite("Hist-Taxomy/SequenceMIX-Total.svg",width = 14.00,height = 10.00)
      b_plot <- plot_grid(iy,jy,ky,ly, legendSequence, ncol = 5, nrow = 1, rel_widths = c(3,3,3,6,3),rel_heights = c(3))
      print(b_plot)
      dev.off()
#
  # Hist ASV Total Mix ---------------------------------------------------------------------

    # Cycle -------------------------------------------------------------------
    ## Day
      totalDayASV <- data_asv_tax %>% select(ASV_Id,TotalDay,RangInterest,RangInterestMix)
      row.names(totalDayASV)<-totalDayASV$ASV_Id ; totalDayASV <- totalDayASV %>% select(-ASV_Id)
      totalDayASV <- totalDayASV %>% group_by(RangInterest,RangInterestMix) %>% summarise_all(sum)
      for (i in rownames(totalDayASV)) { if (totalDayASV[i,"RangInterestMix"]=="Not Affiliated" && totalDayASV[i,"RangInterest"]!="Not Affiliated") { totalDayASV[i,"RangInterestMix"] <- paste("Unaffiliated",totalDayASV[i,"RangInterest"],sep=" ")}
        if (totalDayASV[i,"RangInterestMix"]=="Alveolata_X" ) { totalDayASV[i,"RangInterestMix"] <- paste("Unaffiliated",totalDayASV[i,"RangInterest"],sep=" ")}}
      totalDayASV<-as.data.frame(totalDayASV) %>% select(-RangInterest)
      totalDayASV$Count <- totalDayASV$TotalDay*100/sum(totalDayASV$TotalDay)
      totalDayASV <- totalDayASV %>%
        arrange(desc(RangInterestMix)) %>%
        mutate(lab.ypos = cumsum(Count) - 0.5*Count)
      totalDayASV$label <- paste(round(totalDayASV$Count,1), "%", sep = "")
      for (i in rownames(totalDayASV)) {
        if (totalDayASV[i,"label"] == "0%") { totalDayASV[i,"label"] <- NA}}
      for (i in rownames(totalDayASV)) {
        if (is.na(totalDayASV[i,"label"]) == FALSE) { totalDayASV[i,"label"] <- paste(totalDayASV[i,"RangInterestMix"]," : ",totalDayASV[i,"label"], sep = "")}}
      totalDayASV$Cycle<- rep("Day", each = nrow(totalDayASV))
    ## Night
      totalNightASV <- data_asv_tax %>% select(ASV_Id,TotalNight,RangInterest,RangInterestMix)
      row.names(totalNightASV)<-totalNightASV$ASV_Id ; totalNightASV <- totalNightASV %>% select(-ASV_Id)
      totalNightASV <- totalNightASV %>% group_by(RangInterest,RangInterestMix) %>% summarise_all(sum)
      for (i in rownames(totalNightASV)) { if (totalNightASV[i,"RangInterestMix"]=="Not Affiliated" && totalNightASV[i,"RangInterest"]!="Not Affiliated") { totalNightASV[i,"RangInterestMix"] <- paste("Unaffiliated",totalNightASV[i,"RangInterest"],sep=" ")}
        if (totalNightASV[i,"RangInterestMix"]=="Alveolata_X" ) { totalNightASV[i,"RangInterestMix"] <- paste("Unaffiliated",totalNightASV[i,"RangInterest"],sep=" ")}}
      totalNightASV<-as.data.frame(totalNightASV) %>% select(-RangInterest)
      totalNightASV$Count <- totalNightASV$TotalNight*100/sum(totalNightASV$TotalNight)
      totalNightASV <- totalNightASV %>%
        arrange(desc(RangInterestMix)) %>%
        mutate(lab.ypos = cumsum(Count) - 0.5*Count)
      totalNightASV$label <- paste(round(totalNightASV$Count,1), "%", sep = "")
      for (i in rownames(totalNightASV)) {
        if (totalNightASV[i,"label"] == "0%") { totalNightASV[i,"label"] <- NA}}
      for (i in rownames(totalNightASV)) {
        if (is.na(totalNightASV[i,"label"]) == FALSE) { totalNightASV[i,"label"] <- paste(totalNightASV[i,"RangInterestMix"]," : ",totalNightASV[i,"label"], sep = "")}}
      totalNightASV$Cycle<- rep("Night", each = nrow(totalNightASV))
    ## Cycle
      colnames(totalDayASV)[2]  <- "value"
      totalDayASV$Sum <- rep(0, each = nrow(totalDayASV))
      totalDayASV$Sum <- sum(totalDayASV$value)
      colnames(totalNightASV)[2]  <- "value"
      totalNightASV$Sum <- rep(0, each = nrow(totalNightASV))
      totalNightASV$Sum <- sum(totalNightASV$value)
      totalCycleASV <- rbind(totalDayASV,totalNightASV)
    ## Plot
      iy <- ggplot(totalCycleASV, mapping = aes(y= Count, x = Cycle, fill = factor(RangInterestMix,level=Orderspe), group = factor(RangInterestMix,level=Orderspe)), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity",color="black") +
        scale_fill_manual(values = paletspe) + theme_unique_art()  + guides(fill=guide_legend(ncol=1)) +
        geom_label(aes(y = 106,label = paste(Sum,"ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") + labs(x="Cycle",y="ASVs (%)",fill="Division")
      legendASV <- get_legend(iy)
      iy <- iy + theme(legend.position = "none")
      print(iy)
#
    # Zone -------------------------------------------------------------------
    ## Mixolimnion
      totalMixolimnionASV <- data_asv_tax %>% select(ASV_Id,TotalMixolimnion,RangInterest,RangInterestMix)
      row.names(totalMixolimnionASV)<-totalMixolimnionASV$ASV_Id ; totalMixolimnionASV <- totalMixolimnionASV %>% select(-ASV_Id)
      totalMixolimnionASV <- totalMixolimnionASV %>% group_by(RangInterest,RangInterestMix) %>% summarise_all(sum)
      for (i in rownames(totalMixolimnionASV)) { if (totalMixolimnionASV[i,"RangInterestMix"]=="Not Affiliated" && totalMixolimnionASV[i,"RangInterest"]!="Not Affiliated") { totalMixolimnionASV[i,"RangInterestMix"] <- paste("Unaffiliated",totalMixolimnionASV[i,"RangInterest"],sep=" ")}
        if (totalMixolimnionASV[i,"RangInterestMix"]=="Alveolata_X" ) { totalMixolimnionASV[i,"RangInterestMix"] <- paste("Unaffiliated",totalMixolimnionASV[i,"RangInterest"],sep=" ")}}
      totalMixolimnionASV<-as.data.frame(totalMixolimnionASV) %>% select(-RangInterest)
      totalMixolimnionASV$Count <- totalMixolimnionASV$TotalMixolimnion*100/sum(totalMixolimnionASV$TotalMixolimnion)
      totalMixolimnionASV <- totalMixolimnionASV %>%
        arrange(desc(RangInterestMix)) %>%
        mutate(lab.ypos = cumsum(Count) - 0.5*Count)
      totalMixolimnionASV$label <- paste(round(totalMixolimnionASV$Count,1), "%", sep = "")
      for (i in rownames(totalMixolimnionASV)) {
        if (totalMixolimnionASV[i,"label"] == "0%") { totalMixolimnionASV[i,"label"] <- NA}}
      for (i in rownames(totalMixolimnionASV)) {
        if (is.na(totalMixolimnionASV[i,"label"]) == FALSE) { totalMixolimnionASV[i,"label"] <- paste(totalMixolimnionASV[i,"RangInterestMix"]," : ",totalMixolimnionASV[i,"label"], sep = "")}}
      totalMixolimnionASV$Zone<- rep("Mixolimnion", each = nrow(totalMixolimnionASV))
    ## Monimolimnion
      totalMonimolimnionASV <- data_asv_tax %>% select(ASV_Id,TotalMonimolimnion,RangInterest,RangInterestMix)
      row.names(totalMonimolimnionASV)<-totalMonimolimnionASV$ASV_Id ; totalMonimolimnionASV <- totalMonimolimnionASV %>% select(-ASV_Id)
      totalMonimolimnionASV <- totalMonimolimnionASV %>% group_by(RangInterest,RangInterestMix) %>% summarise_all(sum)
      for (i in rownames(totalMonimolimnionASV)) { if (totalMonimolimnionASV[i,"RangInterestMix"]=="Not Affiliated" && totalMonimolimnionASV[i,"RangInterest"]!="Not Affiliated") { totalMonimolimnionASV[i,"RangInterestMix"] <- paste("Unaffiliated",totalMonimolimnionASV[i,"RangInterest"],sep=" ")}
        if (totalMonimolimnionASV[i,"RangInterestMix"]=="Alveolata_X" ) { totalMonimolimnionASV[i,"RangInterestMix"] <- paste("Unaffiliated",totalMonimolimnionASV[i,"RangInterest"],sep=" ")}}
      totalMonimolimnionASV<-as.data.frame(totalMonimolimnionASV) %>% select(-RangInterest)
      totalMonimolimnionASV$Count <- totalMonimolimnionASV$TotalMonimolimnion*100/sum(totalMonimolimnionASV$TotalMonimolimnion)
      totalMonimolimnionASV <- totalMonimolimnionASV %>%
        arrange(desc(RangInterestMix)) %>%
        mutate(lab.ypos = cumsum(Count) - 0.5*Count)
      totalMonimolimnionASV$label <- paste(round(totalMonimolimnionASV$Count,1), "%", sep = "")
      for (i in rownames(totalMonimolimnionASV)) {
        if (totalMonimolimnionASV[i,"label"] == "0%") { totalMonimolimnionASV[i,"label"] <- NA}}
      for (i in rownames(totalMonimolimnionASV)) {
        if (is.na(totalMonimolimnionASV[i,"label"]) == FALSE) { totalMonimolimnionASV[i,"label"] <- paste(totalMonimolimnionASV[i,"RangInterestMix"]," : ",totalMonimolimnionASV[i,"label"], sep = "")}}
      totalMonimolimnionASV$Zone<- rep("Monimolimnion", each = nrow(totalMonimolimnionASV))
    ## Zone
      colnames(totalMixolimnionASV)[2]  <- "value"
      totalMixolimnionASV$Sum <- rep(0, each = nrow(totalMixolimnionASV))
      totalMixolimnionASV$Sum <- sum(totalMixolimnionASV$value)
      colnames(totalMonimolimnionASV)[2]  <- "value"
      totalMonimolimnionASV$Sum <- rep(0, each = nrow(totalMonimolimnionASV))
      totalMonimolimnionASV$Sum <- sum(totalMonimolimnionASV$value)
      totalZoneASV <- rbind(totalMixolimnionASV,totalMonimolimnionASV)
    ## plot
      jy <- ggplot(totalZoneASV, mapping = aes(y= Count, x = Zone, fill = factor(RangInterestMix,level=Orderspe), group = factor(RangInterestMix,level=Orderspe)), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity",color="black") +
        scale_fill_manual(values = paletspe) + theme_unique_art() +
        geom_label(aes(y = 106,label = paste(Sum,"ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
      jy <- jy + labs(x="Zone",y="ASVs (%)") + theme(legend.position = "none") + theme(axis.title.y = element_blank())
      print(jy)
#
    # Fraction -------------------------------------------------------------------
    ## Small
      totalSmallASV <- data_asv_tax %>% select(ASV_Id,TotalSmall,RangInterest,RangInterestMix)
      row.names(totalSmallASV)<-totalSmallASV$ASV_Id ; totalSmallASV <- totalSmallASV %>% select(-ASV_Id)
      totalSmallASV <- totalSmallASV %>% group_by(RangInterest,RangInterestMix) %>% summarise_all(sum)
      for (i in rownames(totalSmallASV)) { if (totalSmallASV[i,"RangInterestMix"]=="Not Affiliated" && totalSmallASV[i,"RangInterest"]!="Not Affiliated") { totalSmallASV[i,"RangInterestMix"] <- paste("Unaffiliated",totalSmallASV[i,"RangInterest"],sep=" ")}
        if (totalSmallASV[i,"RangInterestMix"]=="Alveolata_X" ) { totalSmallASV[i,"RangInterestMix"] <- paste("Unaffiliated",totalSmallASV[i,"RangInterest"],sep=" ")}}
      totalSmallASV<-as.data.frame(totalSmallASV) %>% select(-RangInterest)
      totalSmallASV$Count <- totalSmallASV$TotalSmall*100/sum(totalSmallASV$TotalSmall)
      totalSmallASV <- totalSmallASV %>%
        arrange(desc(RangInterestMix)) %>%
        mutate(lab.ypos = cumsum(Count) - 0.5*Count)
      totalSmallASV$label <- paste(round(totalSmallASV$Count,1), "%", sep = "")
      for (i in rownames(totalSmallASV)) {
        if (totalSmallASV[i,"label"] == "0%") { totalSmallASV[i,"label"] <- NA}}
      for (i in rownames(totalSmallASV)) {
        if (is.na(totalSmallASV[i,"label"]) == FALSE) { totalSmallASV[i,"label"] <- paste(totalSmallASV[i,"RangInterestMix"]," : ",totalSmallASV[i,"label"], sep = "")}}
      totalSmallASV$Fraction<- rep("Small", each = nrow(totalSmallASV))
    ## Large
      totalLargeASV <- data_asv_tax %>% select(ASV_Id,TotalLarge,RangInterest,RangInterestMix)
      row.names(totalLargeASV)<-totalLargeASV$ASV_Id ; totalLargeASV <- totalLargeASV %>% select(-ASV_Id)
      totalLargeASV <- totalLargeASV %>% group_by(RangInterest,RangInterestMix) %>% summarise_all(sum)
      for (i in rownames(totalLargeASV)) { if (totalLargeASV[i,"RangInterestMix"]=="Not Affiliated" && totalLargeASV[i,"RangInterest"]!="Not Affiliated") { totalLargeASV[i,"RangInterestMix"] <- paste("Unaffiliated",totalLargeASV[i,"RangInterest"],sep=" ")}
        if (totalLargeASV[i,"RangInterestMix"]=="Alveolata_X" ) { totalLargeASV[i,"RangInterestMix"] <- paste("Unaffiliated",totalLargeASV[i,"RangInterest"],sep=" ")}}
      totalLargeASV<-as.data.frame(totalLargeASV) %>% select(-RangInterest)
      totalLargeASV$Count <- totalLargeASV$TotalLarge*100/sum(totalLargeASV$TotalLarge)
      totalLargeASV <- totalLargeASV %>%
        arrange(desc(RangInterestMix)) %>%
        mutate(lab.ypos = cumsum(Count) - 0.5*Count)
      totalLargeASV$label <- paste(round(totalLargeASV$Count,1), "%", sep = "")
      for (i in rownames(totalLargeASV)) {
        if (totalLargeASV[i,"label"] == "0%") { totalLargeASV[i,"label"] <- NA}}
      for (i in rownames(totalLargeASV)) {
        if (is.na(totalLargeASV[i,"label"]) == FALSE) { totalLargeASV[i,"label"] <- paste(totalLargeASV[i,"RangInterestMix"]," : ",totalLargeASV[i,"label"], sep = "")}}
      totalLargeASV$Fraction<- rep("Large", each = nrow(totalLargeASV))
    ## Fraction
      colnames(totalSmallASV)[2]  <- "value"
      totalSmallASV$Sum <- rep(0, each = nrow(totalSmallASV))
      totalSmallASV$Sum <- sum(totalSmallASV$value)
      colnames(totalLargeASV)[2]  <- "value"
      totalLargeASV$Sum <- rep(0, each = nrow(totalLargeASV))
      totalLargeASV$Sum <- sum(totalLargeASV$value)
      totalFractionASV <- rbind(totalSmallASV,totalLargeASV)
    ## plot
      ky <- ggplot(totalFractionASV, mapping = aes(y= Count, x = Fraction, fill = factor(RangInterestMix,level=Orderspe), group = factor(RangInterestMix,level=Orderspe)), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity",color="black") +
        scale_fill_manual(values = paletspe) + theme_unique_art() +
        geom_label(aes(y = 106,label = paste(Sum,"ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
      ky <- ky + labs(x="Fraction",y="ASVs (%)") + theme(legend.position = "none") + theme(axis.title.y = element_blank())
      print(ky)
#
    # Period -------------------------------------------------------------------
    ## Avril
      total04ASV <- data_asv_tax %>% select(ASV_Id,Total04,RangInterest,RangInterestMix)
      row.names(total04ASV)<-total04ASV$ASV_Id ; total04ASV <- total04ASV %>% select(-ASV_Id)
      total04ASV <- total04ASV %>% group_by(RangInterest,RangInterestMix) %>% summarise_all(sum)
      for (i in rownames(total04ASV)) { if (total04ASV[i,"RangInterestMix"]=="Not Affiliated" && total04ASV[i,"RangInterest"]!="Not Affiliated") { total04ASV[i,"RangInterestMix"] <- paste("Unaffiliated",total04ASV[i,"RangInterest"],sep=" ")}
        if (total04ASV[i,"RangInterestMix"]=="Alveolata_X" ) { total04ASV[i,"RangInterestMix"] <- paste("Unaffiliated",total04ASV[i,"RangInterest"],sep=" ")}}
      total04ASV<-as.data.frame(total04ASV) %>% select(-RangInterest)
      total04ASV$Count <- total04ASV$Total04*100/sum(total04ASV$Total04)
      total04ASV <- total04ASV %>%
        arrange(desc(RangInterestMix)) %>%
        mutate(lab.ypos = cumsum(Count) - 0.5*Count)
      total04ASV$label <- paste(round(total04ASV$Count,1), "%", sep = "")
      for (i in rownames(total04ASV)) {
        if (total04ASV[i,"label"] == "0%") { total04ASV[i,"label"] <- NA}}
      for (i in rownames(total04ASV)) {
        if (is.na(total04ASV[i,"label"]) == FALSE) { total04ASV[i,"label"] <- paste(total04ASV[i,"RangInterestMix"]," : ",total04ASV[i,"label"], sep = "")}}
      total04ASV$Period<- rep("04", each = nrow(total04ASV))
    ## Juin
      total06ASV <- data_asv_tax %>% select(ASV_Id,Total06,RangInterest,RangInterestMix)
      row.names(total06ASV)<-total06ASV$ASV_Id ; total06ASV <- total06ASV %>% select(-ASV_Id)
      total06ASV <- total06ASV %>% group_by(RangInterest,RangInterestMix) %>% summarise_all(sum)
      for (i in rownames(total06ASV)) { if (total06ASV[i,"RangInterestMix"]=="Not Affiliated" && total06ASV[i,"RangInterest"]!="Not Affiliated") { total06ASV[i,"RangInterestMix"] <- paste("Unaffiliated",total06ASV[i,"RangInterest"],sep=" ")}
        if (total06ASV[i,"RangInterestMix"]=="Alveolata_X" ) { total06ASV[i,"RangInterestMix"] <- paste("Unaffiliated",total06ASV[i,"RangInterest"],sep=" ")}}
      total06ASV<-as.data.frame(total06ASV) %>% select(-RangInterest)
      total06ASV$Count <- total06ASV$Total06*100/sum(total06ASV$Total06)
      total06ASV <- total06ASV %>%
        arrange(desc(RangInterestMix)) %>%
        mutate(lab.ypos = cumsum(Count) - 0.5*Count)
      total06ASV$label <- paste(round(total06ASV$Count,1), "%", sep = "")
      for (i in rownames(total06ASV)) {
        if (total06ASV[i,"label"] == "0%") { total06ASV[i,"label"] <- NA}}
      for (i in rownames(total06ASV)) {
        if (is.na(total06ASV[i,"label"]) == FALSE) { total06ASV[i,"label"] <- paste(total06ASV[i,"RangInterestMix"]," : ",total06ASV[i,"label"], sep = "")}}
      total06ASV$Period<- rep("06", each = nrow(total06ASV))
    ## Septembre
      total09ASV <- data_asv_tax %>% select(ASV_Id,Total09,RangInterest,RangInterestMix)
      row.names(total09ASV)<-total09ASV$ASV_Id ; total09ASV <- total09ASV %>% select(-ASV_Id)
      total09ASV <- total09ASV %>% group_by(RangInterest,RangInterestMix) %>% summarise_all(sum)
      for (i in rownames(total09ASV)) { if (total09ASV[i,"RangInterestMix"]=="Not Affiliated" && total09ASV[i,"RangInterest"]!="Not Affiliated") { total09ASV[i,"RangInterestMix"] <- paste("Unaffiliated",total09ASV[i,"RangInterest"],sep=" ")}
        if (total09ASV[i,"RangInterestMix"]=="Alveolata_X" ) { total09ASV[i,"RangInterestMix"] <- paste("Unaffiliated",total09ASV[i,"RangInterest"],sep=" ")}}
      total09ASV<-as.data.frame(total09ASV) %>% select(-RangInterest)
      total09ASV$Count <- total09ASV$Total09*100/sum(total09ASV$Total09)
      total09ASV <- total09ASV %>%
        arrange(desc(RangInterestMix)) %>%
        mutate(lab.ypos = cumsum(Count) - 0.5*Count)
      total09ASV$label <- paste(round(total09ASV$Count,1), "%", sep = "")
      for (i in rownames(total09ASV)) {
        if (total09ASV[i,"label"] == "0%") { total09ASV[i,"label"] <- NA}}
      for (i in rownames(total09ASV)) {
        if (is.na(total09ASV[i,"label"]) == FALSE) { total09ASV[i,"label"] <- paste(total09ASV[i,"RangInterestMix"]," : ",total09ASV[i,"label"], sep = "")}}
      total09ASV$Period<- rep("09", each = nrow(total09ASV))
    ## Octobre
      total11ASV <- data_asv_tax %>% select(ASV_Id,Total11,RangInterest,RangInterestMix)
      row.names(total11ASV)<-total11ASV$ASV_Id ; total11ASV <- total11ASV %>% select(-ASV_Id)
      total11ASV <- total11ASV %>% group_by(RangInterest,RangInterestMix) %>% summarise_all(sum)
      for (i in rownames(total11ASV)) { if (total11ASV[i,"RangInterestMix"]=="Not Affiliated" && total11ASV[i,"RangInterest"]!="Not Affiliated") { total11ASV[i,"RangInterestMix"] <- paste("Unaffiliated",total11ASV[i,"RangInterest"],sep=" ")}
        if (total11ASV[i,"RangInterestMix"]=="Alveolata_X" ) { total11ASV[i,"RangInterestMix"] <- paste("Unaffiliated",total11ASV[i,"RangInterest"],sep=" ")}}
      total11ASV<-as.data.frame(total11ASV) %>% select(-RangInterest)
      total11ASV$Count <- total11ASV$Total11*100/sum(total11ASV$Total11)
      total11ASV <- total11ASV %>%
        arrange(desc(RangInterestMix)) %>%
        mutate(lab.ypos = cumsum(Count) - 0.5*Count)
      total11ASV$label <- paste(round(total11ASV$Count,1), "%", sep = "")
      for (i in rownames(total11ASV)) {
        if (total11ASV[i,"label"] == "0%") { total11ASV[i,"label"] <- NA}}
      for (i in rownames(total11ASV)) {
        if (is.na(total11ASV[i,"label"]) == FALSE) { total11ASV[i,"label"] <- paste(total11ASV[i,"RangInterestMix"]," : ",total11ASV[i,"label"], sep = "")}}
      total11ASV$Period<- rep("11", each = nrow(total11ASV))
    ## Periods
      colnames(total04ASV)[2]  <- "value"
      total04ASV$Sum <- rep(0, each = nrow(total04ASV))
      total04ASV$Sum <- sum(total04ASV$value)
      colnames(total06ASV)[2]  <- "value"
      total06ASV$Sum <- rep(0, each = nrow(total06ASV))
      total06ASV$Sum <- sum(total06ASV$value)
      colnames(total09ASV)[2]  <- "value"
      total09ASV$Sum <- rep(0, each = nrow(total09ASV))
      total09ASV$Sum <- sum(total09ASV$value)
      colnames(total11ASV)[2]  <- "value"
      total11ASV$Sum <- rep(0, each = nrow(total11ASV))
      total11ASV$Sum <- sum(total11ASV$value)
      totalPeriodASV <- rbind(total04ASV,total06ASV,total09ASV,total11ASV)
    #Figure
      ly <- ggplot(totalPeriodASV, mapping = aes(y= Count, x = Period, fill = factor(RangInterestMix,level=Orderspe), group = factor(RangInterestMix,level=Orderspe)), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity",color="black") +
        scale_fill_manual(values = paletspe) + theme_unique_art() +
        geom_label(aes(y = 106,label = paste(Sum,"ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
      ly <- ly + labs(x="Period",y="ASVs (%)") + theme(legend.position = "none") + theme(axis.title.y = element_blank())
      print(ly) 
#
    # Coplot -------------------------------------------------------------------
      svglite("Hist-Taxomy/ASVMIX-Total.svg",width = 14.00,height = 10.00)
      b_plot <- plot_grid(iy,jy,ky,ly, legendASV, ncol = 5, nrow = 1, rel_widths = c(3,3,3,6,3),rel_heights = c(3))
      print(b_plot)
      dev.off()
#
  # Table Mix -------------------------------------------------------------------
  ## Sequence
    write.table(totalCycleSequence, file = "dataTables/TotalMix_Cycles_Sequence.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(totalZoneSequence, file = "dataTables/TotalMix_Zone_Sequence.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(totalFractionSequence, file = "dataTables/TotalMix_Fraction_Sequence.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(totalPeriodSequence, file = "dataTables/TotalMix_Periods_Sequence.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  ## ASV
    write.table(totalCycleASV, file = "dataTables/TotalMix_Cycles_ASV.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(totalZoneASV, file = "dataTables/TotalMix_Zone_ASV.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(totalFractionASV, file = "dataTables/TotalMix_Fraction_ASV.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(totalPeriodASV, file = "dataTables/TotalMix_Periods_ASV.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
#
# Balance plot ------------------------------------------------------------
  # Sequence ---------------------------------------------------------------------
  ## Mixolimnion
    balance_totalZoneSequence_Mixolimnion <- totalZoneSequence %>% filter(Zone == "Mixolimnion") %>% select(RangInterestMix,value,Zone) %>% group_by(RangInterestMix,Zone) %>% summarise_all(sum)
    balance_totalZoneSequence_Mixolimnion <- balance_totalZoneSequence_Mixolimnion %>% arrange(totalZoneSequence[1:nrow(balance_totalZoneSequence_Mixolimnion),"RangInterestMix"])
  ## Monimolimnion
    balance_totalZoneSequence_Monimolimnion <- totalZoneSequence %>% filter(Zone == "Monimolimnion") %>% select(RangInterestMix,value,Zone) %>% group_by(RangInterestMix,Zone) %>% summarise_all(sum)
    balance_totalZoneSequence_Monimolimnion <- balance_totalZoneSequence_Monimolimnion %>% arrange(totalZoneSequence[1:nrow(balance_totalZoneSequence_Monimolimnion),"RangInterestMix"])
  ## Bind
    balance_totalZoneSequence<-balance_totalZoneSequence_Mixolimnion
    for (i in row.names(balance_totalZoneSequence)) { balance_totalZoneSequence[i,"value"] <- foldchange(balance_totalZoneSequence_Mixolimnion[i,"value"],balance_totalZoneSequence_Monimolimnion[i,"value"])
    balance_totalZoneSequence[i,"value"] <- foldchange2logratio(balance_totalZoneSequence[i,"value"][[1]])
    if (balance_totalZoneSequence_Mixolimnion[i,"value"] < balance_totalZoneSequence_Monimolimnion[i,"value"]) {balance_totalZoneSequence[i,"Zone"] <- "Monimolimnion"}
    if (balance_totalZoneSequence_Mixolimnion[i,"value"] > balance_totalZoneSequence_Monimolimnion[i,"value"]) {balance_totalZoneSequence[i,"Zone"] <- "Mixolimnion"}
    balance_totalZoneSequence[i,"Label"] <- paste(max(balance_totalZoneSequence_Mixolimnion[i,"value"],balance_totalZoneSequence_Monimolimnion[i,"value"]),"VS",min(balance_totalZoneSequence_Mixolimnion[i,"value"],balance_totalZoneSequence_Monimolimnion[i,"value"]),sep=" ")
    balance_totalZoneSequence[i,"Abundance"] <- max(balance_totalZoneSequence_Mixolimnion[i,"value"],balance_totalZoneSequence_Monimolimnion[i,"value"])}
  ## Plot
    svglite("Hist-Taxomy/Sequence-Balance-Zone.svg",width = 8.00,height = 8.00)
    ggplot(balance_totalZoneSequence, aes(x = factor(RangInterestMix,level=rev(Orderspe)), y = value)) +
      geom_bar(aes(fill = Zone,alpha=log10(Abundance)), stat = 'identity',color="black") +
      coord_flip() +  # horizontal bars
      geom_text(aes(label = Label ,y = 0,vjust = ifelse(value >= 0, 0.5, 0.5), hjust = ifelse(value>= 0, 1.2, -0.2)),color = "black",size = 3,show.legend = FALSE, fill = "white") +
      labs(x="",y="Sequences abundance log-ratio") +
      theme(axis.text.y = element_blank()) + theme_unique_art() + scale_fill_paletteer_d("ggthemes::Nuriel_Stone") + guides(alpha = FALSE)
    dev.off()
#
  # ASV ---------------------------------------------------------------------
  ## Mixolimnion
    balance_totalZoneASV_Mixolimnion <- totalZoneASV %>% filter(Zone == "Mixolimnion") %>% select(RangInterestMix,value,Zone) %>% group_by(RangInterestMix,Zone) %>% summarise_all(sum)
    balance_totalZoneASV_Mixolimnion <- balance_totalZoneASV_Mixolimnion %>% arrange(totalZoneASV[1:nrow(balance_totalZoneASV_Mixolimnion),"RangInterestMix"])
  ## Monimolimnion
    balance_totalZoneASV_Monimolimnion <- totalZoneASV %>% filter(Zone == "Monimolimnion") %>% select(RangInterestMix,value,Zone) %>% group_by(RangInterestMix,Zone) %>% summarise_all(sum)
    balance_totalZoneASV_Monimolimnion <- balance_totalZoneASV_Monimolimnion %>% arrange(totalZoneASV[1:nrow(balance_totalZoneASV_Monimolimnion),"RangInterestMix"])
  ## Bind
    balance_totalZoneASV<-balance_totalZoneASV_Mixolimnion
    for (i in row.names(balance_totalZoneASV)) { balance_totalZoneASV[i,"value"] <- foldchange(balance_totalZoneASV_Mixolimnion[i,"value"],balance_totalZoneASV_Monimolimnion[i,"value"])
    balance_totalZoneASV[i,"value"] <- foldchange2logratio(balance_totalZoneASV[i,"value"][[1]])
    if (balance_totalZoneASV_Mixolimnion[i,"value"] < balance_totalZoneASV_Monimolimnion[i,"value"]) {balance_totalZoneASV[i,"Zone"] <- "Monimolimnion"}
    if (balance_totalZoneASV_Mixolimnion[i,"value"] > balance_totalZoneASV_Monimolimnion[i,"value"]) {balance_totalZoneASV[i,"Zone"] <- "Mixolimnion"}
    balance_totalZoneASV[i,"Abundance"] <- max(balance_totalZoneASV_Mixolimnion[i,"value"],balance_totalZoneASV_Monimolimnion[i,"value"])}
  ## Plot
    svglite("Hist-Taxomy/ASV-Balance-Zone.svg",width = 8.00,height = 8.00)
    ggplot(balance_totalZoneASV, aes(x = factor(RangInterestMix,level=rev(Orderspe)), y = value)) +
      geom_bar(aes(fill = Zone,alpha=log10(Abundance)), stat = 'identity',color="black") +
      coord_flip() +  # horizontal bars
      geom_text(aes(label = round(Abundance) ,y = 0,vjust = ifelse(value >= 0, 0.5, 0.5), hjust = ifelse(value>= 0, 1.2, -0.2)),color = "black",size = 3,show.legend = FALSE, fill = "white") +
      labs(x="",y="ASVs log-ratio") +
      theme(axis.text.y = element_blank()) + theme_unique_art() + scale_fill_paletteer_d("ggthemes::Nuriel_Stone") + guides(alpha = FALSE)
    dev.off()
#
# Save Rdata --------------------------------------------------------------
  save.image(file = "dataTables/my_work_space_Compo.RData")
#