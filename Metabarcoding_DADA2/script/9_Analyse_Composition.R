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
input <- "V4-unified-out"
region <- "V4"
sortop <- "OnlyOne"
Mode <- "Superphylum"
Group <- "Eukaryota"
RarefyYoN <- "yes"
UnifyYoN <- "no"
}
#
# Set argument if using R
args = commandArgs(trailingOnly=TRUE)
if ( inputmode == FALSE ) {
  input <- args[1]
  region <- args[2]
  sortop <- args[3]
  Mode <- args[4]
  Group <- args[5]
  RarefyYoN <- args[6]
  UnifyYoN <- args[7]
}
#
# Import package and palette -----------------------------------------------------------
pkg <- c("ggplot2", "readxl","dplyr","tidyr","cowplot","FactoMineR","factoextra","reshape2","varhandle","ggrepel","ggpubr","ggsci","scales","hrbrthemes","GUniFrac","svglite","treemap", "VennDiagram","stringr")
lapply(pkg, require, character.only = TRUE)
palette <- c("#AD002ACC","#EEA236CC","#00468BCC","#0099B4CC","#D62768CC","#ADB6B6CC","#42B540CC","#357EBDCC","#1B1919CC","#ED0000CC","#FF7F0ECC","#925E9FCC","#D43F3ACC","#9632B8CC","#46B8DACC","#5CB85CCC","#EFC008CC")
  #show_col(palette)

# Set directory, create file result -----------------------------------------------------------
if (inputmode == TRUE) {
  current <- dirname(rstudioapi::getSourceEditorContext()$path)
  setwd(current)
  input <- paste(current,"../dataDADA2/result",input, sep = "/")
}
if (inputmode == FALSE) {
  input <- paste("../dataDADA2/result",input, sep = "/")
}
# Set working directory
if (dir.exists(input) == TRUE) { setwd(input) }
# Create result files
system("mkdir Stat-Analyse")
system("mkdir Composition")
system("mkdir Hist-Taxomy")
system("mkdir AFC-Distribution")
system("mkdir Biplot")
system("mkdir Venn")
system("mkdir ASV-Table")
system("mkdir Krona")

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


# Input ASV Table ---------------------------------------------------------
tableVinput <- read.csv(file = paste("Table/ASV_table.csv",sep = "/"), sep = "\t")
##0.0005% filter
if (sortop == "Bokulich") {
  amplicon <- grep(pattern = "OSTA", colnames(tableVinput), value = TRUE)
  tableVinput$SUM <- rowSums(tableVinput %>% select(all_of(amplicon)))
  tableVinput <- tableVinput %>% filter(SUM > 0.0005*sum(tableVinput$SUM)/100)
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

# Prepare data inf --------------------------------------------------------
infdataini <- read.table(file = "../../../rawdata/data-inf.txt", sep = "\t", header = FALSE,row.names = "row", col.names = c("row","Id","Condition","Technologie","Region"))
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
  # Prepare seq_mat and Pool (seq_mat and seq_mat_pool) ---------------------------------------------------------------------
## Prepare asv_mat
amplicon <- c(grep(pattern = "OSTA", colnames(tableVinput), value = TRUE))
seq_mat <- tableVinput %>% select(all_of(amplicon))
seq_mat[,all_of(amplicon)] <- lapply(seq_mat[,all_of(amplicon)], as.numeric)
seq_mat$ASV_Id <- row.names(seq_mat)

## Prepare correspondance table
tblc <- samples_df %>% select("Amplicon", "Date", "Replicat", "Condition")
tblc$Condition_Date <- paste(tblc$Condition,tblc$Date, sep = "_")
tblc <- tblc %>% select(-"Date",-"Condition")
tblc1 <- tblc %>% filter(tblc$Replicat == 1)
tblc2 <- tblc %>% filter(tblc$Replicat == 2)
tblcx <- merge(x = tblc1,y = tblc2, by = "Condition_Date", all = TRUE)
Amplicon_without_duplicat <- c(tblcx[,"Amplicon.x"][which(is.na(tblcx[,"Amplicon.y"]) == TRUE)])
for (i in row.names(tblcx)) { if (is.na(tblcx[i,"Amplicon.y"]) == TRUE) { tblcx[i,"Amplicon.y"] <- tblcx[i,"Amplicon.x"]}}
for (i in row.names(tblcx)) { if (is.na(tblcx[i,"Replicat.y"]) == TRUE) { tblcx[i,"Replicat.y"] <- tblcx[i,"Replicat.x"]}}

## Unify
if (UnifyYoN == "yes") {
seq_mat_pool <- as.data.frame(x=seq_mat[,"ASV_Id"],row.names = row.names(seq_mat))
colnames(seq_mat_pool) <- "ASV_Id"
for (h in tblcx[,"Amplicon.x"] ) {
  f <- (tblcx %>% filter(Amplicon.x == h) %>% select(Amplicon.y))[[1]]
  print (h)
  print (f)
  for (r in row.names(seq_mat)) {
    if (seq_mat[r,h] * seq_mat[r,f] == 0) { seq_mat_pool[r,h] <- round((seq_mat[r,h]+seq_mat[r,f]))}
    if (seq_mat[r,h] * seq_mat[r,f] != 0) { seq_mat_pool[r,h] <- round((seq_mat[r,h]+seq_mat[r,f])/2)}
  }}
}
if (UnifyYoN == "no") { 
  seq_mat_pool <- seq_mat 
  samples_df <- samples_df %>% filter(Replicat == 1)
}

  # Rarefy (seq_mat_pool_rare) ---------------------------------------------------------------------
##Rarefy yes or no
if (RarefyYoN == "yes") {
seq_mat_poolt <- seq_mat_pool %>% select(-"ASV_Id")
seq_mat_poolt <- t(seq_mat_poolt)
seq_mat_pool_rare <- as.data.frame(t(Rarefy(seq_mat_poolt)$otu.tab.rff)) ## Ref : Jun Chen et al. (2012). Associating microbiome composition with environmental covariates using generalized UniFrac distances. 28(16): 2106–2113.
rm(seq_mat_poolt)
seq_mat_pool_rare$ASV_Id <- rownames(seq_mat_pool_rare)
#for (i in row.names(seq_mat_pool_rare)) { seq_mat_pool_rare[i,"Taxonomy"] <- seq_mat_pool[i,"Taxonomy"]}
}
if (RarefyYoN == "no") { seq_mat_pool_rare <- seq_mat_pool }

  # Prepare asv_mat (asv_mat, asv_mat_pool and asv_mat_pool_rare) ------------------------------------------------------------------
## asv_mat
asv_mat <- seq_mat
amplicon <- c(grep(pattern = "OSTA", colnames(seq_mat), value = TRUE))
asv_mat[,all_of(amplicon)][asv_mat[,all_of(amplicon)] != 0] <- 1
## asv_mat_pool
asv_mat_pool <- seq_mat_pool
amplicon <- c(grep(pattern = "OSTA", colnames(asv_mat_pool), value = TRUE))
asv_mat_pool[,all_of(amplicon)][asv_mat_pool[,all_of(amplicon)] != 0] <- 1
## asv_mat_pool_rare
asv_mat_pool_rare <- seq_mat_pool_rare
amplicon <- c(grep(pattern = "OSTA", colnames(asv_mat_pool_rare), value = TRUE))
asv_mat_pool_rare[,all_of(amplicon)][asv_mat_pool_rare[,all_of(amplicon)] != 0] <- 1

# Stat Rarefy & Pool -------------------------------------------------------------
## Sequence
amplicon <- c(grep(pattern = "OSTA", colnames(seq_mat_pool), value = TRUE))
avRarefyS <- as.data.frame(colSums(seq_mat_pool[,all_of(amplicon)]))
amplicon <- c(grep(pattern = "OSTA", colnames(seq_mat_pool_rare), value = TRUE))
apRarefyS <- as.data.frame(colSums(seq_mat_pool_rare[,all_of(amplicon)]))
statRarefy <- cbind(avRarefyS,apRarefyS)
## ASV
amplicon <- c(grep(pattern = "OSTA", colnames(asv_mat_pool), value = TRUE))
avRarefyO <- as.data.frame(colSums(asv_mat_pool[,all_of(amplicon)]))
amplicon <- c(grep(pattern = "OSTA", colnames(asv_mat_pool_rare), value = TRUE))
apRarefyO <- as.data.frame(colSums(asv_mat_pool_rare[,all_of(amplicon)]))
statRarefy <- cbind(statRarefy,avRarefyO,apRarefyO)
## Totaux
statRarefy["Total",]<-colSums(statRarefy)
colnames(statRarefy) <- c("avRarefy-Sequence","apRarefy-Sequence","avRarefy-ASV","apRarefy-ASV")
## Remove empty ASVs
### Avant raréfaction
ASVav <- asv_mat_pool
ASVav$ASV_Id <- row.names(ASVav)
amplicon <- c(grep(pattern = "OSTA", colnames(ASVav), value = TRUE))
for ( i in row.names(ASVav)) { if (sum(ASVav[i,all_of(amplicon)]) == 0) { ASVav[i,"ASV_Id"] <- "Uniq"}}
ASVav <- ASVav %>% filter(ASV_Id != "Uniq")
### Après rarefaction
ASVap <- asv_mat_pool_rare
ASVap$ASV_Id <- row.names(ASVap)
amplicon <- c(grep(pattern = "OSTA", colnames(ASVap), value = TRUE))
for ( i in row.names(ASVap)) { if (sum(ASVap[i,all_of(amplicon)]) == 0) { ASVap[i,"ASV_Id"] <- "Uniq"}}
ASVap <- ASVap %>% filter(ASV_Id != "Uniq")
### Final
statRarefy["Total","avRarefy-ASV"] <- nrow(ASVav)
statRarefy["Total","apRarefy-ASV"] <- nrow(ASVap)
write.table(statRarefy, file = "Table/StatRarefy_withoutDuplicat.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

  # Sum files ---------------------------------------------------------------
## Rarefaction stat Fig
statRarefyASV <- statRarefy %>% select("avRarefy-ASV","apRarefy-ASV")
colnames(statRarefyASV) <- c("NON","OUI")
statRarefyASV$Echantillon <- row.names(statRarefyASV)
amplicon <- c(grep(pattern = "OSTA", colnames(ASVap), value = TRUE))
statRarefyASV <- melt(statRarefyASV[all_of(amplicon),], id = "Echantillon")
my_comp <- list(c("OUI","NON"))

svglite("Stat-Analyse/Stat-Rarefy.svg",width = 3.00,height = 4.00)
R <- ggplot(statRarefyASV, aes(y = value, x = variable)) + geom_boxplot() + geom_point(aes(color = Echantillon), size = 2.5) + 
  stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
  stat_compare_means(method="wilcox.test", label.y = max(statRarefyASV$value)+0.12*max(statRarefyASV$value)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  labs(x="Raréfaction",y="Nombre d'ASVs") #+ guides(color = FALSE)
R <- R + theme(legend.position="none")
print(R)
dev.off() 

## Pooling stat Fig
statRarefy <- statRarefy[!(row.names(statRarefy) %in% "Total"), ]
statRarefy[,"Fusion"] <- "OUI"
statRarefy[Amplicon_without_duplicat,"Fusion"] <- "NON"
statRarefy$Echantillon <- rownames(statRarefy)
my_comp <- list(c("OUI","NON"))
### ASV après Rarefy
svglite("Stat-Analyse/Analyse-SumApRare-ASV.svg",width = 3.00,height = 4.00)
J <- ggplot(statRarefy, aes(y = `apRarefy-ASV`, x = Fusion)) + geom_boxplot() + geom_point(aes(color = Echantillon), size = 2.5) + 
  stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
  stat_compare_means(method="wilcox.test", label.y = max(statRarefy$`apRarefy-ASV`)+0.12*max(statRarefy$`apRarefy-ASV`)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  labs(x="Regroupement",y="Nombre d'ASVs") #+ guides(color = FALSE)
J <- J + theme(legend.position="none")
print(J)
dev.off()
### Séquence avant Rarefy
svglite("Stat-Analyse/Analyse-SumAvRare-Sequences.svg",width = 3.00,height = 4.00)
G <- ggplot(statRarefy, aes(y = `avRarefy-Sequence`, x = Fusion)) + geom_boxplot() + geom_point(aes(color = Echantillon), size = 2.5) + 
  stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
  stat_compare_means(method="wilcox.test", label.y = max(statRarefy$`avRarefy-Sequence`)+0.12*max(statRarefy$`avRarefy-Sequence`)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  labs(x="Regroupement",y="Nombre de séquences") #+ guides(color = FALSE)
G <- G + theme(legend.position="none")
print(G)
dev.off()  
### ASV avant Rarefy
svglite("Stat-Analyse/Analyse-SumAvRare-ASV.svg",width = 3.00,height = 4.00)
H <- ggplot(statRarefy, aes(y = `avRarefy-ASV`, x = Fusion)) + geom_boxplot() + geom_point(aes(color = Echantillon), size = 2.5) + 
  stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
  stat_compare_means(method="wilcox.test", label.y = max(statRarefy$`avRarefy-ASV`)+0.12*max(statRarefy$`avRarefy-ASV`)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  labs(x="Regroupement",y="Nombre d'ASVs") #+ guides(color = FALSE)
H <- H + theme(legend.position="none")
print(H)
dev.off()
  # Taxonomy & statistics -----------------------------------------------------
## Create Tax stat table
tax_table <- tableVinput %>% select(Taxonomy)
tax_table$ASV_Id <- row.names(tax_table)
tax_table <- separate(tax_table, Taxonomy, c("Domain","Superphylum","Phylum","Class","Order","Family","Genus","Species"), sep =";")
tax_table[tax_table == ""] <- NA
tax_table[tax_table == " "] <- NA
tax_table[tax_table == "NA"] <- NA
## Continue Treatment tax_table
for (i in c("Domain","Superphylum","Phylum","Class","Order","Family","Genus")) {
  j <- grep(i, colnames(tax_table)) + 1
  for (f in rownames(tax_table)) { 
    if (is.na(tax_table[f,i]) == TRUE) { tax_table[f,j] <- NA }}}

tax_tablemix <- tax_table
tax_tablemix$Taxonomy <- paste(tax_table[,"Domain"],tax_table[,"Superphylum"],tax_table[,"Phylum"],tax_table[,"Class"],tax_table[,"Order"],tax_table[,"Family"],tax_table[,"Genus"],tax_table[,"Species"], sep = ";")
tax_tablemix <- tax_tablemix %>% select(ASV_Id, Taxonomy)

## Add taxonomy to seq and asv matrix
seq_mat <- merge(seq_mat,tax_tablemix,by = "ASV_Id")
rownames(seq_mat) <- seq_mat$ASV_Id ; seq_mat <- seq_mat %>% select(-ASV_Id)
seq_mat_pool <- merge(seq_mat_pool,tax_tablemix,by = "ASV_Id")
rownames(seq_mat_pool) <- seq_mat_pool$ASV_Id ; seq_mat_pool <- seq_mat_pool %>% select(-ASV_Id)
seq_mat_pool_rare <- merge(seq_mat_pool_rare,tax_tablemix,by = "ASV_Id")
rownames(seq_mat_pool_rare) <- seq_mat_pool_rare$ASV_Id ; seq_mat_pool_rare <- seq_mat_pool_rare %>% select(-ASV_Id)
asv_mat <- merge(asv_mat,tax_tablemix,by = "ASV_Id")
rownames(asv_mat) <- asv_mat$ASV_Id ; asv_mat <- asv_mat %>% select(-ASV_Id)
asv_mat_pool <- merge(asv_mat_pool,tax_tablemix,by = "ASV_Id")
rownames(asv_mat_pool) <- asv_mat_pool$ASV_Id ; asv_mat_pool <- asv_mat_pool %>% select(-ASV_Id)
asv_mat_pool_rare <- merge(asv_mat_pool_rare,tax_tablemix,by = "ASV_Id")
rownames(asv_mat_pool_rare) <- asv_mat_pool_rare$ASV_Id ; asv_mat_pool_rare <- asv_mat_pool_rare %>% select(-ASV_Id)


tax_table_TF <- as.data.frame(is.na(tax_table))
for (level in c("Domain","Superphylum","Phylum","Class","Order","Family","Genus","Species")) {
  assign(paste("count", level, sep = "_"), nrow(as.data.frame(tax_table_TF[,level][tax_table_TF[,level] == FALSE])))}
tax_stat <- as.data.frame(c("Domain","Superphylum","Phylum","Class","Order","Family","Genus","Species")) ; colnames(tax_stat) <- "Level"
tax_stat$Count <- c(count_Domain,count_Superphylum,count_Phylum,count_Class,count_Order,count_Family,count_Genus,count_Species)

tax_stat$`Not affiliated` <- nrow(tableVinput) - tax_stat$Count
tax_stat <- melt(tax_stat, id = "Level")
## Plot
tax_stat$Level <- factor(tax_stat$Level , levels=c("Domain","Superphylum","Phylum","Class","Order","Family","Genus","Species"))
tax_stat$variable <- factor(tax_stat$variable , levels=c("Not affiliated","Count"))

#
svglite("Stat-Analyse/Stat-Taxo.svg",width = 9.00,height = 5.00)
tax_plot <- ggplot(tax_stat, mapping = aes(x= Level, y = value, fill = variable, color = variable ,linetype = variable), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
  scale_fill_manual(values = alpha(c("#1212ff","#1212ff"), c(0.25,0.5))) +
  geom_label(aes(y = value + 0.07*max(value),label = paste(value," ASVs","\n", round(value*100/nrow(tableVinput),1)," %",sep =""), color = variable, alpha = variable),size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
  scale_color_manual(values = c("#0000FF00","black")) +
  scale_alpha_manual(values = c(0,0.5)) +
  scale_linetype_manual(values = c("blank","dashed")) +
  guides(fill=guide_legend("Affiliation"), color = FALSE, linetype = FALSE) +
  labs(x="Level",y="ASVs (%)") +
  theme(legend.title = element_text(face="bold"),
        axis.title = element_text(color = "black", face = "bold"))
print(tax_plot)
dev.off()

  # ASV distribution statistics ---------------------------------------------
## Raw
amplicon <- c(grep(pattern = "OSTA", colnames(seq_mat), value = TRUE))
ASV_stat_raw <- as.data.frame(rowSums(seq_mat %>% select(all_of(amplicon)))) ; colnames(ASV_stat_raw) <- "ASV Abuncance"
ASV_stat_raw$type <- "Raw"
## Raw_Pool
amplicon <- c(grep(pattern = "OSTA", colnames(seq_mat_pool), value = TRUE))
ASV_stat_pool <- as.data.frame(rowSums(seq_mat_pool %>% select(all_of(amplicon)))) ; colnames(ASV_stat_pool) <- "ASV Abuncance"
ASV_stat_pool$type <- "Raw_Pool"
## Rarefy_Pool
amplicon <- c(grep(pattern = "OSTA", colnames(seq_mat_pool_rare), value = TRUE))
ASV_stat_pool_rare <- as.data.frame(rowSums(seq_mat_pool_rare %>% select(all_of(amplicon)))) ; colnames(ASV_stat_pool_rare) <- "ASV Abuncance"
ASV_stat_pool_rare$type <- "Rarefy_Pool"
## rbind
ASV_stat <- rbind(ASV_stat_raw,ASV_stat_pool,ASV_stat_pool_rare)
## x scale
w <- 0
i <- 1
while ( i < max(ASV_stat$`ASV Abuncance`)) { print (i) 
  w <- c(w,i)
  i <- i*4}
## plot
ASV_stat$type <- factor(ASV_stat$type , levels=c("Raw","Raw_Pool","Rarefy_Pool"))
svglite("Stat-Analyse/ASV-Stat.svg",width = 10.00,height = 5.00)
ggplot(ASV_stat, aes(x = `ASV Abuncance`)) +facet_grid(.~type) + 
  stat_bin(binwidth = 1,color = "white", fill = "#1212ff", alpha = 0.5) +
  scale_y_continuous() +
  scale_x_continuous(trans = log2_trans(),breaks = w) +
  labs(y="ASV #") +
  theme(legend.title = element_text(face="bold"),
        axis.title = element_text(color = "black", face = "bold"),
        strip.text.x = element_text(color = "black", face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
# Save ASV Tables ------------------------------------------------------------
#seq_mat_pool_rare$Total <- FALSE
#for (i in rownames(seq_mat_pool_rare)) {
#  if (rowSums(seq_mat_pool_rare %>% filter(rownames(seq_mat_pool_rare) == i) %>% select(all_of(amplicon))) == 0) { seq_mat_pool_rare[i,"Total"] <- TRUE }}

write.table(seq_mat, file = "ASV-Table/Seq-mat.csv", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
write.table(seq_mat_pool, file = "ASV-Table/Seq-mat-pool.csv", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
write.table(seq_mat_pool_rare, file = "ASV-Table/Seq-mat-pool-rare.csv", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
write.table(asv_mat, file = "ASV-Table/ASV-mat.csv", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
write.table(asv_mat_pool, file = "ASV-Table/ASV-mat-pool.csv", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
write.table(asv_mat_pool_rare, file = "ASV-Table/ASV-mat-pool-rare.csv", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)


# Create sort Condition pattern -------------------------------------------
#Cycle
CycleJour <- samples_df %>% filter(Cycle == "Jour") %>% filter(Replicat == 1)
CycleJour <- CycleJour$Amplicon
CycleNuit <- samples_df %>% filter(Cycle == "Nuit") %>% filter(Replicat == 1)
CycleNuit <- CycleNuit$Amplicon
#Zone
ZoneOxique <- samples_df %>% filter(`Zone` == "Oxique") %>% filter(Replicat == 1)
ZoneOxique <- ZoneOxique$Amplicon
ZoneAnoxique <- samples_df %>% filter(`Zone` == "Anoxique") %>% filter(Replicat == 1)
ZoneAnoxique <- ZoneAnoxique$Amplicon
#Fraction
FractionPetite <- samples_df %>% filter(`Fraction` == "Petite") %>% filter(Replicat == 1)
FractionPetite <- FractionPetite$Amplicon
FractionGrande <- samples_df %>% filter(`Fraction` == "Grande") %>% filter(Replicat == 1)
FractionGrande <- FractionGrande$Amplicon
#Date
Date04 <- samples_df %>% filter(Date == "04") %>% filter(Replicat == 1)
Date04 <- Date04$Amplicon
Date06 <- samples_df %>% filter(Date == "06") %>% filter(Replicat == 1)
Date06 <- Date06$Amplicon
Date09 <- samples_df %>% filter(Date == "09") %>% filter(Replicat == 1)
Date09 <- Date09$Amplicon
Date11 <- samples_df %>% filter(Date == "11") %>% filter(Replicat == 1)
Date11 <- Date11$Amplicon




# Create dataframe Condition --------------------------------------------
  # Séquence ----------------------------------------------------------------
#Cycle
##Jour
Jour_seq <- as.data.frame(rowSums(seq_mat_pool_rare %>% select(all_of(CycleJour))))
colnames(Jour_seq) <- "TotalJour"
Jour_seq$ASV_Id <- row.names(Jour_seq)
##Nuit
Nuit_seq <- as.data.frame(rowSums(seq_mat_pool_rare %>% select(all_of(CycleNuit))))
colnames(Nuit_seq) <- "TotalNuit"
Nuit_seq$ASV_Id <- row.names(Nuit_seq)
##Cycle
Cycle_seq <- merge(x = Jour_seq,y = Nuit_seq, by = "ASV_Id")

#Zone
##Oxique
Oxique_seq <- as.data.frame(rowSums(seq_mat_pool_rare %>% select(all_of(ZoneOxique))))
colnames(Oxique_seq) <- "TotalOxique"
Oxique_seq$ASV_Id <- row.names(Oxique_seq)
##Anoxique
Anoxique_seq <- as.data.frame(rowSums(seq_mat_pool_rare %>% select(all_of(ZoneAnoxique))))
colnames(Anoxique_seq) <- "TotalAnoxique"
Anoxique_seq$ASV_Id <- row.names(Anoxique_seq)
##Zone
Zone_seq <- merge(x = Oxique_seq,y = Anoxique_seq, by = "ASV_Id")

#Fraction
##Petite
Petite_seq <- as.data.frame(rowSums(seq_mat_pool_rare %>% select(all_of(FractionPetite))))
colnames(Petite_seq) <- "TotalPetite"
Petite_seq$ASV_Id <- row.names(Petite_seq)
##Grande
Grande_seq <- as.data.frame(rowSums(seq_mat_pool_rare %>% select(all_of(FractionGrande))))
colnames(Grande_seq) <- "TotalGrande"
Grande_seq$ASV_Id <- row.names(Grande_seq)
##Fraction
Fraction_seq <- merge(x = Petite_seq,y = Grande_seq, by = "ASV_Id")

#Date
##04
`04_seq` <- as.data.frame(rowSums(seq_mat_pool_rare %>% select(all_of(Date04))))
colnames(`04_seq`) <- "Total04"
`04_seq`$ASV_Id <- row.names(`04_seq`)
##06
`06_seq` <- as.data.frame(rowSums(seq_mat_pool_rare %>% select(all_of(Date06))))
colnames(`06_seq`) <- "Total06"
`06_seq`$ASV_Id <- row.names(`06_seq`)
##09
`09_seq` <- as.data.frame(rowSums(seq_mat_pool_rare %>% select(all_of(Date09))))
colnames(`09_seq`) <- "Total09"
`09_seq`$ASV_Id <- row.names(`09_seq`)
##11
`11_seq` <- as.data.frame(rowSums(seq_mat_pool_rare %>% select(all_of(Date11))))
colnames(`11_seq`) <- "Total11"
`11_seq`$ASV_Id <- row.names(`11_seq`)
##Date
Date_seq <- merge(x = `04_seq`,y = `06_seq`, by = "ASV_Id")
Date_seq <- merge(x = Date_seq,y = `09_seq`, by = "ASV_Id")
Date_seq <- merge(x = Date_seq,y = `11_seq`, by = "ASV_Id")

  # ASV ----------------------------------------------------------------
#Cycle
##Jour
Jour_asv <- as.data.frame(rowSums(asv_mat_pool_rare %>% select(all_of(CycleJour))))
colnames(Jour_asv) <- "TotalJour"
for (i in row.names(Jour_asv)) {if (Jour_asv[i,"TotalJour"] > 0) {Jour_asv[i,"TotalJour"] <- 1}}
Jour_asv$ASV_Id <- row.names(Jour_asv)
##Nuit
Nuit_asv <- as.data.frame(rowSums(asv_mat_pool_rare %>% select(all_of(CycleNuit))))
colnames(Nuit_asv) <- "TotalNuit"
for (i in row.names(Nuit_asv)) {if (Nuit_asv[i,"TotalNuit"] > 0) {Nuit_asv[i,"TotalNuit"] <- 1}}
Nuit_asv$ASV_Id <- row.names(Nuit_asv)
##Cycle
Cycle_asv <- merge(x = Jour_asv,y = Nuit_asv, by = "ASV_Id")

#Zone
##Oxique
Oxique_asv <- as.data.frame(rowSums(asv_mat_pool_rare %>% select(all_of(ZoneOxique))))
colnames(Oxique_asv) <- "TotalOxique"
for (i in row.names(Oxique_asv)) {if (Oxique_asv[i,"TotalOxique"] > 0) {Oxique_asv[i,"TotalOxique"] <- 1}}
Oxique_asv$ASV_Id <- row.names(Oxique_asv)
##Anoxique
Anoxique_asv <- as.data.frame(rowSums(asv_mat_pool_rare %>% select(all_of(ZoneAnoxique))))
colnames(Anoxique_asv) <- "TotalAnoxique"
for (i in row.names(Anoxique_asv)) {if (Anoxique_asv[i,"TotalAnoxique"] > 0) {Anoxique_asv[i,"TotalAnoxique"] <- 1}}
Anoxique_asv$ASV_Id <- row.names(Anoxique_asv)
##Zone
Zone_asv <- merge(x = Oxique_asv,y = Anoxique_asv, by = "ASV_Id")

#Fraction
##Petite
Petite_asv <- as.data.frame(rowSums(asv_mat_pool_rare %>% select(all_of(FractionPetite))))
colnames(Petite_asv) <- "TotalPetite"
for (i in row.names(Petite_asv)) {if (Petite_asv[i,"TotalPetite"] > 0) {Petite_asv[i,"TotalPetite"] <- 1}}
Petite_asv$ASV_Id <- row.names(Petite_asv)
##Grande
Grande_asv <- as.data.frame(rowSums(asv_mat_pool_rare %>% select(all_of(FractionGrande))))
colnames(Grande_asv) <- "TotalGrande"
for (i in row.names(Grande_asv)) {if (Grande_asv[i,"TotalGrande"] > 0) {Grande_asv[i,"TotalGrande"] <- 1}}
Grande_asv$ASV_Id <- row.names(Grande_asv)
##Fraction
Fraction_asv <- merge(x = Petite_asv,y = Grande_asv, by = "ASV_Id")

#Date
##04
`04_asv` <- as.data.frame(rowSums(asv_mat_pool_rare %>% select(all_of(Date04))))
colnames(`04_asv`) <- "Total04"
for (i in row.names(`04_asv`)) {if (`04_asv`[i,"Total04"] > 0) {`04_asv`[i,"Total04"] <- 1}}
`04_asv`$ASV_Id <- row.names(`04_asv`)
##06
`06_asv` <- as.data.frame(rowSums(asv_mat_pool_rare %>% select(all_of(Date06))))
colnames(`06_asv`) <- "Total06"
for (i in row.names(`06_asv`)) {if (`06_asv`[i,"Total06"] > 0) {`06_asv`[i,"Total06"] <- 1}}
`06_asv`$ASV_Id <- row.names(`06_asv`)
##09
`09_asv` <- as.data.frame(rowSums(asv_mat_pool_rare %>% select(all_of(Date09))))
colnames(`09_asv`) <- "Total09"
for (i in row.names(`09_asv`)) {if (`09_asv`[i,"Total09"] > 0) {`09_asv`[i,"Total09"] <- 1}}
`09_asv`$ASV_Id <- row.names(`09_asv`)
##11
`11_asv` <- as.data.frame(rowSums(asv_mat_pool_rare %>% select(all_of(Date11))))
colnames(`11_asv`) <- "Total11"
for (i in row.names(`11_asv`)) {if (`11_asv`[i,"Total11"] > 0) {`11_asv`[i,"Total11"] <- 1}}
`11_asv`$ASV_Id <- row.names(`11_asv`)
##Date
Date_asv <- merge(x = `04_asv`,y = `06_asv`, by = "ASV_Id")
Date_asv <- merge(x = Date_asv,y = `09_asv`, by = "ASV_Id")
Date_asv <- merge(x = Date_asv,y = `11_asv`, by = "ASV_Id")

# AFC Plot BASE ----------------------------------------------------------------
  # Séquence ----------------------------------------------------------------
dt <- as.data.frame(seq_mat_pool_rare)
dt <- dt %>% select(-"Taxonomy")
#CA
res.ca <- CA(dt, graph = FALSE,ncp = 2 )
fviz_ca_row(res.ca, repel = FALSE, label = "none")
p <- get_ca_row(res.ca)
coord <- p$coord
coord <- as.data.frame(coord)
coord$ASV_Id <- row.names(coord)
##Cycle
data_seq <- merge(x = coord, y = Cycle_seq, by = "ASV_Id")
### 100% - 0%
  data_seq$Jour <- 0
  data_seq$Nuit <- 0
  for (i in rownames(data_seq)) { if ( data_seq[i,"TotalJour"] != 0) { data_seq[i,"Jour"] <- 1}}
  for (i in rownames(data_seq)) { if ( data_seq[i,"TotalNuit"] != 0) { data_seq[i,"Nuit"] <- 2}}
  data_seq$Cycle <- 0
  for (i in rownames(data_seq)) { data_seq[i,"Cycle"] <- data_seq[i,"Jour"] + data_seq[i,"Nuit"] 
  if ( data_seq[i,"Cycle"] == 1) { data_seq[i,"Cycle"] <- "Jour"}
  if ( data_seq[i,"Cycle"] == 2) { data_seq[i,"Cycle"] <- "Nuit"}
  if ( data_seq[i,"Cycle"] == 3) { data_seq[i,"Cycle"] <- "Communs"}}
### 90% - 10%
  data_seq$Cycle90 <- "Communs"
  for (i in rownames(data_seq)) { if ( 0.1*data_seq[i,"TotalJour"] > data_seq[i,"TotalNuit"]) { data_seq[i,"Cycle90"] <- "Jour"}}
  for (i in rownames(data_seq)) { if ( 0.1*data_seq[i,"TotalNuit"] > data_seq[i,"TotalJour"]) { data_seq[i,"Cycle90"] <- "Nuit"}}

##Zone
data_seq <- merge(x = data_seq, y = Zone_seq, by = "ASV_Id")
### 100% - 0%
  data_seq$Oxique <- 0
  data_seq$Anoxique <- 0
  for (i in rownames(data_seq)) { if ( data_seq[i,"TotalOxique"] != 0) { data_seq[i,"Oxique"] <- 1}}
  for (i in rownames(data_seq)) { if ( data_seq[i,"TotalAnoxique"] != 0) { data_seq[i,"Anoxique"] <- 2}}
  data_seq$Zone <- 0
  for (i in rownames(data_seq)) { data_seq[i,"Zone"] <- data_seq[i,"Oxique"] + data_seq[i,"Anoxique"] 
  if ( data_seq[i,"Zone"] == 1) { data_seq[i,"Zone"] <- "Oxique"}
  if ( data_seq[i,"Zone"] == 2) { data_seq[i,"Zone"] <- "Anoxique"}
  if ( data_seq[i,"Zone"] == 3) { data_seq[i,"Zone"] <- "Communs"}}
### 90% - 10%
  data_seq$Zone90 <- "Communs"
  for (i in rownames(data_seq)) { if ( 0.1*data_seq[i,"TotalOxique"] > data_seq[i,"TotalAnoxique"]) { data_seq[i,"Zone90"] <- "Oxique"}}
  for (i in rownames(data_seq)) { if ( 0.1*data_seq[i,"TotalAnoxique"] > data_seq[i,"TotalOxique"]) { data_seq[i,"Zone90"] <- "Anoxique"}}  

##Fraction
data_seq <- merge(x = data_seq, y = Fraction_seq, by = "ASV_Id")
### 100% - 0%
  data_seq$Petite <- 0
  data_seq$Grande <- 0
  for (i in rownames(data_seq)) { if ( data_seq[i,"TotalPetite"] != 0) { data_seq[i,"Petite"] <- 1}}
  for (i in rownames(data_seq)) { if ( data_seq[i,"TotalGrande"] != 0) { data_seq[i,"Grande"] <- 2}}
  data_seq$Fraction <- 0
  for (i in rownames(data_seq)) { data_seq[i,"Fraction"] <- data_seq[i,"Petite"] + data_seq[i,"Grande"] 
  if ( data_seq[i,"Fraction"] == 1) { data_seq[i,"Fraction"] <- "Petite"}
  if ( data_seq[i,"Fraction"] == 2) { data_seq[i,"Fraction"] <- "Grande"}
  if ( data_seq[i,"Fraction"] == 3) { data_seq[i,"Fraction"] <- "Communs"}}
### 90% - 10%
  data_seq$Fraction90 <- "Communs"
  for (i in rownames(data_seq)) { if ( 0.1*data_seq[i,"TotalPetite"] > data_seq[i,"TotalGrande"]) { data_seq[i,"Fraction90"] <- "Petite"}}
  for (i in rownames(data_seq)) { if ( 0.1*data_seq[i,"TotalGrande"] > data_seq[i,"TotalPetite"]) { data_seq[i,"Fraction90"] <- "Grande"}}  
  
#Date
data_seq <- merge(x = data_seq, y = Date_seq, by = "ASV_Id")
data_seq04 <- data_seq %>% filter(Total04 != 0) %>% select(-Total06,-Total09,-Total11)
data_seq06 <- data_seq %>% filter(Total06 != 0) %>% select(-Total04,-Total09,-Total11)
data_seq09 <- data_seq %>% filter(Total09 != 0) %>% select(-Total04,-Total06,-Total11)
data_seq11 <- data_seq %>% filter(Total11 != 0) %>% select(-Total04,-Total06,-Total09)

#data_seq$Avril <- 0
#data_seq$Juin <- 0
#data_seq$Septembre <- 0
#data_seq$Novembre <- 0
#for (i in rownames(data_seq)) { if ( data_seq[i,"Total04"] != 0) { data_seq[i,"Avril"] <- 1}}
#for (i in rownames(data_seq)) { if ( data_seq[i,"Total06"] != 0) { data_seq[i,"Juin"] <- 3}}
#for (i in rownames(data_seq)) { if ( data_seq[i,"Total09"] != 0) { data_seq[i,"Septembre"] <- 5}}
#for (i in rownames(data_seq)) { if ( data_seq[i,"Total11"] != 0) { data_seq[i,"Novembre"] <- 7}}
#data_seq$Date <- 0
#for (i in rownames(data_seq)) { data_seq[i,"Date"] <- data_seq[i,"Avril"] + data_seq[i,"Juin"] + data_seq[i,"Septembre"] + data_seq[i,"Novembre"]
#if ( data_seq[i,"Date"] == 1) { data_seq[i,"Date"] <- "Avril"}
#if ( data_seq[i,"Date"] == 3) { data_seq[i,"Date"] <- "Juin"}
#if ( data_seq[i,"Date"] == 5) { data_seq[i,"Date"] <- "Septembre"}
#if ( data_seq[i,"Date"] == 7) { data_seq[i,"Date"] <- "Novembre"}
#else { data_seq[i,"Date"] <- "Communs" }}

Xseq <- fviz_screeplot (res.ca, addlabels = TRUE, ylim = c(0, 50))
Xseq <- Xseq$data
Dim1seq <- paste("Dim 1 [",round(Xseq %>% filter(dim == 1) %>% select(eig),1),"%]", sep = "")
Dim2seq <- paste("Dim 2 [",round(Xseq %>% filter(dim == 2) %>% select(eig),1),"%]", sep = "")

  # ASV ----------------------------------------------------------------
dt <- as.data.frame(asv_mat_pool_rare)
dt <- dt %>% select(-"Taxonomy")
#CA
res.ca <- CA(dt, graph = FALSE,ncp = 2 )
fviz_ca_row(res.ca, repel = FALSE, label = "none")
p <- get_ca_row(res.ca)
coord <- p$coord
coord <- as.data.frame(coord)
coord$ASV_Id <- row.names(coord)
##Cycle
data_asv <- merge(x = coord, y = Cycle_asv, by = "ASV_Id")
data_asv$Jour <- 0
data_asv$Nuit <- 0
for (i in rownames(data_asv)) { if ( data_asv[i,"TotalJour"] != 0) { data_asv[i,"Jour"] <- 1}}
for (i in rownames(data_asv)) { if ( data_asv[i,"TotalNuit"] != 0) { data_asv[i,"Nuit"] <- 2}}
data_asv$Cycle <- 0
for (i in rownames(data_asv)) { data_asv[i,"Cycle"] <- data_asv[i,"Jour"] + data_asv[i,"Nuit"] 
if ( data_asv[i,"Cycle"] == 1) { data_asv[i,"Cycle"] <- "Jour"}
if ( data_asv[i,"Cycle"] == 2) { data_asv[i,"Cycle"] <- "Nuit"}
if ( data_asv[i,"Cycle"] == 3) { data_asv[i,"Cycle"] <- "Communs"}}
##Zone
data_asv <- merge(x = data_asv, y = Zone_asv, by = "ASV_Id")
data_asv$Oxique <- 0
data_asv$Anoxique <- 0
for (i in rownames(data_asv)) { if ( data_asv[i,"TotalOxique"] != 0) { data_asv[i,"Oxique"] <- 1}}
for (i in rownames(data_asv)) { if ( data_asv[i,"TotalAnoxique"] != 0) { data_asv[i,"Anoxique"] <- 2}}
data_asv$Zone <- 0
for (i in rownames(data_asv)) { data_asv[i,"Zone"] <- data_asv[i,"Oxique"] + data_asv[i,"Anoxique"] 
if ( data_asv[i,"Zone"] == 1) { data_asv[i,"Zone"] <- "Oxique"}
if ( data_asv[i,"Zone"] == 2) { data_asv[i,"Zone"] <- "Anoxique"}
if ( data_asv[i,"Zone"] == 3) { data_asv[i,"Zone"] <- "Communs"}}
##Fraction
data_asv <- merge(x = data_asv, y = Fraction_asv, by = "ASV_Id")
data_asv$Petite <- 0
data_asv$Grande <- 0
for (i in rownames(data_asv)) { if ( data_asv[i,"TotalPetite"] != 0) { data_asv[i,"Petite"] <- 1}}
for (i in rownames(data_asv)) { if ( data_asv[i,"TotalGrande"] != 0) { data_asv[i,"Grande"] <- 2}}
data_asv$Fraction <- 0
for (i in rownames(data_asv)) { data_asv[i,"Fraction"] <- data_asv[i,"Petite"] + data_asv[i,"Grande"] 
if ( data_asv[i,"Fraction"] == 1) { data_asv[i,"Fraction"] <- "Petite"}
if ( data_asv[i,"Fraction"] == 2) { data_asv[i,"Fraction"] <- "Grande"}
if ( data_asv[i,"Fraction"] == 3) { data_asv[i,"Fraction"] <- "Communs"}}
##Date
data_asv <- merge(x = data_asv, y = Date_asv, by = "ASV_Id")
data_asv04 <- data_asv %>% filter(Total04 != 0) %>% select(-Total06,-Total09,-Total11)
data_asv06 <- data_asv %>% filter(Total06 != 0) %>% select(-Total04,-Total09,-Total11)
data_asv09 <- data_asv %>% filter(Total09 != 0) %>% select(-Total04,-Total06,-Total11)
data_asv11 <- data_asv %>% filter(Total11 != 0) %>% select(-Total04,-Total06,-Total09)

Xasv <- fviz_screeplot (res.ca, addlabels = TRUE, ylim = c(0, 50))
Xasv <- Xasv$data
Dim1asv <- paste("Dim 1 [",round(Xasv %>% filter(dim == 1) %>% select(eig),1),"%]", sep = "")
Dim2asv <- paste("Dim 2 [",round(Xasv %>% filter(dim == 2) %>% select(eig),1),"%]", sep = "")




# AFC Taxonomy ------------------------------------------------------------
  # Séquence ----------------------------------------------------------------
data_seq_tax <- data_seq
colnames(tax_tablemix) <- c("ASV_Id","Taxonomy")
data_seq_tax <- merge(data_seq_tax,tax_tablemix, by = "ASV_Id")
data_seq_tax <- separate(data_seq_tax, Taxonomy,c("Domain","Superphylum","Phylum","Class","Order","Family","Genus","Species"), sep =";")
data_seq_tax[data_seq_tax == ""] <- NA
data_seq_tax[data_seq_tax == " "] <- NA
data_seq_tax[data_seq_tax == "NA"] <- NA
data_seq_tax[is.na(data_seq_tax) == TRUE] <- "Not Affiliated"

    # Export data_seq_tax for Krona -------------------------------------------
## Jour
write.table(data_seq_tax %>% select(TotalJour, Domain, Superphylum, Phylum, Class, Order, Family, Genus, Species), file = "Krona/data_seq_Jour.csv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
## Nuit
write.table(data_seq_tax %>% select(TotalNuit, Domain, Superphylum, Phylum, Class, Order, Family, Genus, Species), file = "Krona/data_seq_Nuit.csv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
## Oxique
write.table(data_seq_tax %>% select(TotalOxique, Domain, Superphylum, Phylum, Class, Order, Family, Genus, Species), file = "Krona/data_seq_Oxique.csv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
## Anoxique
write.table(data_seq_tax %>% select(TotalAnoxique, Domain, Superphylum, Phylum, Class, Order, Family, Genus, Species), file = "Krona/data_seq_Anoxique.csv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
## Grande
write.table(data_seq_tax %>% select(TotalGrande, Domain, Superphylum, Phylum, Class, Order, Family, Genus, Species), file = "Krona/data_seq_Grande.csv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
## Petite
write.table(data_seq_tax %>% select(TotalPetite, Domain, Superphylum, Phylum, Class, Order, Family, Genus, Species), file = "Krona/data_seq_Petite.csv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
## Avril
write.table(data_seq_tax %>% select(Total04, Domain, Superphylum, Phylum, Class, Order, Family, Genus, Species), file = "Krona/data_seq_Avril.csv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
## Juin
write.table(data_seq_tax %>% select(Total06, Domain, Superphylum, Phylum, Class, Order, Family, Genus, Species), file = "Krona/data_seq_Juin.csv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
## Septembre
write.table(data_seq_tax %>% select(Total09, Domain, Superphylum, Phylum, Class, Order, Family, Genus, Species), file = "Krona/data_seq_Septembre.csv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
## Novembre
write.table(data_seq_tax %>% select(Total11, Domain, Superphylum, Phylum, Class, Order, Family, Genus, Species), file = "Krona/data_seq_Novembre.csv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
## Total
data_seq_tax_Sum <- data_seq_tax %>% select(TotalJour, TotalNuit, Domain, Superphylum, Phylum, Class, Order, Family, Genus, Species)
data_seq_tax_Sum$Total <- data_seq_tax_Sum$TotalJour + data_seq_tax_Sum$TotalNuit
data_seq_tax_Sum <- data_seq_tax_Sum %>% select(Total, Domain, Superphylum, Phylum, Class, Order, Family, Genus, Species)
write.table(data_seq_tax_Sum %>% select(Total, Domain, Superphylum, Phylum, Class, Order, Family, Genus, Species), file = "Krona/data_seq_Total.csv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

## Launch Krona
system("for t in $(ls Krona/)
do
perl ../../../script/KronaTools-2.8/scripts/ImportText.pl Krona/$t -o Krona/$t.html
done")

# Sort interest division
if (Mode == "Class") {
  data_seq_tax <- data_seq_tax %>% filter(Phylum == all_of(Group))
  data_seq_tax$Division <- data_seq_tax$Class }
if (Mode == "Phylum") {
data_seq_tax <- data_seq_tax %>% filter(Superphylum == all_of(Group))
data_seq_tax$Division <- data_seq_tax$Phylum }
if (Mode == "Superphylum") {
data_seq_tax$Division <- data_seq_tax$Superphylum }


    # 100% - 0% -------------------------------------------------------------
      # All_2018 ----------------------------------------------------------
##Figure Cycles
aCycle <- ggplot(data_seq_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Cycle)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Cycle, alpha = Cycle),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(aCycle)

##Figure Zone
bZone <- ggplot(data_seq_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Zone)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Zone, alpha = Zone),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Zones", alpha = "Zones", linetype = "Zones")
print(bZone)
##Figure Fraction
cFraction <- ggplot(data_seq_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Fraction)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Fraction, alpha = Fraction),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Fractions", alpha = "Fractions", linetype = "Fractions")
print(cFraction)
##Coplot
svglite("AFC-Distribution/AFC-Sequence-100.svg",width = 12.00,height = 8.00)
All_plot <- plot_grid(aCycle, cFraction, bZone, ncol = 2, nrow = 2, rel_widths = c(3,3), rel_heights = c(3,3))
print(All_plot)
dev.off()
      # Avril_2018 ---------------------------------------------------
data_seq04_tax <- merge(x = data_seq04, y =  tax_tablemix, by = "ASV_Id")

a_04 <- ggplot(data_seq04_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Cycle)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Cycle, alpha = Cycle),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(a_04)

#Figure Zone
b_04 <- ggplot(data_seq04_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Zone)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Zone, alpha = Zone),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Zones", alpha = "Zones", linetype = "Zones")
print(b_04)

#Figure Fraction
c_04 <- ggplot(data_seq04_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Fraction)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Fraction, alpha = Fraction),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Fraction", alpha = "Fraction", linetype = "Fraction")
print(c_04)

svglite("AFC-Distribution/AFC-Sequence-Avril-100.svg",width = 12.00,height = 8.00)
Avril_plot <- plot_grid(a_04, c_04, b_04, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(Avril_plot)
dev.off()  
      # Juin_2018 ---------------------------------------------------
data_seq06_tax <- merge(x = data_seq06, y =  tax_tablemix, by = "ASV_Id")

a_06 <- ggplot(data_seq06_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Cycle)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Cycle, alpha = Cycle),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(a_06)

#Figure Zone
b_06 <- ggplot(data_seq06_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Zone)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Zone, alpha = Zone),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Zones", alpha = "Zones", linetype = "Zones")
print(b_06)

#Figure Fraction
c_06 <- ggplot(data_seq06_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Fraction)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Fraction, alpha = Fraction),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Fraction", alpha = "Fraction", linetype = "Fraction")
print(c_06)

svglite("AFC-Distribution/AFC-Sequence-Juin-100.svg",width = 12.00,height = 8.00)
Juin_plot <- plot_grid(a_06, c_06, b_06, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(Juin_plot)
dev.off()  

      # Septembre_2018 ---------------------------------------------------
data_seq09_tax <- merge(x = data_seq09, y =  tax_tablemix, by = "ASV_Id")

a_09 <- ggplot(data_seq09_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Cycle)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Cycle, alpha = Cycle),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(a_09)

#Figure Zone
b_09 <- ggplot(data_seq09_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Zone)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Zone, alpha = Zone),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Zones", alpha = "Zones", linetype = "Zones")
print(b_09)

#Figure Fraction
c_09 <- ggplot(data_seq09_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Fraction)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Fraction, alpha = Fraction),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Fraction", alpha = "Fraction", linetype = "Fraction")
print(c_09)

svglite("AFC-Distribution/AFC-Sequence-Septembre-100.svg",width = 12.00,height = 8.00)
Septembre_plot <- plot_grid(a_09, c_09, b_09, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(Septembre_plot)
dev.off()  


      # Novembre_2018 ---------------------------------------------------
data_seq11_tax <- merge(x = data_seq11, y =  tax_tablemix, by = "ASV_Id")

a_11 <- ggplot(data_seq11_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Cycle)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Cycle, alpha = Cycle),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(a_11)

#Figure Zone
b_11 <- ggplot(data_seq11_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Zone)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Zone, alpha = Zone),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Zones", alpha = "Zones", linetype = "Zones")
print(b_11)

#Figure Fraction
c_11 <- ggplot(data_seq11_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Fraction)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Fraction, alpha = Fraction),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Fraction", alpha = "Fraction", linetype = "Fraction")
print(c_11)

svglite("AFC-Distribution/AFC-Sequence-Novembre-100.svg",width = 12.00,height = 8.00)
Novembre_plot <- plot_grid(a_11, c_11, b_11, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(Novembre_plot)
dev.off()  


    # 90% - 10% -------------------------------------------------------------
      # All_2018 ----------------------------------------------------------
##Figure Cycles
aCycle90 <- ggplot(data_seq_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Cycle90)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Cycle90, alpha = Cycle90),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(aCycle90)

##Figure Zone
bZone90 <- ggplot(data_seq_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Zone90)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Zone90, alpha = Zone90),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Zones", alpha = "Zones", linetype = "Zones")
print(bZone90)
##Figure Fraction
cFraction90 <- ggplot(data_seq_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Fraction90)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Fraction90, alpha = Fraction90),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Fractions", alpha = "Fractions", linetype = "Fractions")
print(cFraction90)
##Coplot
svglite("AFC-Distribution/AFC-Sequence-90.svg",width = 12.00,height = 8.00)
All_plot90 <- plot_grid(aCycle90, cFraction90, bZone90, ncol = 2, nrow = 2, rel_widths = c(3,3), rel_heights = c(3,3))
print(All_plot90)
dev.off()
      # Avril_2018 ---------------------------------------------------
data_seq04_tax <- merge(x = data_seq04, y =  tax_tablemix, by = "ASV_Id")

a_04_90 <- ggplot(data_seq04_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Cycle90)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Cycle90, alpha = Cycle90),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(a_04_90)

#Figure Zone
b_04_90 <- ggplot(data_seq04_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Zone90)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Zone90, alpha = Zone90),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Zones", alpha = "Zones", linetype = "Zones")
print(b_04_90)

#Figure Fraction
c_04_90 <- ggplot(data_seq04_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Fraction90)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Fraction90, alpha = Fraction90),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Fraction", alpha = "Fraction", linetype = "Fraction")
print(c_04_90)

svglite("AFC-Distribution/AFC-Sequence-Avril-90.svg",width = 12.00,height = 8.00)
Avril_plot90 <- plot_grid(a_04_90, c_04_90, b_04_90, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(Avril_plot90)
dev.off()  
      # Juin_2018 ---------------------------------------------------
data_seq06_tax <- merge(x = data_seq06, y =  tax_tablemix, by = "ASV_Id")

a_06_90 <- ggplot(data_seq06_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Cycle90)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Cycle90, alpha = Cycle90),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(a_06_90)

#Figure Zone
b_06_90 <- ggplot(data_seq06_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Zone90)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Zone90, alpha = Zone90),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Zones", alpha = "Zones", linetype = "Zones")
print(b_06_90)

#Figure Fraction
c_06_90 <- ggplot(data_seq06_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Fraction90)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Fraction90, alpha = Fraction90),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Fraction", alpha = "Fraction", linetype = "Fraction")
print(c_06_90)

svglite("AFC-Distribution/AFC-Sequence-Juin-90.svg",width = 12.00,height = 8.00)
Juin_plot_90 <- plot_grid(a_06_90, c_06_90, b_06_90, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(Juin_plot_90)
dev.off()  

      # Septembre_2018 ---------------------------------------------------
data_seq09_tax <- merge(x = data_seq09, y =  tax_tablemix, by = "ASV_Id")

a_09_90 <- ggplot(data_seq09_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Cycle90)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Cycle90, alpha = Cycle90),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(a_09_90)

#Figure Zone
b_09_90 <- ggplot(data_seq09_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Zone90)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Zone90, alpha = Zone90),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Zones", alpha = "Zones", linetype = "Zones")
print(b_09_90)

#Figure Fraction
c_09_90 <- ggplot(data_seq09_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Fraction90)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Fraction90, alpha = Fraction90),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Fraction", alpha = "Fraction", linetype = "Fraction")
print(c_09_90)

svglite("AFC-Distribution/AFC-Sequence-Septembre-90.svg",width = 12.00,height = 8.00)
Septembre_plot_90 <- plot_grid(a_09_90, c_09_90, b_09_90, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(Septembre_plot_90)
dev.off()  


      # Novembre_2018 ---------------------------------------------------
data_seq11_tax <- merge(x = data_seq11, y =  tax_tablemix, by = "ASV_Id")

a_11_90 <- ggplot(data_seq11_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Cycle90)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Cycle90, alpha = Cycle90),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(a_11_90)

#Figure Zone
b_11_90 <- ggplot(data_seq11_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Zone90)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Zone90, alpha = Zone90),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Zones", alpha = "Zones", linetype = "Zones")
print(b_11_90)

#Figure Fraction
c_11_90 <- ggplot(data_seq11_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Fraction90)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Fraction90, alpha = Fraction90),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1seq,y=Dim2seq,color = "Fraction", alpha = "Fraction", linetype = "Fraction")
print(c_11_90)

svglite("AFC-Distribution/AFC-Sequence-Novembre-90.svg",width = 12.00,height = 8.00)
Novembre_plot_90 <- plot_grid(a_11_90, c_11_90, b_11_90, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(Novembre_plot_90)
dev.off()  


  # ASV ----------------------------------------------------------------
data_asv_tax <- data_asv
colnames(tax_tablemix) <- c("ASV_Id","Taxonomy")
data_asv_tax <- merge(x = data_asv_tax, y =  tax_tablemix, by = "ASV_Id")
data_asv_tax <- separate(data_asv_tax, Taxonomy,c("Domain","Superphylum","Phylum","Class","Order","Family","Genus","Species"), sep =";")
data_asv_tax[data_asv_tax == ""] <- NA
data_asv_tax[data_asv_tax == " "] <- NA
data_asv_tax[data_asv_tax == "NA"] <- NA
data_asv_tax[is.na(data_asv_tax) == TRUE] <- "Not Affiliated"

if (Mode == "Class") {
  data_asv_tax <- data_asv_tax %>% filter(Phylum == all_of(Group))
  data_asv_tax$Division <- data_asv_tax$Class }
if (Mode == "Phylum") {
data_asv_tax <- data_asv_tax %>% filter(Superphylum == all_of(Group))
data_asv_tax$Division <- data_asv_tax$Phylum }
if (Mode == "Superphylum") {
data_asv_tax$Division <- data_asv_tax$Superphylum }
    # All_2018 ----------------------------------------------------------
##Figure Cycles
iCycle <- ggplot(data_asv_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Cycle)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Cycle, alpha = Cycle),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1asv,y=Dim2asv,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(iCycle)

##Figure Zone
jZone <- ggplot(data_asv_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Zone)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Zone, alpha = Zone),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1asv,y=Dim2asv,color = "Zones", alpha = "Zones", linetype = "Zones")
print(jZone)

##Figure Fraction
kFraction <- ggplot(data_asv_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Fraction)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Fraction, alpha = Fraction),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1asv,y=Dim2asv,color = "Fractions", alpha = "Fractions", linetype = "Fractions")
print(kFraction)
##Coplot
svglite("AFC-Distribution/AFC-ASV.svg",width = 12.00,height = 8.00)
All_plot <- plot_grid(iCycle, jZone, kFraction, ncol = 2, nrow = 2, rel_widths = c(3,3), rel_heights = c(3,3))
print(All_plot)
dev.off()
    # Avril_2018 ---------------------------------------------------
data_asv04_tax <- merge(x = data_asv04, y =  tax_tablemix, by = "ASV_Id")

i_04 <- ggplot(data_asv04_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Cycle)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Cycle, alpha = Cycle),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1asv,y=Dim2asv,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(i_04)

#Figure Zone
j_04 <- ggplot(data_asv04_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Zone)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Zone, alpha = Zone),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1asv,y=Dim2asv,color = "Zones", alpha = "Zones", linetype = "Zones")
print(j_04)

#Figure Fraction
k_04 <- ggplot(data_asv04_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Fraction)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Fraction, alpha = Fraction),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1asv,y=Dim2asv,color = "Fraction", alpha = "Fraction", linetype = "Fraction")
print(k_04)

svglite("AFC-Distribution/AFC-ASV-Avril.svg",width = 12.00,height = 8.00)
Avril_plot <- plot_grid(i_04, j_04, k_04, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(Avril_plot)
dev.off()  
    # Juin_2018 ---------------------------------------------------
data_asv06_tax <- merge(x = data_asv06, y =  tax_tablemix, by = "ASV_Id")

i_06 <- ggplot(data_asv06_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Cycle)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Cycle, alpha = Cycle),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1asv,y=Dim2asv,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(i_06)

#Figure Zone
j_06 <- ggplot(data_asv06_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Zone)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Zone, alpha = Zone),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1asv,y=Dim2asv,color = "Zones", alpha = "Zones", linetype = "Zones")
print(j_06)

#Figure Fraction
k_06 <- ggplot(data_asv06_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Fraction)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Fraction, alpha = Fraction),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1asv,y=Dim2asv,color = "Fraction", alpha = "Fraction", linetype = "Fraction")
print(k_06)

svglite("AFC-Distribution/AFC-ASV-Juin.svg",width = 12.00,height = 8.00)
Juin_plot <- plot_grid(i_06, j_06, k_06, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(Juin_plot)
dev.off()  

    # Septembre_2018 ---------------------------------------------------
data_asv09_tax <- merge(x = data_asv09, y =  tax_tablemix, by = "ASV_Id")

i_09 <- ggplot(data_asv09_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Cycle)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Cycle, alpha = Cycle),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1asv,y=Dim2asv,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(i_09)

#Figure Zone
j_09 <- ggplot(data_asv09_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Zone)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Zone, alpha = Zone),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1asv,y=Dim2asv,color = "Zones", alpha = "Zones", linetype = "Zones")
print(j_09)

#Figure Fraction
k_09 <- ggplot(data_asv09_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Fraction)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Fraction, alpha = Fraction),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1asv,y=Dim2asv,color = "Fraction", alpha = "Fraction", linetype = "Fraction")
print(k_09)

svglite("AFC-Distribution/AFC-ASV-Septembre.svg",width = 12.00,height = 8.00)
Septembre_plot <- plot_grid(i_09, j_09, k_09, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(Septembre_plot)
dev.off()  


    # Novembre_2018 ---------------------------------------------------
data_asv11_tax <- merge(x = data_asv11, y =  tax_tablemix, by = "ASV_Id")

i_11 <- ggplot(data_asv11_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Cycle)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Cycle, alpha = Cycle),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1asv,y=Dim2asv,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(i_11)

#Figure Zone
j_11 <- ggplot(data_asv11_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Zone)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Zone, alpha = Zone),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1asv,y=Dim2asv,color = "Zones", alpha = "Zones", linetype = "Zones")
print(j_11)

#Figure Fraction
k_11 <- ggplot(data_asv11_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Fraction)) + geom_point(size = 2) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  stat_ellipse(aes(linetype = Fraction, alpha = Fraction),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1asv,y=Dim2asv,color = "Fraction", alpha = "Fraction", linetype = "Fraction")
print(k_11)

svglite("AFC-Distribution/AFC-ASV-Novembre.svg",width = 12.00,height = 8.00)
Novembre_plot <- plot_grid(i_11, j_11, k_11, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(Novembre_plot)
dev.off()  


# Venn diagramm -----------------------------------------------------------
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  # Cycle
  CycleVenn <- venn.diagram(x = list(data_asv_tax %>% filter(Jour == 1) %>% select(ASV_Id) %>% unlist(), data_asv_tax %>% filter(Nuit == 2) %>% select(ASV_Id) %>% unlist()), 
                            category.names = c("Jour" , "Nuit"),
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
  # Zone
  ZoneVenn <- venn.diagram(x = list(data_asv_tax %>% filter(Oxique == 1) %>% select(ASV_Id) %>% unlist(), data_asv_tax %>% filter(Anoxique == 2) %>% select(ASV_Id) %>% unlist()),
                          category.names = c("Oxique" , "Anoxique"),
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
  # Fraction
  FractionVenn <- venn.diagram(x = list(data_asv_tax %>% filter(Petite == 1) %>% select(ASV_Id) %>% unlist(), data_asv_tax %>% filter(Grande == 2) %>% select(ASV_Id) %>% unlist()),
                           category.names = c("Petite" , "Grande"),
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
  # Date
  DateVenn <- venn.diagram(x = list(data_asv_tax %>% filter(Total04 == 1) %>% select(ASV_Id) %>% unlist(), data_asv_tax %>% filter(Total06 == 1) %>% select(ASV_Id) %>% unlist(),data_asv_tax %>% filter(Total09 == 1) %>% select(ASV_Id) %>% unlist(),data_asv_tax %>% filter(Total11 == 1) %>% select(ASV_Id) %>% unlist()),
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
  ggsave(DateVenn, file="Venn/Date.svg", device = "svg", width = 3, height = 3)
  

# HIST  --------------------------------------------------------------------
  # ASV TOTAL---------------------------------------------------------------------
    # Cycle -------------------------------------------------------------------
## Jour
totalJourASV <- data_asv_tax %>% select(ASV_Id,TotalJour,Division)
row.names(totalJourASV)<-totalJourASV$ASV_Id ; totalJourASV <- totalJourASV %>% select(-ASV_Id)
totalJourASV <- totalJourASV %>% group_by(Division) %>% summarise_all(sum)
totalJourASV$TotalJour <- totalJourASV$TotalJour*100/sum(totalJourASV$TotalJour)
totalJourASV <- totalJourASV %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalJour) - 0.5*TotalJour)
totalJourASV$label <- paste(round(totalJourASV$TotalJour,1), "%", sep = "")
for (i in rownames(totalJourASV)) {
  if (totalJourASV[i,"label"] == "0%") { totalJourASV[i,"label"] <- NA}}
for (i in rownames(totalJourASV)) {
  if (is.na(totalJourASV[i,"label"]) == FALSE) { totalJourASV[i,"label"] <- paste(totalJourASV[i,"Division"]," : ",totalJourASV[i,"label"], sep = "")}}
totalJourASV$Cycle<- rep("Jour", each = nrow(totalJourASV))
##Nuit
totalNuitASV <- data_asv_tax %>% select(ASV_Id,TotalNuit,Division)
row.names(totalNuitASV)<-totalNuitASV$ASV_Id ; totalNuitASV <- totalNuitASV %>% select(-ASV_Id)
totalNuitASV <- totalNuitASV %>% group_by(Division) %>% summarise_all(sum)
totalNuitASV$TotalNuit <- totalNuitASV$TotalNuit*100/sum(totalNuitASV$TotalNuit)
totalNuitASV <- totalNuitASV %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalNuit) - 0.5*TotalNuit)
totalNuitASV$label <- paste(round(totalNuitASV$TotalNuit,1), "%", sep = "")
for (i in rownames(totalNuitASV)) {
  if (totalNuitASV[i,"label"] == "0%") { totalNuitASV[i,"label"] <- NA}}
for (i in rownames(totalNuitASV)) {
  if (is.na(totalNuitASV[i,"label"]) == FALSE) { totalNuitASV[i,"label"] <- paste(totalNuitASV[i,"Division"]," : ",totalNuitASV[i,"label"], sep = "")}}
totalNuitASV$Cycle<- rep("Nuit", each = nrow(totalNuitASV))
##Cycle
colnames(totalJourASV)[2]  <- "value"
totalJourASV$Sum <- rep(0, each = nrow(totalJourASV))
totalJourASV$Count <- rep(0, each = nrow(totalJourASV))
for (i in totalJourASV$Division) { totalJourASV$Count[which(totalJourASV$Division == i)] <- sum(data_asv_tax  %>% filter(Division == i) %>% select(TotalJour))}
totalJourASV$Sum <- sum(totalJourASV$Count)
colnames(totalNuitASV)[2]  <- "value"
totalNuitASV$Sum <- rep(0, each = nrow(totalNuitASV))
totalNuitASV$Count <- rep(0, each = nrow(totalNuitASV))
for (i in totalNuitASV$Division) { totalNuitASV$Count[which(totalNuitASV$Division == i)] <- sum(data_asv_tax  %>% filter(Division == i) %>% select(TotalNuit))}
totalNuitASV$Sum <- sum(totalNuitASV$Count)
totalCycleASV <- rbind(totalJourASV,totalNuitASV)

#Figure
ay <- ggplot(totalCycleASV, mapping = aes(y= value, x = Cycle, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + scale_fill_manual(values = palette) + 
  geom_label(aes(y = 106,label = paste(Sum," ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
legendASV <- get_legend(ay)
ay <- ay + labs(x="Cycles",y="ASVs (%)") + theme(legend.position = "none")
print(ay)

    # Zone -------------------------------------------------------------------
## Oxique
totalOxiqueASV <- data_asv_tax %>% select(ASV_Id,TotalOxique,Division)
row.names(totalOxiqueASV)<-totalOxiqueASV$ASV_Id ; totalOxiqueASV <- totalOxiqueASV %>% select(-ASV_Id)
totalOxiqueASV <- totalOxiqueASV %>% group_by(Division) %>% summarise_all(sum)
totalOxiqueASV$TotalOxique <- totalOxiqueASV$TotalOxique*100/sum(totalOxiqueASV$TotalOxique)
totalOxiqueASV <- totalOxiqueASV %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalOxique) - 0.5*TotalOxique)
totalOxiqueASV$label <- paste(round(totalOxiqueASV$TotalOxique,1), "%", sep = "")
for (i in rownames(totalOxiqueASV)) {
  if (totalOxiqueASV[i,"label"] == "0%") { totalOxiqueASV[i,"label"] <- NA}}
for (i in rownames(totalOxiqueASV)) {
  if (is.na(totalOxiqueASV[i,"label"]) == FALSE) { totalOxiqueASV[i,"label"] <- paste(totalOxiqueASV[i,"Division"]," : ",totalOxiqueASV[i,"label"], sep = "")}}
totalOxiqueASV$Zone<- rep("Oxique", each = nrow(totalOxiqueASV))
## Anoxique
totalAnoxiqueASV <- data_asv_tax %>% select(ASV_Id,TotalAnoxique,Division)
row.names(totalAnoxiqueASV)<-totalAnoxiqueASV$ASV_Id ; totalAnoxiqueASV <- totalAnoxiqueASV %>% select(-ASV_Id)
totalAnoxiqueASV <- totalAnoxiqueASV %>% group_by(Division) %>% summarise_all(sum)
totalAnoxiqueASV$TotalAnoxique <- totalAnoxiqueASV$TotalAnoxique*100/sum(totalAnoxiqueASV$TotalAnoxique)
totalAnoxiqueASV <- totalAnoxiqueASV %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalAnoxique) - 0.5*TotalAnoxique)
totalAnoxiqueASV$label <- paste(round(totalAnoxiqueASV$TotalAnoxique,1), "%", sep = "")
for (i in rownames(totalAnoxiqueASV)) {
  if (totalAnoxiqueASV[i,"label"] == "0%") { totalAnoxiqueASV[i,"label"] <- NA}}
for (i in rownames(totalAnoxiqueASV)) {
  if (is.na(totalAnoxiqueASV[i,"label"]) == FALSE) { totalAnoxiqueASV[i,"label"] <- paste(totalAnoxiqueASV[i,"Division"]," : ",totalAnoxiqueASV[i,"label"], sep = "")}}
totalAnoxiqueASV$Zone<- rep("Anoxique", each = nrow(totalAnoxiqueASV))
## Zone
colnames(totalOxiqueASV)[2]  <- "value"
totalOxiqueASV$Sum <- rep(0, each = nrow(totalOxiqueASV))
totalOxiqueASV$Count <- rep(0, each = nrow(totalOxiqueASV))
for (i in totalOxiqueASV$Division) { totalOxiqueASV$Count[which(totalOxiqueASV$Division == i)] <- sum(data_asv_tax  %>% filter(Division == i) %>% select(TotalOxique))}
totalOxiqueASV$Sum <- sum(totalOxiqueASV$Count)
colnames(totalAnoxiqueASV)[2]  <- "value"
totalAnoxiqueASV$Sum <- rep(0, each = nrow(totalAnoxiqueASV))
totalAnoxiqueASV$Count <- rep(0, each = nrow(totalAnoxiqueASV))
for (i in totalAnoxiqueASV$Division) { totalAnoxiqueASV$Count[which(totalAnoxiqueASV$Division == i)] <- sum(data_asv_tax  %>% filter(Division == i) %>% select(TotalAnoxique))}
totalAnoxiqueASV$Sum <- sum(totalAnoxiqueASV$Count)
totalZoneASV <- rbind(totalOxiqueASV,totalAnoxiqueASV)
#Figure
by <- ggplot(totalZoneASV, mapping = aes(y= value, x = Zone, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + scale_fill_manual(values = palette) + 
  geom_label(aes(y = 106,label = paste(Sum," ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
by <- by + labs(x="Zones",y="ASVs (%)") + theme(legend.position = "none")
print(by)

    # Fraction -------------------------------------------------------------------
## Petite
totalPetiteASV <- data_asv_tax %>% select(ASV_Id,TotalPetite,Division)
row.names(totalPetiteASV)<-totalPetiteASV$ASV_Id ; totalPetiteASV <- totalPetiteASV %>% select(-ASV_Id)
totalPetiteASV <- totalPetiteASV %>% group_by(Division) %>% summarise_all(sum)
totalPetiteASV$TotalPetite <- totalPetiteASV$TotalPetite*100/sum(totalPetiteASV$TotalPetite)
totalPetiteASV <- totalPetiteASV %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalPetite) - 0.5*TotalPetite)
totalPetiteASV$label <- paste(round(totalPetiteASV$TotalPetite,1), "%", sep = "")
for (i in rownames(totalPetiteASV)) {
  if (totalPetiteASV[i,"label"] == "0%") { totalPetiteASV[i,"label"] <- NA}}
for (i in rownames(totalPetiteASV)) {
  if (is.na(totalPetiteASV[i,"label"]) == FALSE) { totalPetiteASV[i,"label"] <- paste(totalPetiteASV[i,"Division"]," : ",totalPetiteASV[i,"label"], sep = "")}}
totalPetiteASV$Fraction<- rep("Petite", each = nrow(totalPetiteASV))
## Grande
totalGrandeASV <- data_asv_tax %>% select(ASV_Id,TotalGrande,Division)
row.names(totalGrandeASV)<-totalGrandeASV$ASV_Id ; totalGrandeASV <- totalGrandeASV %>% select(-ASV_Id)
totalGrandeASV <- totalGrandeASV %>% group_by(Division) %>% summarise_all(sum)
totalGrandeASV$TotalGrande <- totalGrandeASV$TotalGrande*100/sum(totalGrandeASV$TotalGrande)
totalGrandeASV <- totalGrandeASV %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalGrande) - 0.5*TotalGrande)
totalGrandeASV$label <- paste(round(totalGrandeASV$TotalGrande,1), "%", sep = "")
for (i in rownames(totalGrandeASV)) {
  if (totalGrandeASV[i,"label"] == "0%") { totalGrandeASV[i,"label"] <- NA}}
for (i in rownames(totalGrandeASV)) {
  if (is.na(totalGrandeASV[i,"label"]) == FALSE) { totalGrandeASV[i,"label"] <- paste(totalGrandeASV[i,"Division"]," : ",totalGrandeASV[i,"label"], sep = "")}}
totalGrandeASV$Fraction<- rep("Grande", each = nrow(totalGrandeASV))
## Fraction
colnames(totalPetiteASV)[2]  <- "value"
totalPetiteASV$Sum <- rep(0, each = nrow(totalPetiteASV))
totalPetiteASV$Count <- rep(0, each = nrow(totalPetiteASV))
for (i in totalPetiteASV$Division) { totalPetiteASV$Count[which(totalPetiteASV$Division == i)] <- sum(data_asv_tax  %>% filter(Division == i) %>% select(TotalPetite))}
totalPetiteASV$Sum <- sum(totalPetiteASV$Count)
colnames(totalGrandeASV)[2]  <- "value"
totalGrandeASV$Sum <- rep(0, each = nrow(totalGrandeASV))
totalGrandeASV$Count <- rep(0, each = nrow(totalGrandeASV))
for (i in totalGrandeASV$Division) { totalGrandeASV$Count[which(totalGrandeASV$Division == i)] <- sum(data_asv_tax  %>% filter(Division == i) %>% select(TotalGrande))}
totalGrandeASV$Sum <- sum(totalGrandeASV$Count)
totalFractionASV <- rbind(totalPetiteASV,totalGrandeASV)
#Figure
cy <- ggplot(totalFractionASV, mapping = aes(y= value, x = Fraction, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + scale_fill_manual(values = palette)+ 
  geom_label(aes(y = 106,label = paste(Sum," ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
cy <- cy + labs(x="Fractions",y="ASVs (%)") + theme(legend.position = "none")
print(cy) 
    # Date -------------------------------------------------------------------
## Avril
total04ASV <- data_asv_tax %>% select(ASV_Id,Total04,Division)
row.names(total04ASV)<-total04ASV$ASV_Id ; total04ASV <- total04ASV %>% select(-ASV_Id)
total04ASV <- total04ASV %>% group_by(Division) %>% summarise_all(sum)
total04ASV$Total04 <- total04ASV$Total04*100/sum(total04ASV$Total04)
total04ASV <- total04ASV %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(Total04) - 0.5*Total04)
total04ASV$label <- paste(round(total04ASV$Total04,1), "%", sep = "")
for (i in rownames(total04ASV)) {
  if (total04ASV[i,"label"] == "0%") { total04ASV[i,"label"] <- NA}}
for (i in rownames(total04ASV)) {
  if (is.na(total04ASV[i,"label"]) == FALSE) { total04ASV[i,"label"] <- paste(total04ASV[i,"Division"]," : ",total04ASV[i,"label"], sep = "")}}
total04ASV$Date<- rep("04", each = nrow(total04ASV))
## Juin
total06ASV <- data_asv_tax %>% select(ASV_Id,Total06,Division)
row.names(total06ASV)<-total06ASV$ASV_Id ; total06ASV <- total06ASV %>% select(-ASV_Id)
total06ASV <- total06ASV %>% group_by(Division) %>% summarise_all(sum)
total06ASV$Total06 <- total06ASV$Total06*100/sum(total06ASV$Total06)
total06ASV <- total06ASV %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(Total06) - 0.5*Total06)
total06ASV$label <- paste(round(total06ASV$Total06,1), "%", sep = "")
for (i in rownames(total06ASV)) {
  if (total06ASV[i,"label"] == "0%") { total06ASV[i,"label"] <- NA}}
for (i in rownames(total06ASV)) {
  if (is.na(total06ASV[i,"label"]) == FALSE) { total06ASV[i,"label"] <- paste(total06ASV[i,"Division"]," : ",total06ASV[i,"label"], sep = "")}}
total06ASV$Date<- rep("06", each = nrow(total06ASV))
## Septembre
total09ASV <- data_asv_tax %>% select(ASV_Id,Total09,Division)
row.names(total09ASV)<-total09ASV$ASV_Id ; total09ASV <- total09ASV %>% select(-ASV_Id)
total09ASV <- total09ASV %>% group_by(Division) %>% summarise_all(sum)
total09ASV$Total09 <- total09ASV$Total09*100/sum(total09ASV$Total09)
total09ASV <- total09ASV %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(Total09) - 0.5*Total09)
total09ASV$label <- paste(round(total09ASV$Total09,1), "%", sep = "")
for (i in rownames(total09ASV)) {
  if (total09ASV[i,"label"] == "0%") { total09ASV[i,"label"] <- NA}}
for (i in rownames(total09ASV)) {
  if (is.na(total09ASV[i,"label"]) == FALSE) { total09ASV[i,"label"] <- paste(total09ASV[i,"Division"]," : ",total09ASV[i,"label"], sep = "")}}
total09ASV$Date<- rep("09", each = nrow(total09ASV))
## Octobre
total11ASV <- data_asv_tax %>% select(ASV_Id,Total11,Division)
row.names(total11ASV)<-total11ASV$ASV_Id ; total11ASV <- total11ASV %>% select(-ASV_Id)
total11ASV <- total11ASV %>% group_by(Division) %>% summarise_all(sum)
total11ASV$Total11 <- total11ASV$Total11*100/sum(total11ASV$Total11)
total11ASV <- total11ASV %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(Total11) - 0.5*Total11)
total11ASV$label <- paste(round(total11ASV$Total11,1), "%", sep = "")
for (i in rownames(total11ASV)) {
  if (total11ASV[i,"label"] == "0%") { total11ASV[i,"label"] <- NA}}
for (i in rownames(total11ASV)) {
  if (is.na(total11ASV[i,"label"]) == FALSE) { total11ASV[i,"label"] <- paste(total11ASV[i,"Division"]," : ",total11ASV[i,"label"], sep = "")}}
total11ASV$Date<- rep("11", each = nrow(total11ASV))
## Date
colnames(total04ASV)[2]  <- "value"
total04ASV$Sum <- rep(0, each = nrow(total04ASV))
total04ASV$Count <- rep(0, each = nrow(total04ASV))
for (i in total04ASV$Division) { total04ASV$Count[which(total04ASV$Division == i)] <- sum(data_asv_tax  %>% filter(Division == i) %>% select(Total04))}
total04ASV$Sum <- sum(total04ASV$Count)
colnames(total06ASV)[2]  <- "value"
total06ASV$Sum <- rep(0, each = nrow(total06ASV))
total06ASV$Count <- rep(0, each = nrow(total06ASV))
for (i in total06ASV$Division) { total06ASV$Count[which(total06ASV$Division == i)] <- sum(data_asv_tax  %>% filter(Division == i) %>% select(Total06))}
total06ASV$Sum <- sum(total06ASV$Count)
colnames(total09ASV)[2]  <- "value"
total09ASV$Sum <- rep(0, each = nrow(total09ASV))
total09ASV$Count <- rep(0, each = nrow(total09ASV))
for (i in total09ASV$Division) { total09ASV$Count[which(total09ASV$Division == i)] <- sum(data_asv_tax  %>% filter(Division == i) %>% select(Total09))}
total09ASV$Sum <- sum(total09ASV$Count)
colnames(total11ASV)[2]  <- "value"
total11ASV$Sum <- rep(0, each = nrow(total11ASV))
total11ASV$Count <- rep(0, each = nrow(total11ASV))
for (i in total11ASV$Division) { total11ASV$Count[which(total11ASV$Division == i)] <- sum(data_asv_tax  %>% filter(Division == i) %>% select(Total11))}
total11ASV$Sum <- sum(total11ASV$Count)
totalDateASV <- rbind(total04ASV,total06ASV,total09ASV,total11ASV)
#Figure
dy <- ggplot(totalDateASV, mapping = aes(y= value, x = Date, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + scale_fill_manual(values = palette)+ 
  geom_label(aes(y = 106,label = paste(Sum," ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
dy <- dy + theme(legend.position = "none") + labs(x="Dates",y="ASVs (%)")
print(dy) 
    # Coplot -------------------------------------------------------------------
svglite("Hist-Taxomy/ASV-Total.svg",width = 12.00,height = 6.00)
b_plot <- plot_grid(ay,by,cy,dy,legendASV, ncol = 5, nrow = 1, rel_widths = c(3,3,3,5,3),rel_heights = c(3))
print(b_plot)
dev.off()
  # ASV ONLY ---------------------------------------------------------------------
    # Cycles -------------------------------------------------------------------
## Jour
onlyJourASV <- data_asv_tax %>% filter(Cycle == "Jour") %>% select(ASV_Id,TotalJour,Division)
row.names(onlyJourASV) <- onlyJourASV$ASV_Id ; onlyJourASV <- onlyJourASV %>% select(-ASV_Id)
onlyJourASV <- onlyJourASV %>% group_by(Division) %>% summarise_all(sum)
onlyJourASV$TotalJour <- onlyJourASV$TotalJour*100/sum(onlyJourASV$TotalJour)
onlyJourASV <- onlyJourASV %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalJour) - 0.5*TotalJour)
onlyJourASV$label <- paste(round(onlyJourASV$TotalJour,1), "%", sep = "")
for (i in rownames(onlyJourASV)) {
  if (onlyJourASV[i,"label"] == "0%") { onlyJourASV[i,"label"] <- NA}}
for (i in rownames(onlyJourASV)) {
  if (is.na(onlyJourASV[i,"label"]) == FALSE) { onlyJourASV[i,"label"] <- paste(onlyJourASV[i,"Division"]," : ",onlyJourASV[i,"label"], sep = "")}}
onlyJourASV$Cycle<- rep("Jour", each = nrow(onlyJourASV))
## Nuit
onlyNuitASV <- data_asv_tax %>% filter(Cycle == "Nuit") %>% select(ASV_Id,TotalNuit,Division)
row.names(onlyNuitASV)<-onlyNuitASV$ASV_Id ; onlyNuitASV <- onlyNuitASV %>% select(-ASV_Id)
onlyNuitASV <- onlyNuitASV %>% group_by(Division) %>% summarise_all(sum)
onlyNuitASV$TotalNuit <- onlyNuitASV$TotalNuit*100/sum(onlyNuitASV$TotalNuit)
onlyNuitASV <- onlyNuitASV %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalNuit) - 0.5*TotalNuit)
onlyNuitASV$label <- paste(round(onlyNuitASV$TotalNuit,1), "%", sep = "")
for (i in rownames(onlyNuitASV)) {
  if (onlyNuitASV[i,"label"] == "0%") { onlyNuitASV[i,"label"] <- NA}}
for (i in rownames(onlyNuitASV)) {
  if (is.na(onlyNuitASV[i,"label"]) == FALSE) { onlyNuitASV[i,"label"] <- paste(onlyNuitASV[i,"Division"]," : ",onlyNuitASV[i,"label"], sep = "")}}
onlyNuitASV$Cycle<- rep("Nuit", each = nrow(onlyNuitASV))
## Cycle
colnames(onlyJourASV)[2]  <- "value"
onlyJourASV$Sum <- rep(0, each = nrow(onlyJourASV))
onlyJourASV$Count <- rep(0, each = nrow(onlyJourASV))
for (i in onlyJourASV$Division) { onlyJourASV$Count[which(onlyJourASV$Division == i)] <- nrow(data_asv_tax %>% filter(Cycle == "Jour") %>% filter(Division == i))}
onlyJourASV$Sum <- sum(onlyJourASV$Count)
colnames(onlyNuitASV)[2]  <- "value"
onlyNuitASV$Sum <- rep(0, each = nrow(onlyNuitASV))
onlyNuitASV$Count <- rep(0, each = nrow(onlyNuitASV))
for (i in onlyNuitASV$Division) { onlyNuitASV$Count[which(onlyNuitASV$Division == i)] <- nrow(data_asv_tax %>% filter(Cycle == "Nuit") %>% filter(Division == i))}
onlyNuitASV$Sum <- sum(onlyNuitASV$Count)
onlyCycleASV <- rbind(onlyJourASV,onlyNuitASV)
#Figure
az <- ggplot(onlyCycleASV, mapping = aes(y= value, x = Cycle, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
  geom_label(aes(y = 106,label = paste(Sum," ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") + scale_fill_manual(values = palette)
#legendASV <- get_legend(az)
az <- az + labs(x="Cycles",y="ASVs (%)") #+ theme(legend.position = "none")
print(az)

    # Zone -------------------------------------------------------------------
## Oxique
onlyOxiqueASV <- data_asv_tax %>% filter(Zone == "Oxique") %>% select(ASV_Id,TotalOxique,Division)
row.names(onlyOxiqueASV)<-onlyOxiqueASV$ASV_Id ; onlyOxiqueASV <- onlyOxiqueASV %>% select(-ASV_Id)
onlyOxiqueASV <- onlyOxiqueASV %>% group_by(Division) %>% summarise_all(sum)
onlyOxiqueASV$TotalOxique <- onlyOxiqueASV$TotalOxique*100/sum(onlyOxiqueASV$TotalOxique)
onlyOxiqueASV <- onlyOxiqueASV %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalOxique) - 0.5*TotalOxique)
onlyOxiqueASV$label <- paste(round(onlyOxiqueASV$TotalOxique,1), "%", sep = "")
for (i in rownames(onlyOxiqueASV)) {
  if (onlyOxiqueASV[i,"label"] == "0%") { onlyOxiqueASV[i,"label"] <- NA}}
for (i in rownames(onlyOxiqueASV)) {
  if (is.na(onlyOxiqueASV[i,"label"]) == FALSE) { onlyOxiqueASV[i,"label"] <- paste(onlyOxiqueASV[i,"Division"]," : ",onlyOxiqueASV[i,"label"], sep = "")}}
onlyOxiqueASV$Zone<- rep("Oxique", each = nrow(onlyOxiqueASV))
## Anoxique
onlyAnoxiqueASV <- data_asv_tax %>% filter(Zone == "Anoxique") %>% select(ASV_Id,TotalAnoxique,Division)
row.names(onlyAnoxiqueASV)<-onlyAnoxiqueASV$ASV_Id ; onlyAnoxiqueASV <- onlyAnoxiqueASV %>% select(-ASV_Id)
onlyAnoxiqueASV <- onlyAnoxiqueASV %>% group_by(Division) %>% summarise_all(sum)
onlyAnoxiqueASV$TotalAnoxique <- onlyAnoxiqueASV$TotalAnoxique*100/sum(onlyAnoxiqueASV$TotalAnoxique)
onlyAnoxiqueASV <- onlyAnoxiqueASV %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalAnoxique) - 0.5*TotalAnoxique)
onlyAnoxiqueASV$label <- paste(round(onlyAnoxiqueASV$TotalAnoxique,1), "%", sep = "")
for (i in rownames(onlyAnoxiqueASV)) {
  if (onlyAnoxiqueASV[i,"label"] == "0%") { onlyAnoxiqueASV[i,"label"] <- NA}}
for (i in rownames(onlyAnoxiqueASV)) {
  if (is.na(onlyAnoxiqueASV[i,"label"]) == FALSE) { onlyAnoxiqueASV[i,"label"] <- paste(onlyAnoxiqueASV[i,"Division"]," : ",onlyAnoxiqueASV[i,"label"], sep = "")}}
onlyAnoxiqueASV$Zone<- rep("Anoxique", each = nrow(onlyAnoxiqueASV))
## Zone
colnames(onlyOxiqueASV)[2]  <- "value"
onlyOxiqueASV$Sum <- rep(0, each = nrow(onlyOxiqueASV))
onlyOxiqueASV$Count <- rep(0, each = nrow(onlyOxiqueASV))
for (i in onlyOxiqueASV$Division) { onlyOxiqueASV$Count[which(onlyOxiqueASV$Division == i)] <- nrow(data_asv_tax %>% filter(Zone == "Oxique") %>% filter(Division == i))}
onlyOxiqueASV$Sum <- sum(onlyOxiqueASV$Count)
colnames(onlyAnoxiqueASV)[2]  <- "value"
onlyAnoxiqueASV$Sum <- rep(0, each = nrow(onlyAnoxiqueASV))
onlyAnoxiqueASV$Count <- rep(0, each = nrow(onlyAnoxiqueASV))
for (i in onlyAnoxiqueASV$Division) { onlyAnoxiqueASV$Count[which(onlyAnoxiqueASV$Division == i)] <- nrow(data_asv_tax %>% filter(Zone == "Anoxique") %>% filter(Division == i))}
onlyAnoxiqueASV$Sum <- sum(onlyAnoxiqueASV$Count)
onlyZoneASV <- rbind(onlyOxiqueASV,onlyAnoxiqueASV)
#Figure
bz <- ggplot(onlyZoneASV, mapping = aes(y= value, x = Zone, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
  geom_label(aes(y = 106,label = paste(Sum," ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") + scale_fill_manual(values = palette)
bz <- bz + labs(x="Zones",y="ASVs (%)") #+ theme(legend.position = "none")
print(bz)

    # Fraction -------------------------------------------------------------------
## Petite
onlyPetiteASV <- data_asv_tax %>% filter(Fraction == "Petite") %>% select(ASV_Id,TotalPetite,Division)
row.names(onlyPetiteASV)<-onlyPetiteASV$ASV_Id ; onlyPetiteASV <- onlyPetiteASV %>% select(-ASV_Id)
onlyPetiteASV <- onlyPetiteASV %>% group_by(Division) %>% summarise_all(sum)
onlyPetiteASV$TotalPetite <- onlyPetiteASV$TotalPetite*100/sum(onlyPetiteASV$TotalPetite)
onlyPetiteASV <- onlyPetiteASV %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalPetite) - 0.5*TotalPetite)
onlyPetiteASV$label <- paste(round(onlyPetiteASV$TotalPetite,1), "%", sep = "")
for (i in rownames(onlyPetiteASV)) {
  if (onlyPetiteASV[i,"label"] == "0%") { onlyPetiteASV[i,"label"] <- NA}}
for (i in rownames(onlyPetiteASV)) {
  if (is.na(onlyPetiteASV[i,"label"]) == FALSE) { onlyPetiteASV[i,"label"] <- paste(onlyPetiteASV[i,"Division"]," : ",onlyPetiteASV[i,"label"], sep = "")}}
onlyPetiteASV$Fraction<- rep("Petite", each = nrow(onlyPetiteASV))
## Grande
onlyGrandeASV <- data_asv_tax %>% filter(Fraction == "Grande") %>% select(ASV_Id,TotalGrande,Division)
row.names(onlyGrandeASV)<-onlyGrandeASV$ASV_Id
onlyGrandeASV <- onlyGrandeASV %>% select(-ASV_Id)
onlyGrandeASV <- onlyGrandeASV %>% group_by(Division) %>% summarise_all(sum)
onlyGrandeASV$TotalGrande <- onlyGrandeASV$TotalGrande*100/sum(onlyGrandeASV$TotalGrande)
onlyGrandeASV <- onlyGrandeASV %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalGrande) - 0.5*TotalGrande)
onlyGrandeASV$label <- paste(round(onlyGrandeASV$TotalGrande,1), "%", sep = "")
for (i in rownames(onlyGrandeASV)) {
  if (onlyGrandeASV[i,"label"] == "0%") { onlyGrandeASV[i,"label"] <- NA}}
for (i in rownames(onlyGrandeASV)) {
  if (is.na(onlyGrandeASV[i,"label"]) == FALSE) { onlyGrandeASV[i,"label"] <- paste(onlyGrandeASV[i,"Division"]," : ",onlyGrandeASV[i,"label"], sep = "")}}
onlyGrandeASV$Fraction<- rep("Grande", each = nrow(onlyGrandeASV))
## Fraction
colnames(onlyPetiteASV)[2]  <- "value"
onlyPetiteASV$Sum <- rep(0, each = nrow(onlyPetiteASV))
onlyPetiteASV$Count <- rep(0, each = nrow(onlyPetiteASV))
for (i in onlyPetiteASV$Division) { onlyPetiteASV$Count[which(onlyPetiteASV$Division == i)] <- nrow(data_asv_tax %>% filter(Fraction == "Petite") %>% filter(Division == i))}
onlyPetiteASV$Sum <- sum(onlyPetiteASV$Count)
colnames(onlyGrandeASV)[2]  <- "value"
onlyGrandeASV$Sum <- rep(0, each = nrow(onlyGrandeASV))
onlyGrandeASV$Count <- rep(0, each = nrow(onlyGrandeASV))
for (i in onlyGrandeASV$Division) { onlyGrandeASV$Count[which(onlyGrandeASV$Division == i)] <- nrow(data_asv_tax %>% filter(Fraction == "Grande") %>% filter(Division == i))}
onlyGrandeASV$Sum <- sum(onlyGrandeASV$Count)
onlyFractionASV <- rbind(onlyPetiteASV,onlyGrandeASV)
#Figure
cz <- ggplot(onlyFractionASV, mapping = aes(y= value, x = Fraction, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
  geom_label(aes(y = 106,label = paste(Sum," ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
  scale_fill_manual(values = palette)
cz <- cz + labs(x="Fractions",y="ASVs (%)") #+ theme(legend.position = "none")
print(cz) 
    # Coplot -------------------------------------------------------------------
svglite("Hist-Taxomy/ASV-only.svg",width = 12.00,height = 6.00)
b_plot <- plot_grid(az,bz,cz,  ncol = 3, nrow = 1, rel_widths = c(3,3,3),rel_heights = c(3))
print(b_plot)
dev.off()  
    # Coplot X AFC -------------------------------------------------------------------
f <- ggplot() + theme_void()
print(f)
#AFC = séquence
svglite("Biplot/ASV-onlyX-AFCseq.svg",width = 12.00,height = 14.00)
b_plot <- plot_grid(aCycle,az,bZone,bz,cFraction,cz, ncol = 2, nrow = 3, rel_widths = c(3,2),rel_heights = c(3,3,3))
print(b_plot)
dev.off()  
#AFC = ASV
svglite("Biplot/ASV-onlyX-AFCasv.svg",width = 12.00,height = 14.00)
c_plot <- plot_grid(iCycle,az,jZone,bz,kFraction,cz, ncol = 2, nrow = 3, rel_widths = c(3,2),rel_heights = c(3,3,3))
print(c_plot)
dev.off()  



  # Séquence TOTAL---------------------------------------------------------------------
    # Cycle -------------------------------------------------------------------
## Jour
totalJourSequence <- data_seq_tax %>% select(ASV_Id,TotalJour,Division)
row.names(totalJourSequence)<-totalJourSequence$ASV_Id ; totalJourSequence <- totalJourSequence %>% select(-ASV_Id)
totalJourSequence <- totalJourSequence %>% group_by(Division) %>% summarise_all(sum)
totalJourSequence$TotalJour <- totalJourSequence$TotalJour*100/sum(totalJourSequence$TotalJour)
totalJourSequence <- totalJourSequence %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalJour) - 0.5*TotalJour)
totalJourSequence$label <- paste(round(totalJourSequence$TotalJour,1), "%", sep = "")
for (i in rownames(totalJourSequence)) {
  if (totalJourSequence[i,"label"] == "0%") { totalJourSequence[i,"label"] <- NA}}
for (i in rownames(totalJourSequence)) {
  if (is.na(totalJourSequence[i,"label"]) == FALSE) { totalJourSequence[i,"label"] <- paste(totalJourSequence[i,"Division"]," : ",totalJourSequence[i,"label"], sep = "")}}
totalJourSequence$Cycle<- rep("Jour", each = nrow(totalJourSequence))
## Nuit
totalNuitSequence <- data_seq_tax %>% select(ASV_Id,TotalNuit,Division)
row.names(totalNuitSequence)<-totalNuitSequence$ASV_Id ; totalNuitSequence <- totalNuitSequence %>% select(-ASV_Id)
totalNuitSequence <- totalNuitSequence %>% group_by(Division) %>% summarise_all(sum)
totalNuitSequence$TotalNuit <- totalNuitSequence$TotalNuit*100/sum(totalNuitSequence$TotalNuit)
totalNuitSequence <- totalNuitSequence %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalNuit) - 0.5*TotalNuit)
totalNuitSequence$label <- paste(round(totalNuitSequence$TotalNuit,1), "%", sep = "")
for (i in rownames(totalNuitSequence)) {
  if (totalNuitSequence[i,"label"] == "0%") { totalNuitSequence[i,"label"] <- NA}}
for (i in rownames(totalNuitSequence)) {
  if (is.na(totalNuitSequence[i,"label"]) == FALSE) { totalNuitSequence[i,"label"] <- paste(totalNuitSequence[i,"Division"]," : ",totalNuitSequence[i,"label"], sep = "")}}
totalNuitSequence$Cycle<- rep("Nuit", each = nrow(totalNuitSequence))
## Cycle
colnames(totalJourSequence)[2]  <- "value"
totalJourSequence$Sum <- rep(0, each = nrow(totalJourSequence))
totalJourSequence$Count <- rep(0, each = nrow(totalJourSequence))
for (i in totalJourSequence$Division) { totalJourSequence$Count[which(totalJourSequence$Division == i)] <- sum(data_seq_tax  %>% filter(Division == i) %>% select(TotalJour))}
totalJourSequence$Sum <- sum(totalJourSequence$Count)
colnames(totalNuitSequence)[2]  <- "value"
totalNuitSequence$Sum <- rep(0, each = nrow(totalNuitSequence))
totalNuitSequence$Count <- rep(0, each = nrow(totalNuitSequence))
for (i in totalNuitSequence$Division) { totalNuitSequence$Count[which(totalNuitSequence$Division == i)] <- sum(data_seq_tax  %>% filter(Division == i) %>% select(TotalNuit))}
totalNuitSequence$Sum <- sum(totalNuitSequence$Count)
totalCycleSequence <- rbind(totalJourSequence,totalNuitSequence)
totalCycleSequence$percent <- paste("(",round(totalCycleSequence$Sum*100/(colSums(statRarefy %>% select(`apRarefy-Sequence`))/2),1)," %)",sep ="")

#Figure
iy <- ggplot(totalCycleSequence, mapping = aes(y= value, x = Cycle, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
  scale_fill_manual(values = palette) + 
  geom_label(aes(y = 108,label = paste(Sum,"séquences",percent,sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
legendSequence <- get_legend(iy)
iy <- iy + labs(x="Cycles",y="Séquences (%)") + theme(legend.position = "none")
print(iy)

    # Zone -------------------------------------------------------------------
## Oxique
totalOxiqueSequence <- data_seq_tax %>% select(ASV_Id,TotalOxique,Division)
row.names(totalOxiqueSequence)<-totalOxiqueSequence$ASV_Id ; totalOxiqueSequence <- totalOxiqueSequence %>% select(-ASV_Id)
totalOxiqueSequence <- totalOxiqueSequence %>% group_by(Division) %>% summarise_all(sum)
totalOxiqueSequence$TotalOxique <- totalOxiqueSequence$TotalOxique*100/sum(totalOxiqueSequence$TotalOxique)
totalOxiqueSequence <- totalOxiqueSequence %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalOxique) - 0.5*TotalOxique)
totalOxiqueSequence$label <- paste(round(totalOxiqueSequence$TotalOxique,1), "%", sep = "")
for (i in rownames(totalOxiqueSequence)) {
  if (totalOxiqueSequence[i,"label"] == "0%") { totalOxiqueSequence[i,"label"] <- NA}}
for (i in rownames(totalOxiqueSequence)) {
  if (is.na(totalOxiqueSequence[i,"label"]) == FALSE) { totalOxiqueSequence[i,"label"] <- paste(totalOxiqueSequence[i,"Division"]," : ",totalOxiqueSequence[i,"label"], sep = "")}}
totalOxiqueSequence$Zone<- rep("Oxique", each = nrow(totalOxiqueSequence))
## Anoxique
totalAnoxiqueSequence <- data_seq_tax %>% select(ASV_Id,TotalAnoxique,Division)
row.names(totalAnoxiqueSequence)<-totalAnoxiqueSequence$ASV_Id ; totalAnoxiqueSequence <- totalAnoxiqueSequence %>% select(-ASV_Id)
totalAnoxiqueSequence <- totalAnoxiqueSequence %>% group_by(Division) %>% summarise_all(sum)
totalAnoxiqueSequence$TotalAnoxique <- totalAnoxiqueSequence$TotalAnoxique*100/sum(totalAnoxiqueSequence$TotalAnoxique)
totalAnoxiqueSequence <- totalAnoxiqueSequence %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalAnoxique) - 0.5*TotalAnoxique)
totalAnoxiqueSequence$label <- paste(round(totalAnoxiqueSequence$TotalAnoxique,1), "%", sep = "")
for (i in rownames(totalAnoxiqueSequence)) {
  if (totalAnoxiqueSequence[i,"label"] == "0%") { totalAnoxiqueSequence[i,"label"] <- NA}}
for (i in rownames(totalAnoxiqueSequence)) {
  if (is.na(totalAnoxiqueSequence[i,"label"]) == FALSE) { totalAnoxiqueSequence[i,"label"] <- paste(totalAnoxiqueSequence[i,"Division"]," : ",totalAnoxiqueSequence[i,"label"], sep = "")}}
totalAnoxiqueSequence$Zone<- rep("Anoxique", each = nrow(totalAnoxiqueSequence))
## Zone
colnames(totalOxiqueSequence)[2]  <- "value"
totalOxiqueSequence$Sum <- rep(0, each = nrow(totalOxiqueSequence))
totalOxiqueSequence$Count <- rep(0, each = nrow(totalOxiqueSequence))
for (i in totalOxiqueSequence$Division) { totalOxiqueSequence$Count[which(totalOxiqueSequence$Division == i)] <- sum(data_seq_tax  %>% filter(Division == i) %>% select(TotalOxique))}
totalOxiqueSequence$Sum <- sum(totalOxiqueSequence$Count)
colnames(totalAnoxiqueSequence)[2]  <- "value"
totalAnoxiqueSequence$Sum <- rep(0, each = nrow(totalAnoxiqueSequence))
totalAnoxiqueSequence$Count <- rep(0, each = nrow(totalAnoxiqueSequence))
for (i in totalAnoxiqueSequence$Division) { totalAnoxiqueSequence$Count[which(totalAnoxiqueSequence$Division == i)] <- sum(data_seq_tax  %>% filter(Division == i) %>% select(TotalAnoxique))}
totalAnoxiqueSequence$Sum <- sum(totalAnoxiqueSequence$Count)
totalZoneSequence <- rbind(totalOxiqueSequence,totalAnoxiqueSequence)
totalZoneSequence$percent <- paste("(",round(totalZoneSequence$Sum*100/(colSums(statRarefy %>% select(`apRarefy-Sequence`))/2),1)," %)",sep ="")

#Figure
jy <- ggplot(totalZoneSequence, mapping = aes(y= value, x = Zone, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
  scale_fill_manual(values = palette) + 
  geom_label(aes(y = 108,label = paste(Sum,"séquences",percent,sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
jy <- jy + labs(x="Zones",y="Séquences (%)") + theme(legend.position = "none")
print(jy)    
    # Fraction -------------------------------------------------------------------
## Petite
totalPetiteSequence <- data_seq_tax %>% select(ASV_Id,TotalPetite,Division)
row.names(totalPetiteSequence)<-totalPetiteSequence$ASV_Id ; totalPetiteSequence <- totalPetiteSequence %>% select(-ASV_Id)
totalPetiteSequence <- totalPetiteSequence %>% group_by(Division) %>% summarise_all(sum)
totalPetiteSequence$TotalPetite <- totalPetiteSequence$TotalPetite*100/sum(totalPetiteSequence$TotalPetite)
totalPetiteSequence <- totalPetiteSequence %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalPetite) - 0.5*TotalPetite)
totalPetiteSequence$label <- paste(round(totalPetiteSequence$TotalPetite,1), "%", sep = "")
for (i in rownames(totalPetiteSequence)) {
  if (totalPetiteSequence[i,"label"] == "0%") { totalPetiteSequence[i,"label"] <- NA}}
for (i in rownames(totalPetiteSequence)) {
  if (is.na(totalPetiteSequence[i,"label"]) == FALSE) { totalPetiteSequence[i,"label"] <- paste(totalPetiteSequence[i,"Division"]," : ",totalPetiteSequence[i,"label"], sep = "")}}
totalPetiteSequence$Fraction<- rep("Petite", each = nrow(totalPetiteSequence))
## Grande
totalGrandeSequence <- data_seq_tax %>% select(ASV_Id,TotalGrande,Division)
row.names(totalGrandeSequence)<-totalGrandeSequence$ASV_Id ; totalGrandeSequence <- totalGrandeSequence %>% select(-ASV_Id)
totalGrandeSequence <- totalGrandeSequence %>% group_by(Division) %>% summarise_all(sum)
totalGrandeSequence$TotalGrande <- totalGrandeSequence$TotalGrande*100/sum(totalGrandeSequence$TotalGrande)
totalGrandeSequence <- totalGrandeSequence %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalGrande) - 0.5*TotalGrande)
totalGrandeSequence$label <- paste(round(totalGrandeSequence$TotalGrande,1), "%", sep = "")
for (i in rownames(totalGrandeSequence)) {
  if (totalGrandeSequence[i,"label"] == "0%") { totalGrandeSequence[i,"label"] <- NA}}
for (i in rownames(totalGrandeSequence)) {
  if (is.na(totalGrandeSequence[i,"label"]) == FALSE) { totalGrandeSequence[i,"label"] <- paste(totalGrandeSequence[i,"Division"]," : ",totalGrandeSequence[i,"label"], sep = "")}}
totalGrandeSequence$Fraction<- rep("Grande", each = nrow(totalGrandeSequence))
## Fraction
colnames(totalPetiteSequence)[2]  <- "value"
totalPetiteSequence$Sum <- rep(0, each = nrow(totalPetiteSequence))
totalPetiteSequence$Count <- rep(0, each = nrow(totalPetiteSequence))
for (i in totalPetiteSequence$Division) { totalPetiteSequence$Count[which(totalPetiteSequence$Division == i)] <- sum(data_seq_tax  %>% filter(Division == i) %>% select(TotalPetite))}
totalPetiteSequence$Sum <- sum(totalPetiteSequence$Count)
colnames(totalGrandeSequence)[2]  <- "value"
totalGrandeSequence$Sum <- rep(0, each = nrow(totalGrandeSequence))
totalGrandeSequence$Count <- rep(0, each = nrow(totalGrandeSequence))
for (i in totalGrandeSequence$Division) { totalGrandeSequence$Count[which(totalGrandeSequence$Division == i)] <- sum(data_seq_tax  %>% filter(Division == i) %>% select(TotalGrande))}
totalGrandeSequence$Sum <- sum(totalGrandeSequence$Count)
totalFractionSequence <- rbind(totalPetiteSequence,totalGrandeSequence)
totalFractionSequence$percent <- paste("(",round(totalFractionSequence$Sum*100/(colSums(statRarefy %>% select(`apRarefy-Sequence`))/2),1)," %)",sep ="")

#Figure
ky <- ggplot(totalFractionSequence, mapping = aes(y= value, x = Fraction, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
  scale_fill_manual(values = palette) + 
  geom_label(aes(y = 108,label = paste(Sum,"séquences",percent,sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
ky <- ky + labs(x="Fractions",y="Séquences (%)") + theme(legend.position = "none")
print(ky) 
    # Date -------------------------------------------------------------------
## Avril
total04Sequence <- data_seq_tax %>% select(ASV_Id,Total04,Division)
row.names(total04Sequence)<-total04Sequence$ASV_Id ; total04Sequence <- total04Sequence %>% select(-ASV_Id)
total04Sequence <- total04Sequence %>% group_by(Division) %>% summarise_all(sum)
total04Sequence$Total04 <- total04Sequence$Total04*100/sum(total04Sequence$Total04)
total04Sequence <- total04Sequence %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(Total04) - 0.5*Total04)
total04Sequence$label <- paste(round(total04Sequence$Total04,1), "%", sep = "")
for (i in rownames(total04Sequence)) {
  if (total04Sequence[i,"label"] == "0%") { total04Sequence[i,"label"] <- NA}}
for (i in rownames(total04Sequence)) {
  if (is.na(total04Sequence[i,"label"]) == FALSE) { total04Sequence[i,"label"] <- paste(total04Sequence[i,"Division"]," : ",total04Sequence[i,"label"], sep = "")}}
total04Sequence$Date<- rep("04", each = nrow(total04Sequence))
## Juin
total06Sequence <- data_seq_tax %>% select(ASV_Id,Total06,Division)
row.names(total06Sequence)<-total06Sequence$ASV_Id ; total06Sequence <- total06Sequence %>% select(-ASV_Id)
total06Sequence <- total06Sequence %>% group_by(Division) %>% summarise_all(sum)
total06Sequence$Total06 <- total06Sequence$Total06*100/sum(total06Sequence$Total06)
total06Sequence <- total06Sequence %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(Total06) - 0.5*Total06)
total06Sequence$label <- paste(round(total06Sequence$Total06,1), "%", sep = "")
for (i in rownames(total06Sequence)) {
  if (total06Sequence[i,"label"] == "0%") { total06Sequence[i,"label"] <- NA}}
for (i in rownames(total06Sequence)) {
  if (is.na(total06Sequence[i,"label"]) == FALSE) { total06Sequence[i,"label"] <- paste(total06Sequence[i,"Division"]," : ",total06Sequence[i,"label"], sep = "")}}
total06Sequence$Date<- rep("06", each = nrow(total06Sequence))
## Septembre
total09Sequence <- data_seq_tax %>% select(ASV_Id,Total09,Division)
row.names(total09Sequence)<-total09Sequence$ASV_Id ; total09Sequence <- total09Sequence %>% select(-ASV_Id)
total09Sequence <- total09Sequence %>% group_by(Division) %>% summarise_all(sum)
total09Sequence$Total09 <- total09Sequence$Total09*100/sum(total09Sequence$Total09)
total09Sequence <- total09Sequence %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(Total09) - 0.5*Total09)
total09Sequence$label <- paste(round(total09Sequence$Total09,1), "%", sep = "")
for (i in rownames(total09Sequence)) {
  if (total09Sequence[i,"label"] == "0%") { total09Sequence[i,"label"] <- NA}}
for (i in rownames(total09Sequence)) {
  if (is.na(total09Sequence[i,"label"]) == FALSE) { total09Sequence[i,"label"] <- paste(total09Sequence[i,"Division"]," : ",total09Sequence[i,"label"], sep = "")}}
total09Sequence$Date<- rep("09", each = nrow(total09Sequence))
## Octobre
total11Sequence <- data_seq_tax %>% select(ASV_Id,Total11,Division)
row.names(total11Sequence)<-total11Sequence$ASV_Id ; total11Sequence <- total11Sequence %>% select(-ASV_Id)
total11Sequence <- total11Sequence %>% group_by(Division) %>% summarise_all(sum)
total11Sequence$Total11 <- total11Sequence$Total11*100/sum(total11Sequence$Total11)
total11Sequence <- total11Sequence %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(Total11) - 0.5*Total11)
total11Sequence$label <- paste(round(total11Sequence$Total11,1), "%", sep = "")
for (i in rownames(total11Sequence)) {
  if (total11Sequence[i,"label"] == "0%") { total11Sequence[i,"label"] <- NA}}
for (i in rownames(total11Sequence)) {
  if (is.na(total11Sequence[i,"label"]) == FALSE) { total11Sequence[i,"label"] <- paste(total11Sequence[i,"Division"]," : ",total11Sequence[i,"label"], sep = "")}}
total11Sequence$Date<- rep("11", each = nrow(total11Sequence))
## Dates
colnames(total04Sequence)[2]  <- "value"
total04Sequence$Sum <- rep(0, each = nrow(total04Sequence))
total04Sequence$Count <- rep(0, each = nrow(total04Sequence))
for (i in total04Sequence$Division) { total04Sequence$Count[which(total04Sequence$Division == i)] <- sum(data_seq_tax  %>% filter(Division == i) %>% select(Total04))}
total04Sequence$Sum <- sum(total04Sequence$Count)
colnames(total06Sequence)[2]  <- "value"
total06Sequence$Sum <- rep(0, each = nrow(total06Sequence))
total06Sequence$Count <- rep(0, each = nrow(total06Sequence))
for (i in total06Sequence$Division) { total06Sequence$Count[which(total06Sequence$Division == i)] <- sum(data_seq_tax  %>% filter(Division == i) %>% select(Total06))}
total06Sequence$Sum <- sum(total06Sequence$Count)
colnames(total09Sequence)[2]  <- "value"
total09Sequence$Sum <- rep(0, each = nrow(total09Sequence))
total09Sequence$Count <- rep(0, each = nrow(total09Sequence))
for (i in total09Sequence$Division) { total09Sequence$Count[which(total09Sequence$Division == i)] <- sum(data_seq_tax  %>% filter(Division == i) %>% select(Total09))}
total09Sequence$Sum <- sum(total09Sequence$Count)
colnames(total11Sequence)[2]  <- "value"
total11Sequence$Sum <- rep(0, each = nrow(total11Sequence))
total11Sequence$Count <- rep(0, each = nrow(total11Sequence))
for (i in total11Sequence$Division) { total11Sequence$Count[which(total11Sequence$Division == i)] <- sum(data_seq_tax  %>% filter(Division == i) %>% select(Total11))}
total11Sequence$Sum <- sum(total11Sequence$Count)
totalDateSequence <- rbind(total04Sequence,total06Sequence,total09Sequence,total11Sequence)
totalDateSequence$percent <- paste("(",round(totalDateSequence$Sum*100/(colSums(statRarefy %>% select(`apRarefy-Sequence`))/4),1)," %)",sep ="")

#Figure
ly <- ggplot(totalDateSequence, mapping = aes(y= value, x = Date, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
  scale_fill_manual(values = palette) + 
  geom_label(aes(y = 108,label = paste(Sum,"séquences",percent,sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
ly <- ly + labs(x="Dates",y="Séquences (%)") + theme(legend.position = "none")
print(ly) 
    # Coplot -------------------------------------------------------------------
svglite("Hist-Taxomy/Sequence-Total.svg",width = 12.00,height = 7.00)
b_plot <- plot_grid(iy,jy,ky,ly, legendSequence, ncol = 5, nrow = 1, rel_widths = c(3,3,3,5,2.5),rel_heights = c(3))
print(b_plot)
dev.off()
  # Séquence ONLY---------------------------------------------------------------------
    # 100% - 0% ---------------------------------------------------------------
      # Cycle -------------------------------------------------------------------
## Jour
onlyJourSequence <- data_seq_tax %>% filter(Cycle == "Jour") %>% select(ASV_Id,TotalJour,Division)
row.names(onlyJourSequence)<-onlyJourSequence$ASV_Id ; onlyJourSequence <- onlyJourSequence %>% select(-ASV_Id)
onlyJourSequence <- onlyJourSequence %>% group_by(Division) %>% summarise_all(sum)
onlyJourSequence$TotalJour <- onlyJourSequence$TotalJour*100/sum(onlyJourSequence$TotalJour)
onlyJourSequence <- onlyJourSequence %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalJour) - 0.5*TotalJour)
onlyJourSequence$label <- paste(round(onlyJourSequence$TotalJour,1), "%", sep = "")
for (i in rownames(onlyJourSequence)) {
  if (onlyJourSequence[i,"label"] == "0%") { onlyJourSequence[i,"label"] <- NA}}
for (i in rownames(onlyJourSequence)) {
  if (is.na(onlyJourSequence[i,"label"]) == FALSE) { onlyJourSequence[i,"label"] <- paste(onlyJourSequence[i,"Division"]," : ",onlyJourSequence[i,"label"], sep = "")}}
onlyJourSequence$Cycle<- rep("Jour", each = nrow(onlyJourSequence))
## Nuit
onlyNuitSequence <- data_seq_tax %>% filter(Cycle == "Nuit") %>% select(ASV_Id,TotalNuit,Division)
row.names(onlyNuitSequence)<-onlyNuitSequence$ASV_Id ; onlyNuitSequence <- onlyNuitSequence %>% select(-ASV_Id)
onlyNuitSequence <- onlyNuitSequence %>% group_by(Division) %>% summarise_all(sum)
onlyNuitSequence$TotalNuit <- onlyNuitSequence$TotalNuit*100/sum(onlyNuitSequence$TotalNuit)
onlyNuitSequence <- onlyNuitSequence %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalNuit) - 0.5*TotalNuit)
onlyNuitSequence$label <- paste(round(onlyNuitSequence$TotalNuit,1), "%", sep = "")
for (i in rownames(onlyNuitSequence)) {
  if (onlyNuitSequence[i,"label"] == "0%") { onlyNuitSequence[i,"label"] <- NA}}
for (i in rownames(onlyNuitSequence)) {
  if (is.na(onlyNuitSequence[i,"label"]) == FALSE) { onlyNuitSequence[i,"label"] <- paste(onlyNuitSequence[i,"Division"]," : ",onlyNuitSequence[i,"label"], sep = "")}}
onlyNuitSequence$Cycle<- rep("Nuit", each = nrow(onlyNuitSequence))
## Cycle
colnames(onlyJourSequence)[2]  <- "value"
onlyJourSequence$Sum <- rep(0, each = nrow(onlyJourSequence))
onlyJourSequence$Count <- rep(0, each = nrow(onlyJourSequence))
for (i in onlyJourSequence$Division) { onlyJourSequence$Count[which(onlyJourSequence$Division == i)] <- sum(data_seq_tax %>% filter(Cycle == "Jour") %>% filter(Division == i) %>% select(TotalJour))}
onlyJourSequence$Sum <- sum(onlyJourSequence$Count)
colnames(onlyNuitSequence)[2]  <- "value"
onlyNuitSequence$Sum <- rep(0, each = nrow(onlyNuitSequence))
onlyNuitSequence$Count <- rep(0, each = nrow(onlyNuitSequence))
for (i in onlyNuitSequence$Division) { onlyNuitSequence$Count[which(onlyNuitSequence$Division == i)] <- sum(data_seq_tax %>% filter(Cycle == "Nuit") %>% filter(Division == i) %>% select(TotalNuit))}
onlyNuitSequence$Sum <- sum(onlyNuitSequence$Count)
onlyCycleSequence <- rbind(onlyJourSequence,onlyNuitSequence)
#Figure
iz <- ggplot(onlyCycleSequence, mapping = aes(y= value, x = Cycle, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
  geom_label(aes(y = 106,label = paste(Sum,"séquences",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
  scale_fill_manual(values = palette)
#legendSequence <- get_legend(iz)
iz <- iz + labs(x="Cycles",y="Séquences (%)") #+ theme(legend.position = "none")
print(iz)

      # Zone -------------------------------------------------------------------
## Oxique
onlyOxiqueSequence <- data_seq_tax %>% filter(Zone == "Oxique") %>% select(ASV_Id,TotalOxique,Division)
row.names(onlyOxiqueSequence)<-onlyOxiqueSequence$ASV_Id ; onlyOxiqueSequence <- onlyOxiqueSequence %>% select(-ASV_Id)
onlyOxiqueSequence <- onlyOxiqueSequence %>% group_by(Division) %>% summarise_all(sum)
onlyOxiqueSequence$TotalOxique <- onlyOxiqueSequence$TotalOxique*100/sum(onlyOxiqueSequence$TotalOxique)
onlyOxiqueSequence <- onlyOxiqueSequence %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalOxique) - 0.5*TotalOxique)
onlyOxiqueSequence$label <- paste(round(onlyOxiqueSequence$TotalOxique,1), "%", sep = "")
for (i in rownames(onlyOxiqueSequence)) {
  if (onlyOxiqueSequence[i,"label"] == "0%") { onlyOxiqueSequence[i,"label"] <- NA}}
for (i in rownames(onlyOxiqueSequence)) {
  if (is.na(onlyOxiqueSequence[i,"label"]) == FALSE) { onlyOxiqueSequence[i,"label"] <- paste(onlyOxiqueSequence[i,"Division"]," : ",onlyOxiqueSequence[i,"label"], sep = "")}}
onlyOxiqueSequence$Zone<- rep("Oxique", each = nrow(onlyOxiqueSequence))
## Anoxique
onlyAnoxiqueSequence <- data_seq_tax %>% filter(Zone == "Anoxique") %>% select(ASV_Id,TotalAnoxique,Division)
row.names(onlyAnoxiqueSequence)<-onlyAnoxiqueSequence$ASV_Id ; onlyAnoxiqueSequence <- onlyAnoxiqueSequence %>% select(-ASV_Id)
onlyAnoxiqueSequence <- onlyAnoxiqueSequence %>% group_by(Division) %>% summarise_all(sum)
onlyAnoxiqueSequence$TotalAnoxique <- onlyAnoxiqueSequence$TotalAnoxique*100/sum(onlyAnoxiqueSequence$TotalAnoxique)
onlyAnoxiqueSequence <- onlyAnoxiqueSequence %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalAnoxique) - 0.5*TotalAnoxique)
onlyAnoxiqueSequence$label <- paste(round(onlyAnoxiqueSequence$TotalAnoxique,1), "%", sep = "")
for (i in rownames(onlyAnoxiqueSequence)) {
  if (onlyAnoxiqueSequence[i,"label"] == "0%") { onlyAnoxiqueSequence[i,"label"] <- NA}}
for (i in rownames(onlyAnoxiqueSequence)) {
  if (is.na(onlyAnoxiqueSequence[i,"label"]) == FALSE) { onlyAnoxiqueSequence[i,"label"] <- paste(onlyAnoxiqueSequence[i,"Division"]," : ",onlyAnoxiqueSequence[i,"label"], sep = "")}}
onlyAnoxiqueSequence$Zone<- rep("Anoxique", each = nrow(onlyAnoxiqueSequence))
## Zone
colnames(onlyOxiqueSequence)[2]  <- "value"
onlyOxiqueSequence$Sum <- rep(0, each = nrow(onlyOxiqueSequence))
onlyOxiqueSequence$Count <- rep(0, each = nrow(onlyOxiqueSequence))
for (i in onlyOxiqueSequence$Division) { onlyOxiqueSequence$Count[which(onlyOxiqueSequence$Division == i)] <- sum(data_seq_tax %>% filter(Zone == "Oxique") %>% filter(Division == i) %>% select(TotalOxique))}
onlyOxiqueSequence$Sum <- sum(onlyOxiqueSequence$Count)
colnames(onlyAnoxiqueSequence)[2]  <- "value"
onlyAnoxiqueSequence$Sum <- rep(0, each = nrow(onlyAnoxiqueSequence))
onlyAnoxiqueSequence$Count <- rep(0, each = nrow(onlyAnoxiqueSequence))
for (i in onlyAnoxiqueSequence$Division) { onlyAnoxiqueSequence$Count[which(onlyAnoxiqueSequence$Division == i)] <- sum(data_seq_tax %>% filter(Zone == "Anoxique") %>% filter(Division == i) %>% select(TotalAnoxique))}
onlyAnoxiqueSequence$Sum <- sum(onlyAnoxiqueSequence$Count)
onlyZoneSequence <- rbind(onlyOxiqueSequence,onlyAnoxiqueSequence)
#Figure
jz <- ggplot(onlyZoneSequence, mapping = aes(y= value, x = Zone, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
  geom_label(aes(y = 106,label = paste(Sum,"séquences",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
  scale_fill_manual(values = palette)
jz <- jz + labs(x="Zones",y="Séquences (%)") #+ theme(legend.position = "none") 
print(jz)    
      # Fraction -------------------------------------------------------------------
## Petite
onlyPetiteSequence <- data_seq_tax %>% filter(Fraction == "Petite") %>% select(ASV_Id,TotalPetite,Division)
row.names(onlyPetiteSequence)<-onlyPetiteSequence$ASV_Id ; onlyPetiteSequence <- onlyPetiteSequence %>% select(-ASV_Id)
onlyPetiteSequence <- onlyPetiteSequence %>% group_by(Division) %>% summarise_all(sum)
onlyPetiteSequence$TotalPetite <- onlyPetiteSequence$TotalPetite*100/sum(onlyPetiteSequence$TotalPetite)
onlyPetiteSequence <- onlyPetiteSequence %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalPetite) - 0.5*TotalPetite)
onlyPetiteSequence$label <- paste(round(onlyPetiteSequence$TotalPetite,1), "%", sep = "")
for (i in rownames(onlyPetiteSequence)) {
  if (onlyPetiteSequence[i,"label"] == "0%") { onlyPetiteSequence[i,"label"] <- NA}}
for (i in rownames(onlyPetiteSequence)) {
  if (is.na(onlyPetiteSequence[i,"label"]) == FALSE) { onlyPetiteSequence[i,"label"] <- paste(onlyPetiteSequence[i,"Division"]," : ",onlyPetiteSequence[i,"label"], sep = "")}}
onlyPetiteSequence$Fraction<- rep("Petite", each = nrow(onlyPetiteSequence))
## Grande
onlyGrandeSequence <- data_seq_tax %>% filter(Fraction == "Grande") %>% select(ASV_Id,TotalGrande,Division)
row.names(onlyGrandeSequence)<-onlyGrandeSequence$ASV_Id ; onlyGrandeSequence <- onlyGrandeSequence %>% select(-ASV_Id)
onlyGrandeSequence <- onlyGrandeSequence %>% group_by(Division) %>% summarise_all(sum)
onlyGrandeSequence$TotalGrande <- onlyGrandeSequence$TotalGrande*100/sum(onlyGrandeSequence$TotalGrande)
onlyGrandeSequence <- onlyGrandeSequence %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalGrande) - 0.5*TotalGrande)
onlyGrandeSequence$label <- paste(round(onlyGrandeSequence$TotalGrande,1), "%", sep = "")
for (i in rownames(onlyGrandeSequence)) {
  if (onlyGrandeSequence[i,"label"] == "0%") { onlyGrandeSequence[i,"label"] <- NA}}
for (i in rownames(onlyGrandeSequence)) {
  if (is.na(onlyGrandeSequence[i,"label"]) == FALSE) { onlyGrandeSequence[i,"label"] <- paste(onlyGrandeSequence[i,"Division"]," : ",onlyGrandeSequence[i,"label"], sep = "")}}
onlyGrandeSequence$Fraction<- rep("Grande", each = nrow(onlyGrandeSequence))
## Fraction
colnames(onlyPetiteSequence)[2]  <- "value"
onlyPetiteSequence$Sum <- rep(0, each = nrow(onlyPetiteSequence))
onlyPetiteSequence$Count <- rep(0, each = nrow(onlyPetiteSequence))
for (i in onlyPetiteSequence$Division) { onlyPetiteSequence$Count[which(onlyPetiteSequence$Division == i)] <- sum(data_seq_tax %>% filter(Fraction == "Petite") %>% filter(Division == i) %>% select(TotalPetite))}
onlyPetiteSequence$Sum <- sum(onlyPetiteSequence$Count)
colnames(onlyGrandeSequence)[2]  <- "value"
onlyGrandeSequence$Sum <- rep(0, each = nrow(onlyGrandeSequence))
onlyGrandeSequence$Count <- rep(0, each = nrow(onlyGrandeSequence))
for (i in onlyGrandeSequence$Division) { onlyGrandeSequence$Count[which(onlyGrandeSequence$Division == i)] <- sum(data_seq_tax %>% filter(Fraction == "Grande") %>% filter(Division == i) %>% select(TotalGrande))}
onlyGrandeSequence$Sum <- sum(onlyGrandeSequence$Count)
onlyFractionSequence <- rbind(onlyPetiteSequence,onlyGrandeSequence)
#Figure
kz <- ggplot(onlyFractionSequence, mapping = aes(y= value, x = Fraction, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
  geom_label(aes(y = 106,label = paste(Sum,"séquences",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
  scale_fill_manual(values = palette)
kz <- kz + labs(x="Fractions",y="Séquences (%)") #+ theme(legend.position = "none")
print(kz) 
      # Coplot -------------------------------------------------------------------
svglite("Hist-Taxomy/Sequence-only-100.svg",width = 12.00,height = 6.00)
b_plot <- plot_grid(iz,jz,kz,  ncol = 3, nrow = 1, rel_widths = c(3,3,3),rel_heights = c(3))
print(b_plot)
dev.off()  

      # Coplot X AFC -------------------------------------------------------------------
f <- ggplot() + theme_void()
print(f)
#AFC = séquence
svglite("Biplot/Sequence-onlyX-AFCseq-100.svg",width = 12.00,height = 14.00)
b_plot <- plot_grid(aCycle,iz,bZone,jz,cFraction,kz, ncol = 2, nrow = 3, rel_widths = c(3,2),rel_heights = c(3,3,3))
print(b_plot)
dev.off()  
#AFC = ASV
svglite("Biplot/Sequence-onlyX-AFCasv-100.svg",width = 12.00,height = 14.00)
c_plot <- plot_grid(iCycle,iz,jZone,jz,kFraction,kz, ncol = 2, nrow = 3, rel_widths = c(3,2),rel_heights = c(3,3,3))
print(c_plot)
dev.off()


    # 90% - 10% ---------------------------------------------------------------
      # Cycle -------------------------------------------------------------------
## Jour
onlyJourSequence90 <- data_seq_tax %>% filter(Cycle90 == "Jour") %>% select(ASV_Id,TotalJour,Division)
row.names(onlyJourSequence90)<-onlyJourSequence90$ASV_Id ; onlyJourSequence90 <- onlyJourSequence90 %>% select(-ASV_Id)
onlyJourSequence90 <- onlyJourSequence90 %>% group_by(Division) %>% summarise_all(sum)
onlyJourSequence90$TotalJour <- onlyJourSequence90$TotalJour*100/sum(onlyJourSequence90$TotalJour)
onlyJourSequence90 <- onlyJourSequence90 %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalJour) - 0.5*TotalJour)
onlyJourSequence90$label <- paste(round(onlyJourSequence90$TotalJour,1), "%", sep = "")
for (i in rownames(onlyJourSequence90)) {
  if (onlyJourSequence90[i,"label"] == "0%") { onlyJourSequence90[i,"label"] <- NA}}
for (i in rownames(onlyJourSequence90)) {
  if (is.na(onlyJourSequence90[i,"label"]) == FALSE) { onlyJourSequence90[i,"label"] <- paste(onlyJourSequence90[i,"Division"]," : ",onlyJourSequence90[i,"label"], sep = "")}}
onlyJourSequence90$Cycle90<- rep("Jour", each = nrow(onlyJourSequence90))
## Nuit
onlyNuitSequence90 <- data_seq_tax %>% filter(Cycle90 == "Nuit") %>% select(ASV_Id,TotalNuit,Division)
row.names(onlyNuitSequence90)<-onlyNuitSequence90$ASV_Id ; onlyNuitSequence90 <- onlyNuitSequence90 %>% select(-ASV_Id)
onlyNuitSequence90 <- onlyNuitSequence90 %>% group_by(Division) %>% summarise_all(sum)
onlyNuitSequence90$TotalNuit <- onlyNuitSequence90$TotalNuit*100/sum(onlyNuitSequence90$TotalNuit)
onlyNuitSequence90 <- onlyNuitSequence90 %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalNuit) - 0.5*TotalNuit)
onlyNuitSequence90$label <- paste(round(onlyNuitSequence90$TotalNuit,1), "%", sep = "")
for (i in rownames(onlyNuitSequence90)) {
  if (onlyNuitSequence90[i,"label"] == "0%") { onlyNuitSequence90[i,"label"] <- NA}}
for (i in rownames(onlyNuitSequence90)) {
  if (is.na(onlyNuitSequence90[i,"label"]) == FALSE) { onlyNuitSequence90[i,"label"] <- paste(onlyNuitSequence90[i,"Division"]," : ",onlyNuitSequence90[i,"label"], sep = "")}}
onlyNuitSequence90$Cycle90<- rep("Nuit", each = nrow(onlyNuitSequence90))
## Cycle
colnames(onlyJourSequence90)[2]  <- "value"
onlyJourSequence90$Sum <- rep(0, each = nrow(onlyJourSequence90))
onlyJourSequence90$Count <- rep(0, each = nrow(onlyJourSequence90))
for (i in onlyJourSequence90$Division) { onlyJourSequence90$Count[which(onlyJourSequence90$Division == i)] <- sum(data_seq_tax %>% filter(Cycle90 == "Jour") %>% filter(Division == i) %>% select(TotalJour))}
onlyJourSequence90$Sum <- sum(onlyJourSequence90$Count)
colnames(onlyNuitSequence90)[2]  <- "value"
onlyNuitSequence90$Sum <- rep(0, each = nrow(onlyNuitSequence90))
onlyNuitSequence90$Count <- rep(0, each = nrow(onlyNuitSequence90))
for (i in onlyNuitSequence90$Division) { onlyNuitSequence90$Count[which(onlyNuitSequence90$Division == i)] <- sum(data_seq_tax %>% filter(Cycle90 == "Nuit") %>% filter(Division == i) %>% select(TotalNuit))}
onlyNuitSequence90$Sum <- sum(onlyNuitSequence90$Count)
onlyCycleSequence90 <- rbind(onlyJourSequence90,onlyNuitSequence90)
#Figure
iz90 <- ggplot(onlyCycleSequence90, mapping = aes(y= value, x = Cycle90, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
  geom_label(aes(y = 106,label = paste(Sum,"séquences",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
  scale_fill_manual(values = palette)
#legendSequence <- get_legend(iz90)
iz90 <- iz90 + labs(x="Cycles",y="Séquences (%)") #+ theme(legend.position = "none")
print(iz90)

      # Zone -------------------------------------------------------------------
## Oxique
onlyOxiqueSequence90 <- data_seq_tax %>% filter(Zone90 == "Oxique") %>% select(ASV_Id,TotalOxique,Division)
row.names(onlyOxiqueSequence90)<-onlyOxiqueSequence90$ASV_Id ; onlyOxiqueSequence90 <- onlyOxiqueSequence90 %>% select(-ASV_Id)
onlyOxiqueSequence90 <- onlyOxiqueSequence90 %>% group_by(Division) %>% summarise_all(sum)
onlyOxiqueSequence90$TotalOxique <- onlyOxiqueSequence90$TotalOxique*100/sum(onlyOxiqueSequence90$TotalOxique)
onlyOxiqueSequence90 <- onlyOxiqueSequence90 %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalOxique) - 0.5*TotalOxique)
onlyOxiqueSequence90$label <- paste(round(onlyOxiqueSequence90$TotalOxique,1), "%", sep = "")
for (i in rownames(onlyOxiqueSequence90)) {
  if (onlyOxiqueSequence90[i,"label"] == "0%") { onlyOxiqueSequence90[i,"label"] <- NA}}
for (i in rownames(onlyOxiqueSequence90)) {
  if (is.na(onlyOxiqueSequence90[i,"label"]) == FALSE) { onlyOxiqueSequence90[i,"label"] <- paste(onlyOxiqueSequence90[i,"Division"]," : ",onlyOxiqueSequence90[i,"label"], sep = "")}}
onlyOxiqueSequence90$Zone90<- rep("Oxique", each = nrow(onlyOxiqueSequence90))
## Anoxique
onlyAnoxiqueSequence90 <- data_seq_tax %>% filter(Zone90 == "Anoxique") %>% select(ASV_Id,TotalAnoxique,Division)
row.names(onlyAnoxiqueSequence90)<-onlyAnoxiqueSequence90$ASV_Id ; onlyAnoxiqueSequence90 <- onlyAnoxiqueSequence90 %>% select(-ASV_Id)
onlyAnoxiqueSequence90 <- onlyAnoxiqueSequence90 %>% group_by(Division) %>% summarise_all(sum)
onlyAnoxiqueSequence90$TotalAnoxique <- onlyAnoxiqueSequence90$TotalAnoxique*100/sum(onlyAnoxiqueSequence90$TotalAnoxique)
onlyAnoxiqueSequence90 <- onlyAnoxiqueSequence90 %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalAnoxique) - 0.5*TotalAnoxique)
onlyAnoxiqueSequence90$label <- paste(round(onlyAnoxiqueSequence90$TotalAnoxique,1), "%", sep = "")
for (i in rownames(onlyAnoxiqueSequence90)) {
  if (onlyAnoxiqueSequence90[i,"label"] == "0%") { onlyAnoxiqueSequence90[i,"label"] <- NA}}
for (i in rownames(onlyAnoxiqueSequence90)) {
  if (is.na(onlyAnoxiqueSequence90[i,"label"]) == FALSE) { onlyAnoxiqueSequence90[i,"label"] <- paste(onlyAnoxiqueSequence90[i,"Division"]," : ",onlyAnoxiqueSequence90[i,"label"], sep = "")}}
onlyAnoxiqueSequence90$Zone90<- rep("Anoxique", each = nrow(onlyAnoxiqueSequence90))
## Zone
colnames(onlyOxiqueSequence90)[2]  <- "value"
onlyOxiqueSequence90$Sum <- rep(0, each = nrow(onlyOxiqueSequence90))
onlyOxiqueSequence90$Count <- rep(0, each = nrow(onlyOxiqueSequence90))
for (i in onlyOxiqueSequence90$Division) { onlyOxiqueSequence90$Count[which(onlyOxiqueSequence90$Division == i)] <- sum(data_seq_tax %>% filter(Zone90 == "Oxique") %>% filter(Division == i) %>% select(TotalOxique))}
onlyOxiqueSequence90$Sum <- sum(onlyOxiqueSequence90$Count)
colnames(onlyAnoxiqueSequence90)[2]  <- "value"
onlyAnoxiqueSequence90$Sum <- rep(0, each = nrow(onlyAnoxiqueSequence90))
onlyAnoxiqueSequence90$Count <- rep(0, each = nrow(onlyAnoxiqueSequence90))
for (i in onlyAnoxiqueSequence90$Division) { onlyAnoxiqueSequence90$Count[which(onlyAnoxiqueSequence90$Division == i)] <- sum(data_seq_tax %>% filter(Zone90 == "Anoxique") %>% filter(Division == i) %>% select(TotalAnoxique))}
onlyAnoxiqueSequence90$Sum <- sum(onlyAnoxiqueSequence90$Count)
onlyZoneSequence90 <- rbind(onlyOxiqueSequence90,onlyAnoxiqueSequence90)
#Figure
jz90 <- ggplot(onlyZoneSequence90, mapping = aes(y= value, x = Zone90, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
  geom_label(aes(y = 106,label = paste(Sum,"séquences",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
  scale_fill_manual(values = palette)
jz90 <- jz90 + labs(x="Zones",y="Séquences (%)") #+ theme(legend.position = "none") 
print(jz90)    
      # Fraction -------------------------------------------------------------------
## Petite
onlyPetiteSequence90 <- data_seq_tax %>% filter(Fraction90 == "Petite") %>% select(ASV_Id,TotalPetite,Division)
row.names(onlyPetiteSequence90)<-onlyPetiteSequence90$ASV_Id ; onlyPetiteSequence90 <- onlyPetiteSequence90 %>% select(-ASV_Id)
onlyPetiteSequence90 <- onlyPetiteSequence90 %>% group_by(Division) %>% summarise_all(sum)
onlyPetiteSequence90$TotalPetite <- onlyPetiteSequence90$TotalPetite*100/sum(onlyPetiteSequence90$TotalPetite)
onlyPetiteSequence90 <- onlyPetiteSequence90 %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalPetite) - 0.5*TotalPetite)
onlyPetiteSequence90$label <- paste(round(onlyPetiteSequence90$TotalPetite,1), "%", sep = "")
for (i in rownames(onlyPetiteSequence90)) {
  if (onlyPetiteSequence90[i,"label"] == "0%") { onlyPetiteSequence90[i,"label"] <- NA}}
for (i in rownames(onlyPetiteSequence90)) {
  if (is.na(onlyPetiteSequence90[i,"label"]) == FALSE) { onlyPetiteSequence90[i,"label"] <- paste(onlyPetiteSequence90[i,"Division"]," : ",onlyPetiteSequence90[i,"label"], sep = "")}}
onlyPetiteSequence90$Fraction90<- rep("Petite", each = nrow(onlyPetiteSequence90))
## Grande
onlyGrandeSequence90 <- data_seq_tax %>% filter(Fraction90 == "Grande") %>% select(ASV_Id,TotalGrande,Division)
row.names(onlyGrandeSequence90)<-onlyGrandeSequence90$ASV_Id ; onlyGrandeSequence90 <- onlyGrandeSequence90 %>% select(-ASV_Id)
onlyGrandeSequence90 <- onlyGrandeSequence90 %>% group_by(Division) %>% summarise_all(sum)
onlyGrandeSequence90$TotalGrande <- onlyGrandeSequence90$TotalGrande*100/sum(onlyGrandeSequence90$TotalGrande)
onlyGrandeSequence90 <- onlyGrandeSequence90 %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalGrande) - 0.5*TotalGrande)
onlyGrandeSequence90$label <- paste(round(onlyGrandeSequence90$TotalGrande,1), "%", sep = "")
for (i in rownames(onlyGrandeSequence90)) {
  if (onlyGrandeSequence90[i,"label"] == "0%") { onlyGrandeSequence90[i,"label"] <- NA}}
for (i in rownames(onlyGrandeSequence90)) {
  if (is.na(onlyGrandeSequence90[i,"label"]) == FALSE) { onlyGrandeSequence90[i,"label"] <- paste(onlyGrandeSequence90[i,"Division"]," : ",onlyGrandeSequence90[i,"label"], sep = "")}}
onlyGrandeSequence90$Fraction90<- rep("Grande", each = nrow(onlyGrandeSequence90))
## Fraction
colnames(onlyPetiteSequence90)[2]  <- "value"
onlyPetiteSequence90$Sum <- rep(0, each = nrow(onlyPetiteSequence90))
onlyPetiteSequence90$Count <- rep(0, each = nrow(onlyPetiteSequence90))
for (i in onlyPetiteSequence90$Division) { onlyPetiteSequence90$Count[which(onlyPetiteSequence90$Division == i)] <- sum(data_seq_tax %>% filter(Fraction90 == "Petite") %>% filter(Division == i) %>% select(TotalPetite))}
onlyPetiteSequence90$Sum <- sum(onlyPetiteSequence90$Count)
colnames(onlyGrandeSequence90)[2]  <- "value"
onlyGrandeSequence90$Sum <- rep(0, each = nrow(onlyGrandeSequence90))
onlyGrandeSequence90$Count <- rep(0, each = nrow(onlyGrandeSequence90))
for (i in onlyGrandeSequence90$Division) { onlyGrandeSequence90$Count[which(onlyGrandeSequence90$Division == i)] <- sum(data_seq_tax %>% filter(Fraction90 == "Grande") %>% filter(Division == i) %>% select(TotalGrande))}
onlyGrandeSequence90$Sum <- sum(onlyGrandeSequence90$Count)
onlyFractionSequence90 <- rbind(onlyPetiteSequence90,onlyGrandeSequence90)
#Figure
kz90 <- ggplot(onlyFractionSequence90, mapping = aes(y= value, x = Fraction90, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
  geom_label(aes(y = 106,label = paste(Sum,"séquences",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
  scale_fill_manual(values = palette)
kz90 <- kz90 + labs(x="Fractions",y="Séquences (%)") #+ theme(legend.position = "none")
print(kz90) 
      # Coplot -------------------------------------------------------------------
svglite("Hist-Taxomy/Sequence-only-90.svg",width = 12.00,height = 6.00)
b_plot <- plot_grid(iz90,jz90,kz90,  ncol = 3, nrow = 1, rel_widths = c(3,3,3),rel_heights = c(3))
print(b_plot)
dev.off()  

      # Coplot X AFC -------------------------------------------------------------------
f <- ggplot() + theme_void()
print(f)
#AFC = séquence
svglite("Biplot/Sequence-onlyX-AFCseq-90.svg",width = 12.00,height = 14.00)
b_plot <- plot_grid(aCycle90,iz90,bZone90,jz90,cFraction90,kz90, ncol = 2, nrow = 3, rel_widths = c(3,2),rel_heights = c(3,3,3))
print(b_plot)
dev.off()  
#AFC = ASV
svglite("Biplot/Sequence-onlyX-AFCasv-90.svg",width = 12.00,height = 14.00)
c_plot <- plot_grid(iCycle,iz90,jZone,jz90,kFraction,kz90, ncol = 2, nrow = 3, rel_widths = c(3,2),rel_heights = c(3,3,3))
print(c_plot)
dev.off()


# Polar ------------------------------------------------------
  # Séquence ----------------------------------------------------------------
##Total column
Polar_seq <- data_seq_tax
Polar_seq$TotalAnnée <- 0
for (i in rownames(Polar_seq)) { Polar_seq[i,"TotalAnnée"] <- Polar_seq[i,"Total04"] + Polar_seq[i,"Total06"] + Polar_seq[i,"Total09"] + Polar_seq[i,"Total11"]}
Polar_seq <- Polar_seq %>% select(TotalAnnée,Division)
## Total Figure
Polar_seq <- Polar_seq %>% group_by(Division) %>% summarise_all(sum)
for (i in rownames(Polar_seq)) { Polar_seq[i,"Total"] <- (Polar_seq[i,"TotalAnnée"] * 100) / sum(Polar_seq$TotalAnnée)}
##Label
Polar_seq <- Polar_seq %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(Total) - 0.5*Total)
Polar_seq$label <- paste(round(Polar_seq$Total,1), "%", sep = "")
for (i in rownames(Polar_seq)) {
  if (Polar_seq[i,"label"] == "0%") { Polar_seq[i,"label"] <- NA}}
for (i in rownames(Polar_seq)) {
  if (is.na(Polar_seq[i,"label"]) == FALSE) { Polar_seq[i,"label"] <- paste(Polar_seq[i,"Division"]," : ",Polar_seq[i,"label"], sep = "")}}
##Figure
svglite("Composition/Polar-Total-seq.svg",width = 4.50,height = 4.50)
ax <- ggplot(Polar_seq, mapping = aes(y= Total, x = 2, fill = Division), Rowv = NA, col = colMain, scale = "column") +
  geom_bar(stat="identity", color = "white", width = 1) + coord_polar("y") + 
  geom_label_repel(aes(y = lab.ypos,label = label),color = "white",size = 3,segment.color = "black",show.legend = FALSE, nudge_x = 0.5) +
  scale_y_continuous(limits=c(0,sum(Polar_seq %>% select(Total))))+
  xlim(1,2.5) +
  theme_unique_darkbis() + 
  #facet_wrap( ~ variable , nrow = 4) +
  labs(y = "Séquences",x="") + scale_fill_manual(values = rev(palette))
print(ax)
dev.off()

  # ASV ----------------------------------------------------------------
##Total column
Polar_asv <- data_asv_tax
Polar_asv$TotalAnnée <- 0
for (i in rownames(Polar_asv)) { if (Polar_asv[i,"Total04"] + Polar_asv[i,"Total06"] + Polar_asv[i,"Total09"] + Polar_asv[i,"Total11"] > 0) {Polar_asv[i,"TotalAnnée"] <- 1}}
Polar_asv <- Polar_asv %>% select(TotalAnnée,Division)
## Total Figure
Polar_asv <- Polar_asv %>% group_by(Division) %>% summarise_all(sum)
for (i in rownames(Polar_asv)) { Polar_asv[i,"Total"] <- (Polar_asv[i,"TotalAnnée"] * 100) / sum(Polar_asv$TotalAnnée)}
##Label
Polar_asv <- Polar_asv %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(Total) - 0.5*Total)
Polar_asv$label <- paste(round(Polar_asv$Total,1), "%", sep = "")
for (i in rownames(Polar_asv)) {
  if (Polar_asv[i,"label"] == "0%") { Polar_asv[i,"label"] <- NA}}
for (i in rownames(Polar_asv)) {
  if (is.na(Polar_asv[i,"label"]) == FALSE) { Polar_asv[i,"label"] <- paste(Polar_asv[i,"Division"]," : ",Polar_asv[i,"label"], sep = "")}}
##Figure
svglite("Composition/Polar-Total-asv.svg",width = 4.50,height = 4.50)
bx <- ggplot(Polar_asv, mapping = aes(y= Total, x = 2, fill = Division), Rowv = NA, col = colMain, scale = "column") +
  geom_bar(stat="identity", color = "white", width = 1) + coord_polar("y") + 
  geom_label_repel(aes(y = lab.ypos,label = label),color = "white",size = 3,segment.color = "black",show.legend = FALSE, nudge_x = 0.5) +
  scale_y_continuous(limits=c(0,sum(Polar_asv %>% select(Total))))+
  xlim(1,2.5) +
  theme_unique_darkbis() + 
  #facet_wrap( ~ variable , nrow = 4) +
  labs(y = "ASVs",x="") + scale_fill_manual(values = rev(palette))
print(bx)
dev.off()

# Table ASVs majoritaires -------------------------------------------------
# Table ONLY---------------------------------------------------------------------
  # Cycle-100 -------------------------------------------------------------------
##Jour
onlyJourTable <- data_seq_tax %>% filter(Cycle == "Jour") %>% select(ASV_Id,TotalJour,Division,Cycle)
onlyJourTable <- onlyJourTable %>% filter(TotalJour > 0.01*sum(onlyJourTable$TotalJour))
colnames(onlyJourTable)[2]  <- "value"
##Nuit
onlyNuitTable <- data_seq_tax %>% filter(Cycle == "Nuit") %>% select(ASV_Id,TotalNuit,Division,Cycle)
onlyNuitTable <- onlyNuitTable %>% filter(TotalNuit > 0.01*sum(onlyNuitTable$TotalNuit))
colnames(onlyNuitTable)[2]  <- "value"
##Bind
onlyCycleTable <- rbind(onlyJourTable,onlyNuitTable)
onlyCycleTable <- merge(onlyCycleTable,tax_tablemix, by = "ASV_Id")
##Table
write.table(onlyCycleTable, file = "Table/Only_Cycles.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(onlyCycleASV, file = "Table/Only_Cycles_ASV.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(onlyCycleSequence, file = "Table/Only_Cycles_Sequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  # Zone-100 -------------------------------------------------------------------
##Oxique
onlyOxiqueTable <- data_seq_tax %>% filter(Zone == "Oxique") %>% select(ASV_Id,TotalOxique,Division,Zone)
onlyOxiqueTable <- onlyOxiqueTable %>% filter(TotalOxique > 0.01*sum(onlyOxiqueTable$TotalOxique))
colnames(onlyOxiqueTable)[2]  <- "value"
##Anoxique
onlyAnoxiqueTable <- data_seq_tax %>% filter(Zone == "Anoxique") %>% select(ASV_Id,TotalAnoxique,Division,Zone)
onlyAnoxiqueTable <- onlyAnoxiqueTable %>% filter(TotalAnoxique > 0.01*sum(onlyAnoxiqueTable$TotalAnoxique))
colnames(onlyAnoxiqueTable)[2]  <- "value"
##bind
onlyZoneTable <- rbind(onlyOxiqueTable,onlyAnoxiqueTable)
onlyZoneTable <- merge(onlyZoneTable,tax_tablemix, by = "ASV_Id")
##Table
write.table(onlyZoneTable, file = "Table/Only_Zone.txt", sep = "\t", col.names = TRUE, row.names = FALSE,quote = FALSE)
write.table(onlyZoneASV, file = "Table/Only_Zone_ASV.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(onlyZoneSequence, file = "Table/Only_Zone_Sequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  # Fraction-100 -------------------------------------------------------------------
##Petite
onlyPetiteTable <- data_seq_tax %>% filter(Fraction == "Petite") %>% select(ASV_Id,TotalPetite,Division,Fraction)
onlyPetiteTable <- onlyPetiteTable %>% filter(TotalPetite > 0.01*sum(onlyPetiteTable$TotalPetite))
colnames(onlyPetiteTable)[2]  <- "value"
##Grande
onlyGrandeTable <- data_seq_tax %>% filter(Fraction == "Grande") %>% select(ASV_Id,TotalGrande,Division,Fraction)
onlyGrandeTable <- onlyGrandeTable %>% filter(TotalGrande > 0.01*sum(onlyGrandeTable$TotalGrande))
colnames(onlyGrandeTable)[2]  <- "value"
##bind
onlyFractionTable <- rbind(onlyPetiteTable,onlyGrandeTable)
onlyFractionTable <- merge(onlyFractionTable,tax_tablemix, by = "ASV_Id")
##Table
write.table(onlyFractionTable, file = "Table/Only_Fraction.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(onlyFractionASV, file = "Table/Only_Fraction_ASV.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(onlyFractionSequence, file = "Table/Only_Fraction_Sequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  # Cycle-90 -------------------------------------------------------------------
##Jour
onlyJourTable90 <- data_seq_tax %>% filter(Cycle90 == "Jour") %>% select(ASV_Id,TotalJour,Division,Cycle90)
onlyJourTable90 <- onlyJourTable90 %>% filter(TotalJour > 0.01*sum(onlyJourTable90$TotalJour))
colnames(onlyJourTable90)[2]  <- "value"
##Nuit
onlyNuitTable90 <- data_seq_tax %>% filter(Cycle90 == "Nuit") %>% select(ASV_Id,TotalNuit,Division,Cycle90)
onlyNuitTable90 <- onlyNuitTable90 %>% filter(TotalNuit > 0.01*sum(onlyNuitTable90$TotalNuit))
colnames(onlyNuitTable90)[2]  <- "value"
##Bind
onlyCycleTable90 <- rbind(onlyJourTable90,onlyNuitTable90)
onlyCycleTable90 <- merge(onlyCycleTable90,tax_tablemix, by = "ASV_Id")
##Table
write.table(onlyCycleTable90, file = "Table/Only_Cycles90.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(onlyCycleSequence90, file = "Table/Only_Cycles_Sequence90.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  # Zone-90 -------------------------------------------------------------------
##Oxique
onlyOxiqueTable90 <- data_seq_tax %>% filter(Zone90 == "Oxique") %>% select(ASV_Id,TotalOxique,Division,Zone90)
onlyOxiqueTable90 <- onlyOxiqueTable90 %>% filter(TotalOxique > 0.01*sum(onlyOxiqueTable90$TotalOxique))
colnames(onlyOxiqueTable90)[2]  <- "value"
##Anoxique
onlyAnoxiqueTable90 <- data_seq_tax %>% filter(Zone90 == "Anoxique") %>% select(ASV_Id,TotalAnoxique,Division,Zone90)
onlyAnoxiqueTable90 <- onlyAnoxiqueTable90 %>% filter(TotalAnoxique > 0.01*sum(onlyAnoxiqueTable90$TotalAnoxique))
colnames(onlyAnoxiqueTable90)[2]  <- "value"
##bind
onlyZoneTable90 <- rbind(onlyOxiqueTable90,onlyAnoxiqueTable90)
onlyZoneTable90 <- merge(onlyZoneTable90,tax_tablemix, by = "ASV_Id")
##Table
write.table(onlyZoneTable90, file = "Table/Only_Zone90.txt", sep = "\t", col.names = TRUE, row.names = FALSE,quote = FALSE)
write.table(onlyZoneSequence90, file = "Table/Only_Zone_Sequence90.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  # Fraction-90 -------------------------------------------------------------------
##Petite
onlyPetiteTable90 <- data_seq_tax %>% filter(Fraction90 == "Petite") %>% select(ASV_Id,TotalPetite,Division,Fraction90)
onlyPetiteTable90 <- onlyPetiteTable90 %>% filter(TotalPetite > 0.01*sum(onlyPetiteTable90$TotalPetite))
colnames(onlyPetiteTable90)[2]  <- "value"
##Grande
onlyGrandeTable90 <- data_seq_tax %>% filter(Fraction90 == "Grande") %>% select(ASV_Id,TotalGrande,Division,Fraction90)
onlyGrandeTable90 <- onlyGrandeTable90 %>% filter(TotalGrande > 0.01*sum(onlyGrandeTable90$TotalGrande))
colnames(onlyGrandeTable90)[2]  <- "value"
##bind
onlyFractionTable90 <- rbind(onlyPetiteTable90,onlyGrandeTable90)
onlyFractionTable90 <- merge(onlyFractionTable90,tax_tablemix, by = "ASV_Id")
##Table
write.table(onlyFractionTable90, file = "Table/Only_Fraction90.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(onlyFractionSequence90, file = "Table/Only_Fraction_Sequence90.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Table Total---------------------------------------------------------------------
  # Cycle Total -------------------------------------------------------------------
##Jour
totalJourTable <- data_seq_tax %>% select(ASV_Id,TotalJour,Division,Cycle)
totalJourTable$Cycle <- "Jour"
totalJourTable <- totalJourTable %>% filter(TotalJour > 0.05*sum(totalJourTable$TotalJour))
colnames(totalJourTable)[2]  <- "value"
##Nuit
totalNuitTable <- data_seq_tax %>% select(ASV_Id,TotalNuit,Division,Cycle)
totalNuitTable$Cycle <- "Nuit"
totalNuitTable <- totalNuitTable %>% filter(TotalNuit > 0.05*sum(totalNuitTable$TotalNuit))
colnames(totalNuitTable)[2]  <- "value"
##Bind
totalCycleTable <- rbind(totalJourTable,totalNuitTable)
totalCycleTable <- merge(totalCycleTable,tax_tablemix, by = "ASV_Id")
##Table
write.table(totalCycleTable, file = "Table/Total_Cycles.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(totalCycleASV, file = "Table/Total_Cycles_ASV.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(totalCycleSequence, file = "Table/Total_Cycles_Sequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  # Zone Total -------------------------------------------------------------------
##Oxique
totalOxiqueTable <- data_seq_tax %>% select(ASV_Id,TotalOxique,Division,Zone)
totalOxiqueTable$Zone <- rep("Oxique",each = nrow(totalOxiqueTable))
totalOxiqueTable <- totalOxiqueTable %>% filter(TotalOxique > 0.05*sum(totalOxiqueTable$TotalOxique))
colnames(totalOxiqueTable)[2]  <- "value"
##Anoxique
totalAnoxiqueTable <- data_seq_tax %>% select(ASV_Id,TotalAnoxique,Division,Zone)
totalAnoxiqueTable$Zone <- rep("Anoxique",each = nrow(totalAnoxiqueTable))
totalAnoxiqueTable <- totalAnoxiqueTable %>% filter(TotalAnoxique > 0.05*sum(totalAnoxiqueTable$TotalAnoxique))
colnames(totalAnoxiqueTable)[2]  <- "value"
##Bind
totalZoneTable <- rbind(totalOxiqueTable,totalAnoxiqueTable)
totalZoneTable <- merge(totalZoneTable,tax_tablemix, by = "ASV_Id")
##Table
write.table(totalZoneTable, file = "Table/Total_Zone.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(totalZoneASV, file = "Table/Total_Zone_ASV.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(totalZoneSequence, file = "Table/Total_Zone_Sequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  # Fraction Total -------------------------------------------------------------------
##Petite
totalPetiteTable <- data_seq_tax %>% select(ASV_Id,TotalPetite,Division,Fraction)
totalPetiteTable$Fraction <- rep("Petite",each = nrow(totalPetiteTable))
totalPetiteTable <- totalPetiteTable %>% filter(TotalPetite > 0.05*sum(totalPetiteTable$TotalPetite))
colnames(totalPetiteTable)[2]  <- "value"
##Grande
totalGrandeTable <- data_seq_tax %>% select(ASV_Id,TotalGrande,Division,Fraction)
totalGrandeTable$Fraction <- rep("Grande",each = nrow(totalGrandeTable))
totalGrandeTable <- totalGrandeTable %>% filter(TotalGrande > 0.05*sum(totalGrandeTable$TotalGrande))
colnames(totalGrandeTable)[2]  <- "value"
##Bind
totalFractionTable <- rbind(totalPetiteTable,totalGrandeTable)
totalFractionTable <- merge(totalFractionTable,tax_tablemix, by = "ASV_Id")
##Table
write.table(totalFractionTable, file = "Table/Total_Fraction.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(totalFractionASV, file = "Table/Total_Fraction_ASV.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(totalFractionSequence, file = "Table/Total_Fraction_Sequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  # Date Total -------------------------------------------------------------------
total04Table <- data_seq_tax %>% select(ASV_Id,Total04,Division)
total04Table$Dates <- rep(04,each = nrow(total04Table))
total04Table <- total04Table %>% filter(Total04 > 0.05*sum(total04Table$Total04))
colnames(total04Table)[2]  <- "value"
##
total06Table <- data_seq_tax %>% select(ASV_Id,Total06,Division)
total06Table$Dates <- rep(06,each = nrow(total06Table))
total06Table <- total06Table %>% filter(Total06 > 0.05*sum(total06Table$Total06))
colnames(total06Table)[2]  <- "value"
##
total09Table <- data_seq_tax %>% select(ASV_Id,Total09,Division)
total09Table$Dates <- rep(09,each = nrow(total09Table))
total09Table <- total09Table %>% filter(Total09 > 0.05*sum(total09Table$Total09))
colnames(total09Table)[2]  <- "value"
##
total11Table <- data_seq_tax %>% select(ASV_Id,Total11,Division)
total11Table$Dates <- rep(11,each = nrow(total11Table))
total11Table <- total11Table %>% filter(Total11 > 0.05*sum(total11Table$Total11))
colnames(total11Table)[2]  <- "value"
##
totalDateTable <- rbind(total04Table,total06Table,total09Table,total11Table)
totalDateTable <- merge(totalDateTable,tax_tablemix, by = "ASV_Id")
##
write.table(totalDateTable, file = "Table/Total_Dates.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(totalDateASV, file = "Table/Total_Dates_ASV.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(totalDateSequence, file = "Table/Total_Dates_Sequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)



