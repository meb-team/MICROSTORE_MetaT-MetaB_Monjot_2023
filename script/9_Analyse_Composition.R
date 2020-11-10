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
#
#output <- "Analyse-Composition-Rarefy-V4-095-199-NOfilter-Total"
#input <- "../dataPANAM/PANAM2/V4-result-095-199/OTU_distribution_tax.txt"
#region <- "V4"
#sortop <- "no"
#Taxonomy <- "NN"
#Mode <- "Superphylum"
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

# Import package and palette -----------------------------------------------------------
pkg <- c("ggplot2", "readxl","dplyr","tidyr","cowplot","FactoMineR","factoextra","reshape2","varhandle","ggrepel","ggpubr","ggsci","scales","hrbrthemes","GUniFrac","svglite")
lapply(pkg, require, character.only = TRUE)

palette <-  sample(c(pal_locuszoom(alpha = 0.8)(7), pal_lancet(alpha = 0.8)(9)))
show_col(palette)
#palette <- c("#D43F3ACC","#EEA236CC","#AD002ACC","#46B8DACC","#357EBDCC","#9632B8CC","#B8B8B8CC","#00468BCC","#ED0000CC","#42B540CC","#0099B4CC","#925E9FCC")

input <- paste("..",input, sep = "/")
output <- paste("../result",output, sep = "/")
if (dir.exists("../result") == FALSE) { dir.create("../result") }
if (dir.exists(output) == FALSE) { dir.create(output) }
if (dir.exists(output) == TRUE) { setwd(output) }

#Enter mode and Group to study

if (length(args)==2) {
  cat("Enter mode (Superphylum or Phylum) : ");
  Mode <- readLines("stdin",n=1);
  cat("You entered")
  str(Mode);
  cat( "\n" )}
if (length(args)==6) {
  Mode <- args[6]
}
if (Mode == "Phylum") {
  cat("Enter Phylum (Fungi, Alveolata, etc) : ");
  Group <- readLines("stdin",n=1);
  cat("You entered")
  str(Group);
  cat( "\n" )}
if (length(args)==7) {
  Group <- args[7]
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


# Input OTU Table ---------------------------------------------------------
tableVinput <- read.csv(file = input, sep = "\t", row.names = "OTU_Id")
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
  # Prepare seq_mat and Pool (seq_mat and seq_mat_pool) ---------------------------------------------------------------------
## Prepare otu_mat
amplicon <- c(grep(pattern = "OSTA", colnames(tableVinput), value = TRUE))
seq_mat <- tableVinput %>% select(all_of(amplicon))
seq_mat[,all_of(amplicon)] <- lapply(seq_mat[,all_of(amplicon)], as.numeric)

## Prepare taxonomy
if (length(args)==2) {
  cat("Enter Taxonomy mod (LCA, NN or Best_HIT) ? : ");
  Taxonomy <- readLines("stdin",n=1);
  cat("You entered")
  str(Taxonomy);
  cat( "\n" )}
if (length(args)==5) {
  Taxonomy <- args[5]
}

if (Taxonomy == "LCA") { Taxonomy <- "LCA.taxonomy"}
if (Taxonomy == "NN") { Taxonomy <- "NN.taxonomy"}
if (Taxonomy == "Best_HIT") { Taxonomy <- "Best_hit_identity"}
seq_mat[,"Taxonomy"] <- tableVinput[,Taxonomy]

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

## Pool
seq_mat_pool <- as.data.frame(x=seq_mat[,"Taxonomy"],row.names = row.names(seq_mat))
colnames(seq_mat_pool) <- "Taxonomy"
for (h in tblcx[,"Amplicon.x"] ) {
  f <- (tblcx %>% filter(Amplicon.x == h) %>% select(Amplicon.y))[[1]]
  print (h)
  print (f)
  for (r in row.names(seq_mat)) {
    if (seq_mat[r,h] * seq_mat[r,f] == 0) { seq_mat_pool[r,h] <- round((seq_mat[r,h]+seq_mat[r,f]))}
    if (seq_mat[r,h] * seq_mat[r,f] != 0) { seq_mat_pool[r,h] <- round((seq_mat[r,h]+seq_mat[r,f])/2)}
  }}

  # Rarefy (seq_mat_pool_rare) ---------------------------------------------------------------------
seq_mat_poolt <- seq_mat_pool %>% select(-"Taxonomy")
seq_mat_poolt <- t(seq_mat_poolt)
seq_mat_pool_rare <- as.data.frame(t(Rarefy(seq_mat_poolt)$otu.tab.rff)) ## Ref : Jun Chen et al. (2012). Associating microbiome composition with environmental covariates using generalized UniFrac distances. 28(16): 2106–2113.
rm(seq_mat_poolt)
for (i in row.names(seq_mat_pool_rare)) { seq_mat_pool_rare[i,"Taxonomy"] <- seq_mat_pool[i,"Taxonomy"]}
  # Prepare otu_mat (otu_mat, otu_mat_pool and otu_mat_pool_rare) ------------------------------------------------------------------
## otu_mat
otu_mat <- seq_mat
amplicon <- c(grep(pattern = "OSTA", colnames(seq_mat), value = TRUE))
otu_mat[,all_of(amplicon)][otu_mat[,all_of(amplicon)] != 0] <- 1
## otu_mat_pool
otu_mat_pool <- seq_mat_pool
amplicon <- c(grep(pattern = "OSTA", colnames(otu_mat_pool), value = TRUE))
otu_mat_pool[,all_of(amplicon)][otu_mat_pool[,all_of(amplicon)] != 0] <- 1
## otu_mat_pool_rare
otu_mat_pool_rare <- seq_mat_pool_rare
amplicon <- c(grep(pattern = "OSTA", colnames(otu_mat_pool_rare), value = TRUE))
otu_mat_pool_rare[,all_of(amplicon)][otu_mat_pool_rare[,all_of(amplicon)] != 0] <- 1

# Stat Rarefy & Pool -------------------------------------------------------------
## Sequence
amplicon <- c(grep(pattern = "OSTA", colnames(seq_mat_pool), value = TRUE))
avRarefyS <- as.data.frame(colSums(seq_mat_pool[,all_of(amplicon)]))
amplicon <- c(grep(pattern = "OSTA", colnames(seq_mat_pool_rare), value = TRUE))
apRarefyS <- as.data.frame(colSums(seq_mat_pool_rare[,all_of(amplicon)]))
statRarefy <- cbind(avRarefyS,apRarefyS)
## OTU
amplicon <- c(grep(pattern = "OSTA", colnames(otu_mat_pool), value = TRUE))
avRarefyO <- as.data.frame(colSums(otu_mat_pool[,all_of(amplicon)]))
amplicon <- c(grep(pattern = "OSTA", colnames(otu_mat_pool_rare), value = TRUE))
apRarefyO <- as.data.frame(colSums(otu_mat_pool_rare[,all_of(amplicon)]))
statRarefy <- cbind(statRarefy,avRarefyO,apRarefyO)
## Totaux
statRarefy["Total",]<-colSums(statRarefy)
colnames(statRarefy) <- c("avRarefy-Sequence","apRarefy-Sequence","avRarefy-OTU","apRarefy-OTU")
## Remove empty OTUs
### Avant raréfaction
OTUav <- otu_mat_pool
OTUav$OTU_Id <- row.names(OTUav)
amplicon <- c(grep(pattern = "OSTA", colnames(OTUav), value = TRUE))
for ( i in row.names(OTUav)) { if (sum(OTUav[i,all_of(amplicon)]) == 0) { OTUav[i,"OTU_Id"] <- "Uniq"}}
OTUav <- OTUav %>% filter(OTU_Id != "Uniq")
### Après rarefaction
OTUap <- otu_mat_pool_rare
OTUap$OTU_Id <- row.names(OTUap)
amplicon <- c(grep(pattern = "OSTA", colnames(OTUap), value = TRUE))
for ( i in row.names(OTUap)) { if (sum(OTUap[i,all_of(amplicon)]) == 0) { OTUap[i,"OTU_Id"] <- "Uniq"}}
OTUap <- OTUap %>% filter(OTU_Id != "Uniq")
### Final
statRarefy["Total","avRarefy-OTU"] <- nrow(OTUav)
statRarefy["Total","apRarefy-OTU"] <- nrow(OTUap)
write.table(statRarefy, file = "TableOnly/StatRarefy_withoutDuplicat.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

  # Sum files ---------------------------------------------------------------
## Rarefaction stat Fig
statRarefyOTU <- statRarefy %>% select("avRarefy-OTU","apRarefy-OTU")
colnames(statRarefyOTU) <- c("NON","OUI")
statRarefyOTU$Echantillon <- row.names(statRarefyOTU)
amplicon <- c(grep(pattern = "OSTA", colnames(OTUap), value = TRUE))
statRarefyOTU <- melt(statRarefyOTU[all_of(amplicon),], id = "Echantillon")
my_comp <- list(c("OUI","NON"))

svglite("Figure-Sum/Stat-Rarefy.svg",width = 3.00,height = 4.00)
R <- ggplot(statRarefyOTU, aes(y = value, x = variable)) + geom_boxplot() + geom_point(aes(color = Echantillon), size = 2.5) + 
  stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
  stat_compare_means(method="wilcox.test", label.y = max(statRarefyOTU$value)+0.12*max(statRarefyOTU$value)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  labs(x="Raréfaction",y="Nombre d'OTUs") #+ guides(color = FALSE)
R <- R + theme(legend.position="none")
print(R)
dev.off() 

## Pooling stat Fig
statRarefy <- statRarefy[!(row.names(statRarefy) %in% "Total"), ]
statRarefy[,"Fusion"] <- "OUI"
statRarefy[Amplicon_without_duplicat,"Fusion"] <- "NON"
statRarefy$Echantillon <- rownames(statRarefy)
my_comp <- list(c("OUI","NON"))
### OTU après Rarefy
svglite("Figure-Sum/Analyse-SumApRare-OTU.svg",width = 3.00,height = 4.00)
J <- ggplot(statRarefy, aes(y = `apRarefy-OTU`, x = Fusion)) + geom_boxplot() + geom_point(aes(color = Echantillon), size = 2.5) + 
  stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
  stat_compare_means(method="wilcox.test", label.y = max(statRarefy$`apRarefy-OTU`)+0.12*max(statRarefy$`apRarefy-OTU`)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  labs(x="Regroupement",y="Nombre d'OTUs") #+ guides(color = FALSE)
J <- J + theme(legend.position="none")
print(J)
dev.off()  
### Séquence avant Rarefy
svglite("Figure-Sum/Analyse-SumAvRare-Sequences.svg",width = 3.00,height = 4.00)
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
### OTU avant Rarefy
svglite("Figure-Sum/Analyse-SumAvRare-OTU.svg",width = 3.00,height = 4.00)
H <- ggplot(statRarefy, aes(y = `avRarefy-OTU`, x = Fusion)) + geom_boxplot() + geom_point(aes(color = Echantillon), size = 2.5) + 
  stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
  stat_compare_means(method="wilcox.test", label.y = max(statRarefy$`avRarefy-OTU`)+0.12*max(statRarefy$`avRarefy-OTU`)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  labs(x="Regroupement",y="Nombre d'OTUs") #+ guides(color = FALSE)
H <- H + theme(legend.position="none")
print(H)
dev.off()  
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
Jour_seq$OTU_Id <- row.names(Jour_seq)
##Nuit
Nuit_seq <- as.data.frame(rowSums(seq_mat_pool_rare %>% select(all_of(CycleNuit))))
colnames(Nuit_seq) <- "TotalNuit"
Nuit_seq$OTU_Id <- row.names(Nuit_seq)
##Cycle
Cycle_seq <- merge(x = Jour_seq,y = Nuit_seq, by = "OTU_Id")

#Zone
##Oxique
Oxique_seq <- as.data.frame(rowSums(seq_mat_pool_rare %>% select(all_of(ZoneOxique))))
colnames(Oxique_seq) <- "TotalOxique"
Oxique_seq$OTU_Id <- row.names(Oxique_seq)
##Anoxique
Anoxique_seq <- as.data.frame(rowSums(seq_mat_pool_rare %>% select(all_of(ZoneAnoxique))))
colnames(Anoxique_seq) <- "TotalAnoxique"
Anoxique_seq$OTU_Id <- row.names(Anoxique_seq)
##Zone
Zone_seq <- merge(x = Oxique_seq,y = Anoxique_seq, by = "OTU_Id")

#Fraction
##Petite
Petite_seq <- as.data.frame(rowSums(seq_mat_pool_rare %>% select(all_of(FractionPetite))))
colnames(Petite_seq) <- "TotalPetite"
Petite_seq$OTU_Id <- row.names(Petite_seq)
##Grande
Grande_seq <- as.data.frame(rowSums(seq_mat_pool_rare %>% select(all_of(FractionGrande))))
colnames(Grande_seq) <- "TotalGrande"
Grande_seq$OTU_Id <- row.names(Grande_seq)
##Fraction
Fraction_seq <- merge(x = Petite_seq,y = Grande_seq, by = "OTU_Id")

#Date
##04
`04_seq` <- as.data.frame(rowSums(seq_mat_pool_rare %>% select(all_of(Date04))))
colnames(`04_seq`) <- "Total04"
`04_seq`$OTU_Id <- row.names(`04_seq`)
##06
`06_seq` <- as.data.frame(rowSums(seq_mat_pool_rare %>% select(all_of(Date06))))
colnames(`06_seq`) <- "Total06"
`06_seq`$OTU_Id <- row.names(`06_seq`)
##09
`09_seq` <- as.data.frame(rowSums(seq_mat_pool_rare %>% select(all_of(Date09))))
colnames(`09_seq`) <- "Total09"
`09_seq`$OTU_Id <- row.names(`09_seq`)
##11
`11_seq` <- as.data.frame(rowSums(seq_mat_pool_rare %>% select(all_of(Date11))))
colnames(`11_seq`) <- "Total11"
`11_seq`$OTU_Id <- row.names(`11_seq`)
##Date
Date_seq <- merge(x = `04_seq`,y = `06_seq`, by = "OTU_Id")
Date_seq <- merge(x = Date_seq,y = `09_seq`, by = "OTU_Id")
Date_seq <- merge(x = Date_seq,y = `11_seq`, by = "OTU_Id")

  # OTU ----------------------------------------------------------------
#Cycle
##Jour
Jour_otu <- as.data.frame(rowSums(otu_mat_pool_rare %>% select(all_of(CycleJour))))
colnames(Jour_otu) <- "TotalJour"
for (i in row.names(Jour_otu)) {if (Jour_otu[i,"TotalJour"] > 0) {Jour_otu[i,"TotalJour"] <- 1}}
Jour_otu$OTU_Id <- row.names(Jour_otu)
##Nuit
Nuit_otu <- as.data.frame(rowSums(otu_mat_pool_rare %>% select(all_of(CycleNuit))))
colnames(Nuit_otu) <- "TotalNuit"
for (i in row.names(Nuit_otu)) {if (Nuit_otu[i,"TotalNuit"] > 0) {Nuit_otu[i,"TotalNuit"] <- 1}}
Nuit_otu$OTU_Id <- row.names(Nuit_otu)
##Cycle
Cycle_otu <- merge(x = Jour_otu,y = Nuit_otu, by = "OTU_Id")

#Zone
##Oxique
Oxique_otu <- as.data.frame(rowSums(otu_mat_pool_rare %>% select(all_of(ZoneOxique))))
colnames(Oxique_otu) <- "TotalOxique"
for (i in row.names(Oxique_otu)) {if (Oxique_otu[i,"TotalOxique"] > 0) {Oxique_otu[i,"TotalOxique"] <- 1}}
Oxique_otu$OTU_Id <- row.names(Oxique_otu)
##Anoxique
Anoxique_otu <- as.data.frame(rowSums(otu_mat_pool_rare %>% select(all_of(ZoneAnoxique))))
colnames(Anoxique_otu) <- "TotalAnoxique"
for (i in row.names(Anoxique_otu)) {if (Anoxique_otu[i,"TotalAnoxique"] > 0) {Anoxique_otu[i,"TotalAnoxique"] <- 1}}
Anoxique_otu$OTU_Id <- row.names(Anoxique_otu)
##Zone
Zone_otu <- merge(x = Oxique_otu,y = Anoxique_otu, by = "OTU_Id")

#Fraction
##Petite
Petite_otu <- as.data.frame(rowSums(otu_mat_pool_rare %>% select(all_of(FractionPetite))))
colnames(Petite_otu) <- "TotalPetite"
for (i in row.names(Petite_otu)) {if (Petite_otu[i,"TotalPetite"] > 0) {Petite_otu[i,"TotalPetite"] <- 1}}
Petite_otu$OTU_Id <- row.names(Petite_otu)
##Grande
Grande_otu <- as.data.frame(rowSums(otu_mat_pool_rare %>% select(all_of(FractionGrande))))
colnames(Grande_otu) <- "TotalGrande"
for (i in row.names(Grande_otu)) {if (Grande_otu[i,"TotalGrande"] > 0) {Grande_otu[i,"TotalGrande"] <- 1}}
Grande_otu$OTU_Id <- row.names(Grande_otu)
##Fraction
Fraction_otu <- merge(x = Petite_otu,y = Grande_otu, by = "OTU_Id")

#Date
##04
`04_otu` <- as.data.frame(rowSums(otu_mat_pool_rare %>% select(all_of(Date04))))
colnames(`04_otu`) <- "Total04"
for (i in row.names(`04_otu`)) {if (`04_otu`[i,"Total04"] > 0) {`04_otu`[i,"Total04"] <- 1}}
`04_otu`$OTU_Id <- row.names(`04_otu`)
##06
`06_otu` <- as.data.frame(rowSums(otu_mat_pool_rare %>% select(all_of(Date06))))
colnames(`06_otu`) <- "Total06"
for (i in row.names(`06_otu`)) {if (`06_otu`[i,"Total06"] > 0) {`06_otu`[i,"Total06"] <- 1}}
`06_otu`$OTU_Id <- row.names(`06_otu`)
##09
`09_otu` <- as.data.frame(rowSums(otu_mat_pool_rare %>% select(all_of(Date09))))
colnames(`09_otu`) <- "Total09"
for (i in row.names(`09_otu`)) {if (`09_otu`[i,"Total09"] > 0) {`09_otu`[i,"Total09"] <- 1}}
`09_otu`$OTU_Id <- row.names(`09_otu`)
##11
`11_otu` <- as.data.frame(rowSums(otu_mat_pool_rare %>% select(all_of(Date11))))
colnames(`11_otu`) <- "Total11"
for (i in row.names(`11_otu`)) {if (`11_otu`[i,"Total11"] > 0) {`11_otu`[i,"Total11"] <- 1}}
`11_otu`$OTU_Id <- row.names(`11_otu`)
##Date
Date_otu <- merge(x = `04_otu`,y = `06_otu`, by = "OTU_Id")
Date_otu <- merge(x = Date_otu,y = `09_otu`, by = "OTU_Id")
Date_otu <- merge(x = Date_otu,y = `11_otu`, by = "OTU_Id")

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
coord$OTU_Id <- row.names(coord)
##Cycle
data_seq <- merge(x = coord, y = Cycle_seq, by = "OTU_Id")
data_seq$Jour <- 0
data_seq$Nuit <- 0
for (i in rownames(data_seq)) { if ( data_seq[i,"TotalJour"] != 0) { data_seq[i,"Jour"] <- 1}}
for (i in rownames(data_seq)) { if ( data_seq[i,"TotalNuit"] != 0) { data_seq[i,"Nuit"] <- 2}}
data_seq$Cycle <- 0
for (i in rownames(data_seq)) { data_seq[i,"Cycle"] <- data_seq[i,"Jour"] + data_seq[i,"Nuit"] 
if ( data_seq[i,"Cycle"] == 1) { data_seq[i,"Cycle"] <- "Jour"}
if ( data_seq[i,"Cycle"] == 2) { data_seq[i,"Cycle"] <- "Nuit"}
if ( data_seq[i,"Cycle"] == 3) { data_seq[i,"Cycle"] <- "Communs"}}
##Zone
data_seq <- merge(x = data_seq, y = Zone_seq, by = "OTU_Id")
data_seq$Oxique <- 0
data_seq$Anoxique <- 0
for (i in rownames(data_seq)) { if ( data_seq[i,"TotalOxique"] != 0) { data_seq[i,"Oxique"] <- 1}}
for (i in rownames(data_seq)) { if ( data_seq[i,"TotalAnoxique"] != 0) { data_seq[i,"Anoxique"] <- 2}}
data_seq$Zone <- 0
for (i in rownames(data_seq)) { data_seq[i,"Zone"] <- data_seq[i,"Oxique"] + data_seq[i,"Anoxique"] 
if ( data_seq[i,"Zone"] == 1) { data_seq[i,"Zone"] <- "Oxique"}
if ( data_seq[i,"Zone"] == 2) { data_seq[i,"Zone"] <- "Anoxique"}
if ( data_seq[i,"Zone"] == 3) { data_seq[i,"Zone"] <- "Communs"}}
##Fraction
data_seq <- merge(x = data_seq, y = Fraction_seq, by = "OTU_Id")
data_seq$Petite <- 0
data_seq$Grande <- 0
for (i in rownames(data_seq)) { if ( data_seq[i,"TotalPetite"] != 0) { data_seq[i,"Petite"] <- 1}}
for (i in rownames(data_seq)) { if ( data_seq[i,"TotalGrande"] != 0) { data_seq[i,"Grande"] <- 2}}
data_seq$Fraction <- 0
for (i in rownames(data_seq)) { data_seq[i,"Fraction"] <- data_seq[i,"Petite"] + data_seq[i,"Grande"] 
if ( data_seq[i,"Fraction"] == 1) { data_seq[i,"Fraction"] <- "Petite"}
if ( data_seq[i,"Fraction"] == 2) { data_seq[i,"Fraction"] <- "Grande"}
if ( data_seq[i,"Fraction"] == 3) { data_seq[i,"Fraction"] <- "Communs"}}
#Date
data_seq <- merge(x = data_seq, y = Date_seq, by = "OTU_Id")
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

  # OTU ----------------------------------------------------------------
dt <- as.data.frame(otu_mat_pool_rare)
dt <- dt %>% select(-"Taxonomy")
#CA
res.ca <- CA(dt, graph = FALSE,ncp = 2 )
fviz_ca_row(res.ca, repel = FALSE, label = "none")
p <- get_ca_row(res.ca)
coord <- p$coord
coord <- as.data.frame(coord)
coord$OTU_Id <- row.names(coord)
##Cycle
data_otu <- merge(x = coord, y = Cycle_otu, by = "OTU_Id")
data_otu$Jour <- 0
data_otu$Nuit <- 0
for (i in rownames(data_otu)) { if ( data_otu[i,"TotalJour"] != 0) { data_otu[i,"Jour"] <- 1}}
for (i in rownames(data_otu)) { if ( data_otu[i,"TotalNuit"] != 0) { data_otu[i,"Nuit"] <- 2}}
data_otu$Cycle <- 0
for (i in rownames(data_otu)) { data_otu[i,"Cycle"] <- data_otu[i,"Jour"] + data_otu[i,"Nuit"] 
if ( data_otu[i,"Cycle"] == 1) { data_otu[i,"Cycle"] <- "Jour"}
if ( data_otu[i,"Cycle"] == 2) { data_otu[i,"Cycle"] <- "Nuit"}
if ( data_otu[i,"Cycle"] == 3) { data_otu[i,"Cycle"] <- "Communs"}}
##Zone
data_otu <- merge(x = data_otu, y = Zone_otu, by = "OTU_Id")
data_otu$Oxique <- 0
data_otu$Anoxique <- 0
for (i in rownames(data_otu)) { if ( data_otu[i,"TotalOxique"] != 0) { data_otu[i,"Oxique"] <- 1}}
for (i in rownames(data_otu)) { if ( data_otu[i,"TotalAnoxique"] != 0) { data_otu[i,"Anoxique"] <- 2}}
data_otu$Zone <- 0
for (i in rownames(data_otu)) { data_otu[i,"Zone"] <- data_otu[i,"Oxique"] + data_otu[i,"Anoxique"] 
if ( data_otu[i,"Zone"] == 1) { data_otu[i,"Zone"] <- "Oxique"}
if ( data_otu[i,"Zone"] == 2) { data_otu[i,"Zone"] <- "Anoxique"}
if ( data_otu[i,"Zone"] == 3) { data_otu[i,"Zone"] <- "Communs"}}
##Fraction
data_otu <- merge(x = data_otu, y = Fraction_otu, by = "OTU_Id")
data_otu$Petite <- 0
data_otu$Grande <- 0
for (i in rownames(data_otu)) { if ( data_otu[i,"TotalPetite"] != 0) { data_otu[i,"Petite"] <- 1}}
for (i in rownames(data_otu)) { if ( data_otu[i,"TotalGrande"] != 0) { data_otu[i,"Grande"] <- 2}}
data_otu$Fraction <- 0
for (i in rownames(data_otu)) { data_otu[i,"Fraction"] <- data_otu[i,"Petite"] + data_otu[i,"Grande"] 
if ( data_otu[i,"Fraction"] == 1) { data_otu[i,"Fraction"] <- "Petite"}
if ( data_otu[i,"Fraction"] == 2) { data_otu[i,"Fraction"] <- "Grande"}
if ( data_otu[i,"Fraction"] == 3) { data_otu[i,"Fraction"] <- "Communs"}}
##Date
data_otu <- merge(x = data_otu, y = Date_otu, by = "OTU_Id")
data_otu04 <- data_otu %>% filter(Total04 != 0) %>% select(-Total06,-Total09,-Total11)
data_otu06 <- data_otu %>% filter(Total06 != 0) %>% select(-Total04,-Total09,-Total11)
data_otu09 <- data_otu %>% filter(Total09 != 0) %>% select(-Total04,-Total06,-Total11)
data_otu11 <- data_otu %>% filter(Total11 != 0) %>% select(-Total04,-Total06,-Total09)

Xotu <- fviz_screeplot (res.ca, addlabels = TRUE, ylim = c(0, 50))
Xotu <- Xotu$data
Dim1otu <- paste("Dim 1 [",round(Xotu %>% filter(dim == 1) %>% select(eig),1),"%]", sep = "")
Dim2otu <- paste("Dim 2 [",round(Xotu %>% filter(dim == 2) %>% select(eig),1),"%]", sep = "")




# AFC Taxonomy ------------------------------------------------------------
  # Séquence ----------------------------------------------------------------
data_seq_tax <- data_seq
taxo <- tableVinput %>% select(all_of(Taxonomy)) ; taxo$OTU_Id <- row.names(taxo)
data_seq_tax <- merge(x = data_seq_tax, y =  taxo, by = "OTU_Id")
data_seq_tax <- separate(data_seq_tax, all_of(Taxonomy),c("Domain","Superphylum","Phylum","Class","Order","Family","Genus","Species"), sep =";")

if (Mode == "Phylum") {
data_seq_tax <- data_seq_tax %>% filter(Superphylum == all_of(Group))
data_seq_tax$Division <- paste(data_seq_tax$Domain, data_seq_tax$Superphylum, data_seq_tax$Phylum,sep=";")}
if (Mode == "Superphylum") {
data_seq_tax$Division <- paste(data_seq_tax$Domain, data_seq_tax$Superphylum,sep=";")}
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
svglite("AFC-Distribution/AFC-Sequence.svg",width = 12.00,height = 8.00)
All_plot <- plot_grid(aCycle, cFraction, bZone, ncol = 2, nrow = 2, rel_widths = c(3,3), rel_heights = c(3,3))
print(All_plot)
dev.off()
    # Avril_2018 ---------------------------------------------------
data_seq04_tax <- merge(x = data_seq04, y =  taxo, by = "OTU_Id")

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

svglite("AFC-Distribution/AFC-Sequence-Avril.svg",width = 12.00,height = 8.00)
Avril_plot <- plot_grid(a_04, c_04, b_04, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(Avril_plot)
dev.off()  
    # Juin_2018 ---------------------------------------------------
data_seq06_tax <- merge(x = data_seq06, y =  taxo, by = "OTU_Id")

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

svglite("AFC-Distribution/AFC-Sequence-Juin.svg",width = 12.00,height = 8.00)
Juin_plot <- plot_grid(a_06, c_06, b_06, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(Juin_plot)
dev.off()  

    # Septembre_2018 ---------------------------------------------------
data_seq09_tax <- merge(x = data_seq09, y =  taxo, by = "OTU_Id")

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

svglite("AFC-Distribution/AFC-Sequence-Septembre.svg",width = 12.00,height = 8.00)
Septembre_plot <- plot_grid(a_09, c_09, b_09, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(Septembre_plot)
dev.off()  


    # Novembre_2018 ---------------------------------------------------
data_seq11_tax <- merge(x = data_seq11, y =  taxo, by = "OTU_Id")

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

svglite("AFC-Distribution/AFC-Sequence-Novembre.svg",width = 12.00,height = 8.00)
Novembre_plot <- plot_grid(a_11, c_11, b_11, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(Novembre_plot)
dev.off()  


  # OTU ----------------------------------------------------------------
data_otu_tax <- data_otu
taxo <- tableVinput %>% select(all_of(Taxonomy)) ; taxo$OTU_Id <- row.names(taxo)
data_otu_tax <- merge(x = data_otu_tax, y =  taxo, by = "OTU_Id")
data_otu_tax <- separate(data_otu_tax, all_of(Taxonomy),c("Domain","Superphylum","Phylum","Class","Order","Family","Genus","Species"), sep =";")
if (Mode == "Phylum") {
data_otu_tax <- data_otu_tax %>% filter(Superphylum == all_of(Group))
data_otu_tax$Division <- paste(data_otu_tax$Domain, data_otu_tax$Superphylum, data_otu_tax$Phylum,sep=";")}
if (Mode == "Superphylum") {
data_otu_tax$Division <- paste(data_otu_tax$Domain, data_otu_tax$Superphylum,sep=";")}
    # All_2018 ----------------------------------------------------------
##Figure Cycles
iCycle <- ggplot(data_otu_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Cycle)) + geom_point(size = 2) +
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
  labs(x=Dim1otu,y=Dim2otu,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(iCycle)

##Figure Zone
jZone <- ggplot(data_otu_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Zone)) + geom_point(size = 2) +
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
  labs(x=Dim1otu,y=Dim2otu,color = "Zones", alpha = "Zones", linetype = "Zones")
print(jZone)

##Figure Fraction
kFraction <- ggplot(data_otu_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Fraction)) + geom_point(size = 2) +
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
  labs(x=Dim1otu,y=Dim2otu,color = "Fractions", alpha = "Fractions", linetype = "Fractions")
print(kFraction)
##Coplot
svglite("AFC-Distribution/AFC-OTU.svg",width = 12.00,height = 8.00)
All_plot <- plot_grid(iCycle, jZone, kFraction, ncol = 2, nrow = 2, rel_widths = c(3,3), rel_heights = c(3,3))
print(All_plot)
dev.off()
    # Avril_2018 ---------------------------------------------------
data_otu04_tax <- merge(x = data_otu04, y =  taxo, by = "OTU_Id")

i_04 <- ggplot(data_otu04_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Cycle)) + geom_point(size = 2) +
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
  labs(x=Dim1otu,y=Dim2otu,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(i_04)

#Figure Zone
j_04 <- ggplot(data_otu04_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Zone)) + geom_point(size = 2) +
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
  labs(x=Dim1otu,y=Dim2otu,color = "Zones", alpha = "Zones", linetype = "Zones")
print(j_04)

#Figure Fraction
k_04 <- ggplot(data_otu04_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Fraction)) + geom_point(size = 2) +
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
  labs(x=Dim1otu,y=Dim2otu,color = "Fraction", alpha = "Fraction", linetype = "Fraction")
print(k_04)

svglite("AFC-Distribution/AFC-OTU-Avril.svg",width = 12.00,height = 8.00)
Avril_plot <- plot_grid(i_04, j_04, k_04, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(Avril_plot)
dev.off()  
    # Juin_2018 ---------------------------------------------------
data_otu06_tax <- merge(x = data_otu06, y =  taxo, by = "OTU_Id")

i_06 <- ggplot(data_otu06_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Cycle)) + geom_point(size = 2) +
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
  labs(x=Dim1otu,y=Dim2otu,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(i_06)

#Figure Zone
j_06 <- ggplot(data_otu06_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Zone)) + geom_point(size = 2) +
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
  labs(x=Dim1otu,y=Dim2otu,color = "Zones", alpha = "Zones", linetype = "Zones")
print(j_06)

#Figure Fraction
k_06 <- ggplot(data_otu06_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Fraction)) + geom_point(size = 2) +
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
  labs(x=Dim1otu,y=Dim2otu,color = "Fraction", alpha = "Fraction", linetype = "Fraction")
print(k_06)

svglite("AFC-Distribution/AFC-OTU-Juin.svg",width = 12.00,height = 8.00)
Juin_plot <- plot_grid(i_06, j_06, k_06, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(Juin_plot)
dev.off()  

    # Septembre_2018 ---------------------------------------------------
data_otu09_tax <- merge(x = data_otu09, y =  taxo, by = "OTU_Id")

i_09 <- ggplot(data_otu09_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Cycle)) + geom_point(size = 2) +
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
  labs(x=Dim1otu,y=Dim2otu,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(i_09)

#Figure Zone
j_09 <- ggplot(data_otu09_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Zone)) + geom_point(size = 2) +
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
  labs(x=Dim1otu,y=Dim2otu,color = "Zones", alpha = "Zones", linetype = "Zones")
print(j_09)

#Figure Fraction
k_09 <- ggplot(data_otu09_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Fraction)) + geom_point(size = 2) +
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
  labs(x=Dim1otu,y=Dim2otu,color = "Fraction", alpha = "Fraction", linetype = "Fraction")
print(k_09)

svglite("AFC-Distribution/AFC-OTU-Septembre.svg",width = 12.00,height = 8.00)
Septembre_plot <- plot_grid(i_09, j_09, k_09, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(Septembre_plot)
dev.off()  


    # Novembre_2018 ---------------------------------------------------
data_otu11_tax <- merge(x = data_otu11, y =  taxo, by = "OTU_Id")

i_11 <- ggplot(data_otu11_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Cycle)) + geom_point(size = 2) +
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
  labs(x=Dim1otu,y=Dim2otu,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(i_11)

#Figure Zone
j_11 <- ggplot(data_otu11_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Zone)) + geom_point(size = 2) +
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
  labs(x=Dim1otu,y=Dim2otu,color = "Zones", alpha = "Zones", linetype = "Zones")
print(j_11)

#Figure Fraction
k_11 <- ggplot(data_otu11_tax, aes(y = `Dim 2`, x = `Dim 1`, color = Fraction)) + geom_point(size = 2) +
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
  labs(x=Dim1otu,y=Dim2otu,color = "Fraction", alpha = "Fraction", linetype = "Fraction")
print(k_11)

svglite("AFC-Distribution/AFC-OTU-Novembre.svg",width = 12.00,height = 8.00)
Novembre_plot <- plot_grid(i_11, j_11, k_11, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(Novembre_plot)
dev.off()  


# HIST  --------------------------------------------------------------------
  # OTU ONLY ---------------------------------------------------------------------
    # Cycles -------------------------------------------------------------------
## Jour
onlyJourOTU <- data_otu_tax %>% filter(Cycle == "Jour") %>% select(OTU_Id,TotalJour,Division)
row.names(onlyJourOTU) <- onlyJourOTU$OTU_Id ; onlyJourOTU <- onlyJourOTU %>% select(-OTU_Id)
onlyJourOTU <- onlyJourOTU %>% group_by(Division) %>% summarise_all(sum)
onlyJourOTU$TotalJour <- onlyJourOTU$TotalJour*100/sum(onlyJourOTU$TotalJour)
onlyJourOTU <- onlyJourOTU %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalJour) - 0.5*TotalJour)
onlyJourOTU$label <- paste(round(onlyJourOTU$TotalJour,1), "%", sep = "")
for (i in rownames(onlyJourOTU)) {
  if (onlyJourOTU[i,"label"] == "0%") { onlyJourOTU[i,"label"] <- NA}}
for (i in rownames(onlyJourOTU)) {
  if (is.na(onlyJourOTU[i,"label"]) == FALSE) { onlyJourOTU[i,"label"] <- paste(onlyJourOTU[i,"Division"]," : ",onlyJourOTU[i,"label"], sep = "")}}
onlyJourOTU$Cycle<- rep("Jour", each = nrow(onlyJourOTU))
## Nuit
onlyNuitOTU <- data_otu_tax %>% filter(Cycle == "Nuit") %>% select(OTU_Id,TotalNuit,Division)
row.names(onlyNuitOTU)<-onlyNuitOTU$OTU_Id ; onlyNuitOTU <- onlyNuitOTU %>% select(-OTU_Id)
onlyNuitOTU <- onlyNuitOTU %>% group_by(Division) %>% summarise_all(sum)
onlyNuitOTU$TotalNuit <- onlyNuitOTU$TotalNuit*100/sum(onlyNuitOTU$TotalNuit)
onlyNuitOTU <- onlyNuitOTU %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalNuit) - 0.5*TotalNuit)
onlyNuitOTU$label <- paste(round(onlyNuitOTU$TotalNuit,1), "%", sep = "")
for (i in rownames(onlyNuitOTU)) {
  if (onlyNuitOTU[i,"label"] == "0%") { onlyNuitOTU[i,"label"] <- NA}}
for (i in rownames(onlyNuitOTU)) {
  if (is.na(onlyNuitOTU[i,"label"]) == FALSE) { onlyNuitOTU[i,"label"] <- paste(onlyNuitOTU[i,"Division"]," : ",onlyNuitOTU[i,"label"], sep = "")}}
onlyNuitOTU$Cycle<- rep("Nuit", each = nrow(onlyNuitOTU))
## Cycle
colnames(onlyJourOTU)[2]  <- "value"
onlyJourOTU$Sum <- rep(0, each = nrow(onlyJourOTU))
onlyJourOTU$Count <- rep(0, each = nrow(onlyJourOTU))
for (i in onlyJourOTU$Division) { onlyJourOTU$Count[which(onlyJourOTU$Division == i)] <- nrow(data_otu_tax %>% filter(Cycle == "Jour") %>% filter(Division == i))}
onlyJourOTU$Sum <- sum(onlyJourOTU$Count)
colnames(onlyNuitOTU)[2]  <- "value"
onlyNuitOTU$Sum <- rep(0, each = nrow(onlyNuitOTU))
onlyNuitOTU$Count <- rep(0, each = nrow(onlyNuitOTU))
for (i in onlyNuitOTU$Division) { onlyNuitOTU$Count[which(onlyNuitOTU$Division == i)] <- nrow(data_otu_tax %>% filter(Cycle == "Nuit") %>% filter(Division == i))}
onlyNuitOTU$Sum <- sum(onlyNuitOTU$Count)
onlyCycleOTU <- rbind(onlyJourOTU,onlyNuitOTU)
#Figure
az <- ggplot(onlyCycleOTU, mapping = aes(y= value, x = Cycle, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
  geom_label(aes(y = 104,label = paste(Sum," OTUs",sep ="")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") + scale_fill_manual(values = palette)
legendOTU <- get_legend(az)
az <- az + labs(x="Cycles",y="OTUs (%)") + theme(legend.position = "none")
print(az)


    # Zone -------------------------------------------------------------------
## Oxique
onlyOxiqueOTU <- data_otu_tax %>% filter(Zone == "Oxique") %>% select(OTU_Id,TotalOxique,Division)
row.names(onlyOxiqueOTU)<-onlyOxiqueOTU$OTU_Id ; onlyOxiqueOTU <- onlyOxiqueOTU %>% select(-OTU_Id)
onlyOxiqueOTU <- onlyOxiqueOTU %>% group_by(Division) %>% summarise_all(sum)
onlyOxiqueOTU$TotalOxique <- onlyOxiqueOTU$TotalOxique*100/sum(onlyOxiqueOTU$TotalOxique)
onlyOxiqueOTU <- onlyOxiqueOTU %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalOxique) - 0.5*TotalOxique)
onlyOxiqueOTU$label <- paste(round(onlyOxiqueOTU$TotalOxique,1), "%", sep = "")
for (i in rownames(onlyOxiqueOTU)) {
  if (onlyOxiqueOTU[i,"label"] == "0%") { onlyOxiqueOTU[i,"label"] <- NA}}
for (i in rownames(onlyOxiqueOTU)) {
  if (is.na(onlyOxiqueOTU[i,"label"]) == FALSE) { onlyOxiqueOTU[i,"label"] <- paste(onlyOxiqueOTU[i,"Division"]," : ",onlyOxiqueOTU[i,"label"], sep = "")}}
onlyOxiqueOTU$Zone<- rep("Oxique", each = nrow(onlyOxiqueOTU))
## Anoxique
onlyAnoxiqueOTU <- data_otu_tax %>% filter(Zone == "Anoxique") %>% select(OTU_Id,TotalAnoxique,Division)
row.names(onlyAnoxiqueOTU)<-onlyAnoxiqueOTU$OTU_Id ; onlyAnoxiqueOTU <- onlyAnoxiqueOTU %>% select(-OTU_Id)
onlyAnoxiqueOTU <- onlyAnoxiqueOTU %>% group_by(Division) %>% summarise_all(sum)
onlyAnoxiqueOTU$TotalAnoxique <- onlyAnoxiqueOTU$TotalAnoxique*100/sum(onlyAnoxiqueOTU$TotalAnoxique)
onlyAnoxiqueOTU <- onlyAnoxiqueOTU %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalAnoxique) - 0.5*TotalAnoxique)
onlyAnoxiqueOTU$label <- paste(round(onlyAnoxiqueOTU$TotalAnoxique,1), "%", sep = "")
for (i in rownames(onlyAnoxiqueOTU)) {
  if (onlyAnoxiqueOTU[i,"label"] == "0%") { onlyAnoxiqueOTU[i,"label"] <- NA}}
for (i in rownames(onlyAnoxiqueOTU)) {
  if (is.na(onlyAnoxiqueOTU[i,"label"]) == FALSE) { onlyAnoxiqueOTU[i,"label"] <- paste(onlyAnoxiqueOTU[i,"Division"]," : ",onlyAnoxiqueOTU[i,"label"], sep = "")}}
onlyAnoxiqueOTU$Zone<- rep("Anoxique", each = nrow(onlyAnoxiqueOTU))
## Zone
colnames(onlyOxiqueOTU)[2]  <- "value"
onlyOxiqueOTU$Sum <- rep(0, each = nrow(onlyOxiqueOTU))
onlyOxiqueOTU$Count <- rep(0, each = nrow(onlyOxiqueOTU))
for (i in onlyOxiqueOTU$Division) { onlyOxiqueOTU$Count[which(onlyOxiqueOTU$Division == i)] <- nrow(data_otu_tax %>% filter(Zone == "Oxique") %>% filter(Division == i))}
onlyOxiqueOTU$Sum <- sum(onlyOxiqueOTU$Count)
colnames(onlyAnoxiqueOTU)[2]  <- "value"
onlyAnoxiqueOTU$Sum <- rep(0, each = nrow(onlyAnoxiqueOTU))
onlyAnoxiqueOTU$Count <- rep(0, each = nrow(onlyAnoxiqueOTU))
for (i in onlyAnoxiqueOTU$Division) { onlyAnoxiqueOTU$Count[which(onlyAnoxiqueOTU$Division == i)] <- nrow(data_otu_tax %>% filter(Zone == "Anoxique") %>% filter(Division == i))}
onlyAnoxiqueOTU$Sum <- sum(onlyAnoxiqueOTU$Count)
onlyZoneOTU <- rbind(onlyOxiqueOTU,onlyAnoxiqueOTU)
#Figure
bz <- ggplot(onlyZoneOTU, mapping = aes(y= value, x = Zone, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
  geom_label(aes(y = 104,label = paste(Sum," OTUs",sep ="")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") + scale_fill_manual(values = palette)
bz <- bz + labs(x="Zones",y="OTUs (%)") + theme(legend.position = "none")
print(bz)

    # Fraction -------------------------------------------------------------------
## Petite
onlyPetiteOTU <- data_otu_tax %>% filter(Fraction == "Petite") %>% select(OTU_Id,TotalPetite,Division)
row.names(onlyPetiteOTU)<-onlyPetiteOTU$OTU_Id ; onlyPetiteOTU <- onlyPetiteOTU %>% select(-OTU_Id)
onlyPetiteOTU <- onlyPetiteOTU %>% group_by(Division) %>% summarise_all(sum)
onlyPetiteOTU$TotalPetite <- onlyPetiteOTU$TotalPetite*100/sum(onlyPetiteOTU$TotalPetite)
onlyPetiteOTU <- onlyPetiteOTU %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalPetite) - 0.5*TotalPetite)
onlyPetiteOTU$label <- paste(round(onlyPetiteOTU$TotalPetite,1), "%", sep = "")
for (i in rownames(onlyPetiteOTU)) {
  if (onlyPetiteOTU[i,"label"] == "0%") { onlyPetiteOTU[i,"label"] <- NA}}
for (i in rownames(onlyPetiteOTU)) {
  if (is.na(onlyPetiteOTU[i,"label"]) == FALSE) { onlyPetiteOTU[i,"label"] <- paste(onlyPetiteOTU[i,"Division"]," : ",onlyPetiteOTU[i,"label"], sep = "")}}
onlyPetiteOTU$Fraction<- rep("Petite", each = nrow(onlyPetiteOTU))
## Grande
onlyGrandeOTU <- data_otu_tax %>% filter(Fraction == "Grande") %>% select(OTU_Id,TotalGrande,Division)
row.names(onlyGrandeOTU)<-onlyGrandeOTU$OTU_Id
onlyGrandeOTU <- onlyGrandeOTU %>% select(-OTU_Id)
onlyGrandeOTU <- onlyGrandeOTU %>% group_by(Division) %>% summarise_all(sum)
onlyGrandeOTU$TotalGrande <- onlyGrandeOTU$TotalGrande*100/sum(onlyGrandeOTU$TotalGrande)
onlyGrandeOTU <- onlyGrandeOTU %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalGrande) - 0.5*TotalGrande)
onlyGrandeOTU$label <- paste(round(onlyGrandeOTU$TotalGrande,1), "%", sep = "")
for (i in rownames(onlyGrandeOTU)) {
  if (onlyGrandeOTU[i,"label"] == "0%") { onlyGrandeOTU[i,"label"] <- NA}}
for (i in rownames(onlyGrandeOTU)) {
  if (is.na(onlyGrandeOTU[i,"label"]) == FALSE) { onlyGrandeOTU[i,"label"] <- paste(onlyGrandeOTU[i,"Division"]," : ",onlyGrandeOTU[i,"label"], sep = "")}}
onlyGrandeOTU$Fraction<- rep("Grande", each = nrow(onlyGrandeOTU))
## Fraction
colnames(onlyPetiteOTU)[2]  <- "value"
onlyPetiteOTU$Sum <- rep(0, each = nrow(onlyPetiteOTU))
onlyPetiteOTU$Count <- rep(0, each = nrow(onlyPetiteOTU))
for (i in onlyPetiteOTU$Division) { onlyPetiteOTU$Count[which(onlyPetiteOTU$Division == i)] <- nrow(data_otu_tax %>% filter(Fraction == "Petite") %>% filter(Division == i))}
onlyPetiteOTU$Sum <- sum(onlyPetiteOTU$Count)
colnames(onlyGrandeOTU)[2]  <- "value"
onlyGrandeOTU$Sum <- rep(0, each = nrow(onlyGrandeOTU))
onlyGrandeOTU$Count <- rep(0, each = nrow(onlyGrandeOTU))
for (i in onlyGrandeOTU$Division) { onlyGrandeOTU$Count[which(onlyGrandeOTU$Division == i)] <- nrow(data_otu_tax %>% filter(Fraction == "Grande") %>% filter(Division == i))}
onlyGrandeOTU$Sum <- sum(onlyGrandeOTU$Count)
onlyFractionOTU <- rbind(onlyPetiteOTU,onlyGrandeOTU)
#Figure
cz <- ggplot(onlyFractionOTU, mapping = aes(y= value, x = Fraction, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
  geom_label(aes(y = 104,label = paste(Sum," OTUs",sep ="")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
  scale_fill_manual(values = palette)
cz <- cz + labs(x="Fractions",y="OTUs (%)") + theme(legend.position = "none")
print(cz) 
    # Coplot -------------------------------------------------------------------
svglite("HistOnly/OTU-only.svg",width = 10.00,height = 8.00)
b_plot <- plot_grid(az,bz,cz,legendOTU,  ncol = 4, nrow = 1, rel_widths = c(3,3,3,3),rel_heights = c(3))
print(b_plot)
dev.off()  
    # Coplot X AFC -------------------------------------------------------------------
f <- ggplot() + theme_void()
print(f)
#AFC = séquence
svglite("Biplot/OTU-onlyX-AFCseq.svg",width = 12.00,height = 11.00)
b_plot <- plot_grid(aCycle,az,f,bZone,bz,legendOTU,cFraction,cz,f, ncol = 3, nrow = 3, rel_widths = c(3,1.5,1),rel_heights = c(3,3,3))
print(b_plot)
dev.off()  
#AFC = OTU
svglite("Biplot/OTU-onlyX-AFCotu.svg",width = 12.00,height = 11.00)
c_plot <- plot_grid(iCycle,az,f,jZone,bz,legendOTU,kFraction,cz,f, ncol = 3, nrow = 3, rel_widths = c(3,1.5,1),rel_heights = c(3,3,3))
print(c_plot)
dev.off()  



  # OTU TOTAL---------------------------------------------------------------------
    # Cycle -------------------------------------------------------------------
## Jour
totalJourOTU <- data_otu_tax %>% select(OTU_Id,TotalJour,Division)
row.names(totalJourOTU)<-totalJourOTU$OTU_Id ; totalJourOTU <- totalJourOTU %>% select(-OTU_Id)
totalJourOTU <- totalJourOTU %>% group_by(Division) %>% summarise_all(sum)
totalJourOTU$TotalJour <- totalJourOTU$TotalJour*100/sum(totalJourOTU$TotalJour)
totalJourOTU <- totalJourOTU %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalJour) - 0.5*TotalJour)
totalJourOTU$label <- paste(round(totalJourOTU$TotalJour,1), "%", sep = "")
for (i in rownames(totalJourOTU)) {
  if (totalJourOTU[i,"label"] == "0%") { totalJourOTU[i,"label"] <- NA}}
for (i in rownames(totalJourOTU)) {
  if (is.na(totalJourOTU[i,"label"]) == FALSE) { totalJourOTU[i,"label"] <- paste(totalJourOTU[i,"Division"]," : ",totalJourOTU[i,"label"], sep = "")}}
totalJourOTU$Cycle<- rep("Jour", each = nrow(totalJourOTU))
##Nuit
totalNuitOTU <- data_otu_tax %>% select(OTU_Id,TotalNuit,Division)
row.names(totalNuitOTU)<-totalNuitOTU$OTU_Id ; totalNuitOTU <- totalNuitOTU %>% select(-OTU_Id)
totalNuitOTU <- totalNuitOTU %>% group_by(Division) %>% summarise_all(sum)
totalNuitOTU$TotalNuit <- totalNuitOTU$TotalNuit*100/sum(totalNuitOTU$TotalNuit)
totalNuitOTU <- totalNuitOTU %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalNuit) - 0.5*TotalNuit)
totalNuitOTU$label <- paste(round(totalNuitOTU$TotalNuit,1), "%", sep = "")
for (i in rownames(totalNuitOTU)) {
  if (totalNuitOTU[i,"label"] == "0%") { totalNuitOTU[i,"label"] <- NA}}
for (i in rownames(totalNuitOTU)) {
  if (is.na(totalNuitOTU[i,"label"]) == FALSE) { totalNuitOTU[i,"label"] <- paste(totalNuitOTU[i,"Division"]," : ",totalNuitOTU[i,"label"], sep = "")}}
totalNuitOTU$Cycle<- rep("Nuit", each = nrow(totalNuitOTU))
##Cycle
colnames(totalJourOTU)[2]  <- "value"
totalJourOTU$Sum <- rep(0, each = nrow(totalJourOTU))
totalJourOTU$Count <- rep(0, each = nrow(totalJourOTU))
for (i in totalJourOTU$Division) { totalJourOTU$Count[which(totalJourOTU$Division == i)] <- sum(data_otu_tax  %>% filter(Division == i) %>% select(TotalJour))}
totalJourOTU$Sum <- sum(totalJourOTU$Count)
colnames(totalNuitOTU)[2]  <- "value"
totalNuitOTU$Sum <- rep(0, each = nrow(totalNuitOTU))
totalNuitOTU$Count <- rep(0, each = nrow(totalNuitOTU))
for (i in totalNuitOTU$Division) { totalNuitOTU$Count[which(totalNuitOTU$Division == i)] <- sum(data_otu_tax  %>% filter(Division == i) %>% select(TotalNuit))}
totalNuitOTU$Sum <- sum(totalNuitOTU$Count)
totalCycleOTU <- rbind(totalJourOTU,totalNuitOTU)

#Figure
ay <- ggplot(totalCycleOTU, mapping = aes(y= value, x = Cycle, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + scale_fill_manual(values = palette) + 
  geom_label(aes(y = 104,label = paste(Sum," OTUs",sep ="")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
legendOTU <- get_legend(ay)
ay <- ay + labs(x="Cycles",y="OTUs (%)") + theme(legend.position = "none")
print(ay)

    # Zone -------------------------------------------------------------------
## Oxique
totalOxiqueOTU <- data_otu_tax %>% select(OTU_Id,TotalOxique,Division)
row.names(totalOxiqueOTU)<-totalOxiqueOTU$OTU_Id ; totalOxiqueOTU <- totalOxiqueOTU %>% select(-OTU_Id)
totalOxiqueOTU <- totalOxiqueOTU %>% group_by(Division) %>% summarise_all(sum)
totalOxiqueOTU$TotalOxique <- totalOxiqueOTU$TotalOxique*100/sum(totalOxiqueOTU$TotalOxique)
totalOxiqueOTU <- totalOxiqueOTU %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalOxique) - 0.5*TotalOxique)
totalOxiqueOTU$label <- paste(round(totalOxiqueOTU$TotalOxique,1), "%", sep = "")
for (i in rownames(totalOxiqueOTU)) {
  if (totalOxiqueOTU[i,"label"] == "0%") { totalOxiqueOTU[i,"label"] <- NA}}
for (i in rownames(totalOxiqueOTU)) {
  if (is.na(totalOxiqueOTU[i,"label"]) == FALSE) { totalOxiqueOTU[i,"label"] <- paste(totalOxiqueOTU[i,"Division"]," : ",totalOxiqueOTU[i,"label"], sep = "")}}
totalOxiqueOTU$Zone<- rep("Oxique", each = nrow(totalOxiqueOTU))
## Anoxique
totalAnoxiqueOTU <- data_otu_tax %>% select(OTU_Id,TotalAnoxique,Division)
row.names(totalAnoxiqueOTU)<-totalAnoxiqueOTU$OTU_Id ; totalAnoxiqueOTU <- totalAnoxiqueOTU %>% select(-OTU_Id)
totalAnoxiqueOTU <- totalAnoxiqueOTU %>% group_by(Division) %>% summarise_all(sum)
totalAnoxiqueOTU$TotalAnoxique <- totalAnoxiqueOTU$TotalAnoxique*100/sum(totalAnoxiqueOTU$TotalAnoxique)
totalAnoxiqueOTU <- totalAnoxiqueOTU %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalAnoxique) - 0.5*TotalAnoxique)
totalAnoxiqueOTU$label <- paste(round(totalAnoxiqueOTU$TotalAnoxique,1), "%", sep = "")
for (i in rownames(totalAnoxiqueOTU)) {
  if (totalAnoxiqueOTU[i,"label"] == "0%") { totalAnoxiqueOTU[i,"label"] <- NA}}
for (i in rownames(totalAnoxiqueOTU)) {
  if (is.na(totalAnoxiqueOTU[i,"label"]) == FALSE) { totalAnoxiqueOTU[i,"label"] <- paste(totalAnoxiqueOTU[i,"Division"]," : ",totalAnoxiqueOTU[i,"label"], sep = "")}}
totalAnoxiqueOTU$Zone<- rep("Anoxique", each = nrow(totalAnoxiqueOTU))
## Zone
colnames(totalOxiqueOTU)[2]  <- "value"
totalOxiqueOTU$Sum <- rep(0, each = nrow(totalOxiqueOTU))
totalOxiqueOTU$Count <- rep(0, each = nrow(totalOxiqueOTU))
for (i in totalOxiqueOTU$Division) { totalOxiqueOTU$Count[which(totalOxiqueOTU$Division == i)] <- sum(data_otu_tax  %>% filter(Division == i) %>% select(TotalOxique))}
totalOxiqueOTU$Sum <- sum(totalOxiqueOTU$Count)
colnames(totalAnoxiqueOTU)[2]  <- "value"
totalAnoxiqueOTU$Sum <- rep(0, each = nrow(totalAnoxiqueOTU))
totalAnoxiqueOTU$Count <- rep(0, each = nrow(totalAnoxiqueOTU))
for (i in totalAnoxiqueOTU$Division) { totalAnoxiqueOTU$Count[which(totalAnoxiqueOTU$Division == i)] <- sum(data_otu_tax  %>% filter(Division == i) %>% select(TotalAnoxique))}
totalAnoxiqueOTU$Sum <- sum(totalAnoxiqueOTU$Count)
totalZoneOTU <- rbind(totalOxiqueOTU,totalAnoxiqueOTU)
#Figure
by <- ggplot(totalZoneOTU, mapping = aes(y= value, x = Zone, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + scale_fill_manual(values = palette) + 
  geom_label(aes(y = 104,label = paste(Sum," OTUs",sep ="")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
by <- by + labs(x="Zones",y="OTUs (%)") + theme(legend.position = "none")
print(by)

    # Fraction -------------------------------------------------------------------
## Petite
totalPetiteOTU <- data_otu_tax %>% select(OTU_Id,TotalPetite,Division)
row.names(totalPetiteOTU)<-totalPetiteOTU$OTU_Id ; totalPetiteOTU <- totalPetiteOTU %>% select(-OTU_Id)
totalPetiteOTU <- totalPetiteOTU %>% group_by(Division) %>% summarise_all(sum)
totalPetiteOTU$TotalPetite <- totalPetiteOTU$TotalPetite*100/sum(totalPetiteOTU$TotalPetite)
totalPetiteOTU <- totalPetiteOTU %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalPetite) - 0.5*TotalPetite)
totalPetiteOTU$label <- paste(round(totalPetiteOTU$TotalPetite,1), "%", sep = "")
for (i in rownames(totalPetiteOTU)) {
  if (totalPetiteOTU[i,"label"] == "0%") { totalPetiteOTU[i,"label"] <- NA}}
for (i in rownames(totalPetiteOTU)) {
  if (is.na(totalPetiteOTU[i,"label"]) == FALSE) { totalPetiteOTU[i,"label"] <- paste(totalPetiteOTU[i,"Division"]," : ",totalPetiteOTU[i,"label"], sep = "")}}
totalPetiteOTU$Fraction<- rep("Petite", each = nrow(totalPetiteOTU))
## Grande
totalGrandeOTU <- data_otu_tax %>% select(OTU_Id,TotalGrande,Division)
row.names(totalGrandeOTU)<-totalGrandeOTU$OTU_Id ; totalGrandeOTU <- totalGrandeOTU %>% select(-OTU_Id)
totalGrandeOTU <- totalGrandeOTU %>% group_by(Division) %>% summarise_all(sum)
totalGrandeOTU$TotalGrande <- totalGrandeOTU$TotalGrande*100/sum(totalGrandeOTU$TotalGrande)
totalGrandeOTU <- totalGrandeOTU %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(TotalGrande) - 0.5*TotalGrande)
totalGrandeOTU$label <- paste(round(totalGrandeOTU$TotalGrande,1), "%", sep = "")
for (i in rownames(totalGrandeOTU)) {
  if (totalGrandeOTU[i,"label"] == "0%") { totalGrandeOTU[i,"label"] <- NA}}
for (i in rownames(totalGrandeOTU)) {
  if (is.na(totalGrandeOTU[i,"label"]) == FALSE) { totalGrandeOTU[i,"label"] <- paste(totalGrandeOTU[i,"Division"]," : ",totalGrandeOTU[i,"label"], sep = "")}}
totalGrandeOTU$Fraction<- rep("Grande", each = nrow(totalGrandeOTU))
## Fraction
colnames(totalPetiteOTU)[2]  <- "value"
totalPetiteOTU$Sum <- rep(0, each = nrow(totalPetiteOTU))
totalPetiteOTU$Count <- rep(0, each = nrow(totalPetiteOTU))
for (i in totalPetiteOTU$Division) { totalPetiteOTU$Count[which(totalPetiteOTU$Division == i)] <- sum(data_otu_tax  %>% filter(Division == i) %>% select(TotalPetite))}
totalPetiteOTU$Sum <- sum(totalPetiteOTU$Count)
colnames(totalGrandeOTU)[2]  <- "value"
totalGrandeOTU$Sum <- rep(0, each = nrow(totalGrandeOTU))
totalGrandeOTU$Count <- rep(0, each = nrow(totalGrandeOTU))
for (i in totalGrandeOTU$Division) { totalGrandeOTU$Count[which(totalGrandeOTU$Division == i)] <- sum(data_otu_tax  %>% filter(Division == i) %>% select(TotalGrande))}
totalGrandeOTU$Sum <- sum(totalGrandeOTU$Count)
totalFractionOTU <- rbind(totalPetiteOTU,totalGrandeOTU)
#Figure
cy <- ggplot(totalFractionOTU, mapping = aes(y= value, x = Fraction, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + scale_fill_manual(values = palette)+ 
  geom_label(aes(y = 104,label = paste(Sum," OTUs",sep ="")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
cy <- cy + labs(x="Fractions",y="OTUs (%)") + theme(legend.position = "none")
print(cy) 
    # Date -------------------------------------------------------------------
## Avril
total04OTU <- data_otu_tax %>% select(OTU_Id,Total04,Division)
row.names(total04OTU)<-total04OTU$OTU_Id ; total04OTU <- total04OTU %>% select(-OTU_Id)
total04OTU <- total04OTU %>% group_by(Division) %>% summarise_all(sum)
total04OTU$Total04 <- total04OTU$Total04*100/sum(total04OTU$Total04)
total04OTU <- total04OTU %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(Total04) - 0.5*Total04)
total04OTU$label <- paste(round(total04OTU$Total04,1), "%", sep = "")
for (i in rownames(total04OTU)) {
  if (total04OTU[i,"label"] == "0%") { total04OTU[i,"label"] <- NA}}
for (i in rownames(total04OTU)) {
  if (is.na(total04OTU[i,"label"]) == FALSE) { total04OTU[i,"label"] <- paste(total04OTU[i,"Division"]," : ",total04OTU[i,"label"], sep = "")}}
total04OTU$Date<- rep("04", each = nrow(total04OTU))
## Juin
total06OTU <- data_otu_tax %>% select(OTU_Id,Total06,Division)
row.names(total06OTU)<-total06OTU$OTU_Id ; total06OTU <- total06OTU %>% select(-OTU_Id)
total06OTU <- total06OTU %>% group_by(Division) %>% summarise_all(sum)
total06OTU$Total06 <- total06OTU$Total06*100/sum(total06OTU$Total06)
total06OTU <- total06OTU %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(Total06) - 0.5*Total06)
total06OTU$label <- paste(round(total06OTU$Total06,1), "%", sep = "")
for (i in rownames(total06OTU)) {
  if (total06OTU[i,"label"] == "0%") { total06OTU[i,"label"] <- NA}}
for (i in rownames(total06OTU)) {
  if (is.na(total06OTU[i,"label"]) == FALSE) { total06OTU[i,"label"] <- paste(total06OTU[i,"Division"]," : ",total06OTU[i,"label"], sep = "")}}
total06OTU$Date<- rep("06", each = nrow(total06OTU))
## Septembre
total09OTU <- data_otu_tax %>% select(OTU_Id,Total09,Division)
row.names(total09OTU)<-total09OTU$OTU_Id ; total09OTU <- total09OTU %>% select(-OTU_Id)
total09OTU <- total09OTU %>% group_by(Division) %>% summarise_all(sum)
total09OTU$Total09 <- total09OTU$Total09*100/sum(total09OTU$Total09)
total09OTU <- total09OTU %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(Total09) - 0.5*Total09)
total09OTU$label <- paste(round(total09OTU$Total09,1), "%", sep = "")
for (i in rownames(total09OTU)) {
  if (total09OTU[i,"label"] == "0%") { total09OTU[i,"label"] <- NA}}
for (i in rownames(total09OTU)) {
  if (is.na(total09OTU[i,"label"]) == FALSE) { total09OTU[i,"label"] <- paste(total09OTU[i,"Division"]," : ",total09OTU[i,"label"], sep = "")}}
total09OTU$Date<- rep("09", each = nrow(total09OTU))
## Octobre
total11OTU <- data_otu_tax %>% select(OTU_Id,Total11,Division)
row.names(total11OTU)<-total11OTU$OTU_Id ; total11OTU <- total11OTU %>% select(-OTU_Id)
total11OTU <- total11OTU %>% group_by(Division) %>% summarise_all(sum)
total11OTU$Total11 <- total11OTU$Total11*100/sum(total11OTU$Total11)
total11OTU <- total11OTU %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(Total11) - 0.5*Total11)
total11OTU$label <- paste(round(total11OTU$Total11,1), "%", sep = "")
for (i in rownames(total11OTU)) {
  if (total11OTU[i,"label"] == "0%") { total11OTU[i,"label"] <- NA}}
for (i in rownames(total11OTU)) {
  if (is.na(total11OTU[i,"label"]) == FALSE) { total11OTU[i,"label"] <- paste(total11OTU[i,"Division"]," : ",total11OTU[i,"label"], sep = "")}}
total11OTU$Date<- rep("11", each = nrow(total11OTU))
## Date
colnames(total04OTU)[2]  <- "value"
total04OTU$Sum <- rep(0, each = nrow(total04OTU))
total04OTU$Count <- rep(0, each = nrow(total04OTU))
for (i in total04OTU$Division) { total04OTU$Count[which(total04OTU$Division == i)] <- sum(data_otu_tax  %>% filter(Division == i) %>% select(Total04))}
total04OTU$Sum <- sum(total04OTU$Count)
colnames(total06OTU)[2]  <- "value"
total06OTU$Sum <- rep(0, each = nrow(total06OTU))
total06OTU$Count <- rep(0, each = nrow(total06OTU))
for (i in total06OTU$Division) { total06OTU$Count[which(total06OTU$Division == i)] <- sum(data_otu_tax  %>% filter(Division == i) %>% select(Total06))}
total06OTU$Sum <- sum(total06OTU$Count)
colnames(total09OTU)[2]  <- "value"
total09OTU$Sum <- rep(0, each = nrow(total09OTU))
total09OTU$Count <- rep(0, each = nrow(total09OTU))
for (i in total09OTU$Division) { total09OTU$Count[which(total09OTU$Division == i)] <- sum(data_otu_tax  %>% filter(Division == i) %>% select(Total09))}
total09OTU$Sum <- sum(total09OTU$Count)
colnames(total11OTU)[2]  <- "value"
total11OTU$Sum <- rep(0, each = nrow(total11OTU))
total11OTU$Count <- rep(0, each = nrow(total11OTU))
for (i in total11OTU$Division) { total11OTU$Count[which(total11OTU$Division == i)] <- sum(data_otu_tax  %>% filter(Division == i) %>% select(Total11))}
total11OTU$Sum <- sum(total11OTU$Count)
totalDateOTU <- rbind(total04OTU,total06OTU,total09OTU,total11OTU)
#Figure
dy <- ggplot(totalDateOTU, mapping = aes(y= value, x = Date, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + scale_fill_manual(values = palette)+ 
  geom_label(aes(y = 104,label = paste(Sum," OTUs",sep ="")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
dy <- dy + theme(legend.position = "none") + labs(x="Dates",y="OTUs (%)")
print(dy) 
    # Coplot -------------------------------------------------------------------
svglite("HistOnly/OTU-Total.svg",width = 14.00,height = 6.00)
b_plot <- plot_grid(ay,by,cy,dy,legendOTU, ncol = 5, nrow = 1, rel_widths = c(3,3,3,5,3),rel_heights = c(3))
print(b_plot)
dev.off()
  # Séquence ONLY---------------------------------------------------------------------
    # Cycle -------------------------------------------------------------------
## Jour
onlyJourSequence <- data_seq_tax %>% filter(Cycle == "Jour") %>% select(OTU_Id,TotalJour,Division)
row.names(onlyJourSequence)<-onlyJourSequence$OTU_Id ; onlyJourSequence <- onlyJourSequence %>% select(-OTU_Id)
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
onlyNuitSequence <- data_seq_tax %>% filter(Cycle == "Nuit") %>% select(OTU_Id,TotalNuit,Division)
row.names(onlyNuitSequence)<-onlyNuitSequence$OTU_Id ; onlyNuitSequence <- onlyNuitSequence %>% select(-OTU_Id)
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
legendSequence <- get_legend(iz)
iz <- iz + labs(x="Cycles",y="Séquences (%)") + theme(legend.position = "none")
print(iz)

    # Zone -------------------------------------------------------------------
## Oxique
onlyOxiqueSequence <- data_seq_tax %>% filter(Zone == "Oxique") %>% select(OTU_Id,TotalOxique,Division)
row.names(onlyOxiqueSequence)<-onlyOxiqueSequence$OTU_Id ; onlyOxiqueSequence <- onlyOxiqueSequence %>% select(-OTU_Id)
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
onlyAnoxiqueSequence <- data_seq_tax %>% filter(Zone == "Anoxique") %>% select(OTU_Id,TotalAnoxique,Division)
row.names(onlyAnoxiqueSequence)<-onlyAnoxiqueSequence$OTU_Id ; onlyAnoxiqueSequence <- onlyAnoxiqueSequence %>% select(-OTU_Id)
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
jz <- jz + theme(legend.position = "none") + labs(x="Zones",y="Séquences (%)")
print(jz)    
    # Fraction -------------------------------------------------------------------
## Petite
onlyPetiteSequence <- data_seq_tax %>% filter(Fraction == "Petite") %>% select(OTU_Id,TotalPetite,Division)
row.names(onlyPetiteSequence)<-onlyPetiteSequence$OTU_Id ; onlyPetiteSequence <- onlyPetiteSequence %>% select(-OTU_Id)
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
onlyGrandeSequence <- data_seq_tax %>% filter(Fraction == "Grande") %>% select(OTU_Id,TotalGrande,Division)
row.names(onlyGrandeSequence)<-onlyGrandeSequence$OTU_Id ; onlyGrandeSequence <- onlyGrandeSequence %>% select(-OTU_Id)
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
kz <- kz + labs(x="Fractions",y="Séquences (%)") + theme(legend.position = "none")
print(kz) 
    # Coplot -------------------------------------------------------------------
svglite("HistOnly/Sequence-only.svg",width = 10.00,height = 8.00)
b_plot <- plot_grid(iz,jz,kz,legendSequence,  ncol = 4, nrow = 1, rel_widths = c(3,3,3,3),rel_heights = c(3))
print(b_plot)
dev.off()  

    # Coplot X AFC -------------------------------------------------------------------
f <- ggplot() + theme_void()
print(f)
#AFC = séquence
svglite("Biplot/Sequence-onlyX-AFCseq.svg",width = 12.00,height = 11.00)
b_plot <- plot_grid(aCycle,iz,f,bZone,jz,legendSequence,cFraction,kz,f, ncol = 3, nrow = 3, rel_widths = c(3,1.5,1),rel_heights = c(3,3,3))
print(b_plot)
dev.off()  
#AFC = OTU
svglite("Biplot/Sequence-onlyX-AFCotu.svg",width = 12.00,height = 11.00)
c_plot <- plot_grid(iCycle,iz,f,jZone,jz,legendSequence,kFraction,kz,f, ncol = 3, nrow = 3, rel_widths = c(3,1.5,1),rel_heights = c(3,3,3))
print(c_plot)
dev.off()


  # Séquence TOTAL---------------------------------------------------------------------
    # Cycle -------------------------------------------------------------------
## Jour
totalJourSequence <- data_seq_tax %>% select(OTU_Id,TotalJour,Division)
row.names(totalJourSequence)<-totalJourSequence$OTU_Id ; totalJourSequence <- totalJourSequence %>% select(-OTU_Id)
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
totalNuitSequence <- data_seq_tax %>% select(OTU_Id,TotalNuit,Division)
row.names(totalNuitSequence)<-totalNuitSequence$OTU_Id ; totalNuitSequence <- totalNuitSequence %>% select(-OTU_Id)
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
  geom_label(aes(y = 106,label = paste(Sum,"séquences",percent,sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
legendSequence <- get_legend(iy)
iy <- iy + labs(x="Cycles",y="Séquences (%)") + theme(legend.position = "none")
print(iy)

    # Zone -------------------------------------------------------------------
## Oxique
totalOxiqueSequence <- data_seq_tax %>% select(OTU_Id,TotalOxique,Division)
row.names(totalOxiqueSequence)<-totalOxiqueSequence$OTU_Id ; totalOxiqueSequence <- totalOxiqueSequence %>% select(-OTU_Id)
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
totalAnoxiqueSequence <- data_seq_tax %>% select(OTU_Id,TotalAnoxique,Division)
row.names(totalAnoxiqueSequence)<-totalAnoxiqueSequence$OTU_Id ; totalAnoxiqueSequence <- totalAnoxiqueSequence %>% select(-OTU_Id)
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
  geom_label(aes(y = 106,label = paste(Sum,"séquences",percent,sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
jy <- jy + labs(x="Zones",y="Séquences (%)") + theme(legend.position = "none")
print(jy)    
    # Fraction -------------------------------------------------------------------
## Petite
totalPetiteSequence <- data_seq_tax %>% select(OTU_Id,TotalPetite,Division)
row.names(totalPetiteSequence)<-totalPetiteSequence$OTU_Id ; totalPetiteSequence <- totalPetiteSequence %>% select(-OTU_Id)
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
totalGrandeSequence <- data_seq_tax %>% select(OTU_Id,TotalGrande,Division)
row.names(totalGrandeSequence)<-totalGrandeSequence$OTU_Id ; totalGrandeSequence <- totalGrandeSequence %>% select(-OTU_Id)
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
  geom_label(aes(y = 106,label = paste(Sum,"séquences",percent,sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
ky <- ky + labs(x="Fractions",y="Séquences (%)") + theme(legend.position = "none")
print(ky) 
    # Date -------------------------------------------------------------------
## Avril
total04Sequence <- data_seq_tax %>% select(OTU_Id,Total04,Division)
row.names(total04Sequence)<-total04Sequence$OTU_Id ; total04Sequence <- total04Sequence %>% select(-OTU_Id)
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
total06Sequence <- data_seq_tax %>% select(OTU_Id,Total06,Division)
row.names(total06Sequence)<-total06Sequence$OTU_Id ; total06Sequence <- total06Sequence %>% select(-OTU_Id)
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
total09Sequence <- data_seq_tax %>% select(OTU_Id,Total09,Division)
row.names(total09Sequence)<-total09Sequence$OTU_Id ; total09Sequence <- total09Sequence %>% select(-OTU_Id)
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
total11Sequence <- data_seq_tax %>% select(OTU_Id,Total11,Division)
row.names(total11Sequence)<-total11Sequence$OTU_Id ; total11Sequence <- total11Sequence %>% select(-OTU_Id)
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
  geom_label(aes(y = 106,label = paste(Sum,"séquences",percent,sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
ly <- ly + labs(x="Dates",y="Séquences (%)") + theme(legend.position = "none")
print(ly) 
    # Coplot -------------------------------------------------------------------
svglite("HistOnly/Sequence-Total.svg",width = 14.00,height = 6.00)
b_plot <- plot_grid(iy,jy,ky,ly, legendSequence, ncol = 5, nrow = 1, rel_widths = c(3,3,3,5,3),rel_heights = c(3))
print(b_plot)
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
svglite("Composition/Polar-Total-seq.svg",width = 4.00,height = 4.00)
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

  # OTU ----------------------------------------------------------------
##Total column
Polar_otu <- data_otu_tax
Polar_otu$TotalAnnée <- 0
for (i in rownames(Polar_otu)) { if (Polar_otu[i,"Total04"] + Polar_otu[i,"Total06"] + Polar_otu[i,"Total09"] + Polar_otu[i,"Total11"] > 0) {Polar_otu[i,"TotalAnnée"] <- 1}}
Polar_otu <- Polar_otu %>% select(TotalAnnée,Division)
## Total Figure
Polar_otu <- Polar_otu %>% group_by(Division) %>% summarise_all(sum)
for (i in rownames(Polar_otu)) { Polar_otu[i,"Total"] <- (Polar_otu[i,"TotalAnnée"] * 100) / sum(Polar_otu$TotalAnnée)}
##Label
Polar_otu <- Polar_otu %>%
  arrange(desc(Division)) %>%
  mutate(lab.ypos = cumsum(Total) - 0.5*Total)
Polar_otu$label <- paste(round(Polar_otu$Total,1), "%", sep = "")
for (i in rownames(Polar_otu)) {
  if (Polar_otu[i,"label"] == "0%") { Polar_otu[i,"label"] <- NA}}
for (i in rownames(Polar_otu)) {
  if (is.na(Polar_otu[i,"label"]) == FALSE) { Polar_otu[i,"label"] <- paste(Polar_otu[i,"Division"]," : ",Polar_otu[i,"label"], sep = "")}}
##Figure
svglite("Composition/Polar-Total-otu.svg",width = 4.00,height = 4.00)
bx <- ggplot(Polar_otu, mapping = aes(y= Total, x = 2, fill = Division), Rowv = NA, col = colMain, scale = "column") +
  geom_bar(stat="identity", color = "white", width = 1) + coord_polar("y") + 
  geom_label_repel(aes(y = lab.ypos,label = label),color = "white",size = 3,segment.color = "black",show.legend = FALSE, nudge_x = 0.5) +
  scale_y_continuous(limits=c(0,sum(Polar_otu %>% select(Total))))+
  xlim(1,2.5) +
  theme_unique_darkbis() + 
  #facet_wrap( ~ variable , nrow = 4) +
  labs(y = "OTUs",x="") + scale_fill_manual(values = rev(palette))
print(bx)
dev.off()

# Table OTUs majoritaires -------------------------------------------------
# Table ONLY---------------------------------------------------------------------
  # Cycle -------------------------------------------------------------------
##Jour
onlyJourTable <- data_seq_tax %>% filter(Cycle == "Jour") %>% select(OTU_Id,TotalJour,Division,Cycle)
onlyJourTable <- onlyJourTable %>% filter(TotalJour > 0.01*sum(onlyJourTable$TotalJour))
colnames(onlyJourTable)[2]  <- "value"
##Nuit
onlyNuitTable <- data_seq_tax %>% filter(Cycle == "Nuit") %>% select(OTU_Id,TotalNuit,Division,Cycle)
onlyNuitTable <- onlyNuitTable %>% filter(TotalNuit > 0.01*sum(onlyNuitTable$TotalNuit))
colnames(onlyNuitTable)[2]  <- "value"
##Bind
onlyCycleTable <- rbind(onlyJourTable,onlyNuitTable)
onlyCycleTable <- merge(onlyCycleTable,taxo, by = "OTU_Id")
##Table
write.table(onlyCycleTable, file = "TableOnly/Only_Cycles.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(onlyCycleOTU, file = "TableOnly/Only_Cycles_OTU.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(onlyCycleSequence, file = "TableOnly/Only_Cycles_Sequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  # Zone -------------------------------------------------------------------
##Oxique
onlyOxiqueTable <- data_seq_tax %>% filter(Zone == "Oxique") %>% select(OTU_Id,TotalOxique,Division,Zone)
onlyOxiqueTable <- onlyOxiqueTable %>% filter(TotalOxique > 0.01*sum(onlyOxiqueTable$TotalOxique))
colnames(onlyOxiqueTable)[2]  <- "value"
##Anoxique
onlyAnoxiqueTable <- data_seq_tax %>% filter(Zone == "Anoxique") %>% select(OTU_Id,TotalAnoxique,Division,Zone)
onlyAnoxiqueTable <- onlyAnoxiqueTable %>% filter(TotalAnoxique > 0.01*sum(onlyAnoxiqueTable$TotalAnoxique))
colnames(onlyAnoxiqueTable)[2]  <- "value"
##bind
onlyZoneTable <- rbind(onlyOxiqueTable,onlyAnoxiqueTable)
onlyZoneTable <- merge(onlyZoneTable,taxo, by = "OTU_Id")
##Table
write.table(onlyZoneTable, file = "TableOnly/Only_Zone.txt", sep = "\t", col.names = TRUE, row.names = FALSE,quote = FALSE)
write.table(onlyZoneOTU, file = "TableOnly/Only_Zone_OTU.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(onlyZoneSequence, file = "TableOnly/Only_Zone_Sequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  # Fraction -------------------------------------------------------------------
##Petite
onlyPetiteTable <- data_seq_tax %>% filter(Fraction == "Petite") %>% select(OTU_Id,TotalPetite,Division,Fraction)
onlyPetiteTable <- onlyPetiteTable %>% filter(TotalPetite > 0.01*sum(onlyPetiteTable$TotalPetite))
colnames(onlyPetiteTable)[2]  <- "value"
##Grande
onlyGrandeTable <- data_seq_tax %>% filter(Fraction == "Grande") %>% select(OTU_Id,TotalGrande,Division,Fraction)
onlyGrandeTable <- onlyGrandeTable %>% filter(TotalGrande > 0.01*sum(onlyGrandeTable$TotalGrande))
colnames(onlyGrandeTable)[2]  <- "value"
##bind
onlyFractionTable <- rbind(onlyPetiteTable,onlyGrandeTable)
onlyFractionTable <- merge(onlyFractionTable,taxo, by = "OTU_Id")
##Table
write.table(onlyFractionTable, file = "TableOnly/Only_Fraction.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(onlyFractionOTU, file = "TableOnly/Only_Fraction_OTU.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(onlyFractionSequence, file = "TableOnly/Only_Fraction_Sequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Table Total---------------------------------------------------------------------
  # Cycle Total -------------------------------------------------------------------
##Jour
totalJourTable <- data_seq_tax %>% select(OTU_Id,TotalJour,Division,Cycle)
totalJourTable$Cycle <- "Jour"
totalJourTable <- totalJourTable %>% filter(TotalJour > 0.05*sum(totalJourTable$TotalJour))
colnames(totalJourTable)[2]  <- "value"
##Nuit
totalNuitTable <- data_seq_tax %>% select(OTU_Id,TotalNuit,Division,Cycle)
totalNuitTable$Cycle <- "Nuit"
totalNuitTable <- totalNuitTable %>% filter(TotalNuit > 0.05*sum(totalNuitTable$TotalNuit))
colnames(totalNuitTable)[2]  <- "value"
##Bind
totalCycleTable <- rbind(totalJourTable,totalNuitTable)
totalCycleTable <- merge(totalCycleTable,taxo, by = "OTU_Id")
##Table
write.table(totalCycleTable, file = "TableOnly/Total_Cycles.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(totalCycleOTU, file = "TableOnly/Total_Cycles_OTU.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(totalCycleSequence, file = "TableOnly/Total_Cycles_Sequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  # Zone Total -------------------------------------------------------------------
##Oxique
totalOxiqueTable <- data_seq_tax %>% select(OTU_Id,TotalOxique,Division,Zone)
totalOxiqueTable$Zone <- rep("Oxique",each = nrow(totalOxiqueTable))
totalOxiqueTable <- totalOxiqueTable %>% filter(TotalOxique > 0.05*sum(totalOxiqueTable$TotalOxique))
colnames(totalOxiqueTable)[2]  <- "value"
##Anoxique
totalAnoxiqueTable <- data_seq_tax %>% select(OTU_Id,TotalAnoxique,Division,Zone)
totalAnoxiqueTable$Zone <- rep("Anoxique",each = nrow(totalAnoxiqueTable))
totalAnoxiqueTable <- totalAnoxiqueTable %>% filter(TotalAnoxique > 0.05*sum(totalAnoxiqueTable$TotalAnoxique))
colnames(totalAnoxiqueTable)[2]  <- "value"
##Bind
totalZoneTable <- rbind(totalOxiqueTable,totalAnoxiqueTable)
totalZoneTable <- merge(totalZoneTable,taxo, by = "OTU_Id")
##Table
write.table(totalZoneTable, file = "TableOnly/Total_Zone.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(totalZoneOTU, file = "TableOnly/Total_Zone_OTU.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(totalZoneSequence, file = "TableOnly/Total_Zone_Sequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  # Fraction Total -------------------------------------------------------------------
##Petite
totalPetiteTable <- data_seq_tax %>% select(OTU_Id,TotalPetite,Division,Fraction)
totalPetiteTable$Fraction <- rep("Petite",each = nrow(totalPetiteTable))
totalPetiteTable <- totalPetiteTable %>% filter(TotalPetite > 0.05*sum(totalPetiteTable$TotalPetite))
colnames(totalPetiteTable)[2]  <- "value"
##Grande
totalGrandeTable <- data_seq_tax %>% select(OTU_Id,TotalGrande,Division,Fraction)
totalGrandeTable$Fraction <- rep("Grande",each = nrow(totalGrandeTable))
totalGrandeTable <- totalGrandeTable %>% filter(TotalGrande > 0.05*sum(totalGrandeTable$TotalGrande))
colnames(totalGrandeTable)[2]  <- "value"
##Bind
totalFractionTable <- rbind(totalPetiteTable,totalGrandeTable)
totalFractionTable <- merge(totalFractionTable,taxo, by = "OTU_Id")
##Table
write.table(totalFractionTable, file = "TableOnly/Total_Fraction.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(totalFractionOTU, file = "TableOnly/Total_Fraction_OTU.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(totalFractionSequence, file = "TableOnly/Total_Fraction_Sequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  # Date Total -------------------------------------------------------------------
total04Table <- data_seq_tax %>% select(OTU_Id,Total04,Division)
total04Table$Dates <- rep(04,each = nrow(total04Table))
total04Table <- total04Table %>% filter(Total04 > 0.05*sum(total04Table$Total04))
colnames(total04Table)[2]  <- "value"
##
total06Table <- data_seq_tax %>% select(OTU_Id,Total06,Division)
total06Table$Dates <- rep(06,each = nrow(total06Table))
total06Table <- total06Table %>% filter(Total06 > 0.05*sum(total06Table$Total06))
colnames(total06Table)[2]  <- "value"
##
total09Table <- data_seq_tax %>% select(OTU_Id,Total09,Division)
total09Table$Dates <- rep(09,each = nrow(total09Table))
total09Table <- total09Table %>% filter(Total09 > 0.05*sum(total09Table$Total09))
colnames(total09Table)[2]  <- "value"
##
total11Table <- data_seq_tax %>% select(OTU_Id,Total11,Division)
total11Table$Dates <- rep(11,each = nrow(total11Table))
total11Table <- total11Table %>% filter(Total11 > 0.05*sum(total11Table$Total11))
colnames(total11Table)[2]  <- "value"
##
totalDateTable <- rbind(total04Table,total06Table,total09Table,total11Table)
totalDateTable <- merge(totalDateTable,taxo, by = "OTU_Id")
##
write.table(totalDateTable, file = "TableOnly/Total_Dates.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(totalDateOTU, file = "TableOnly/Total_Dates_OTU.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(totalDateSequence, file = "TableOnly/Total_Dates_Sequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

