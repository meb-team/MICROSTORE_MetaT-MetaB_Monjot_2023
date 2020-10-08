#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# XX/XX/XXXX
#
# Script Composition

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

# Set directory and import package -----------------------------------------------------------
setwd("..")

pkg <- c("ggplot2", "readxl","dplyr","tidyr","cowplot","FactoMineR","factoextra","reshape2","varhandle","ggrepel","ggpubr","ggsci","scales","hrbrthemes","GUniFrac","svglite")
lapply(pkg, require, character.only = TRUE)

palette <- c(pal_locuszoom(alpha = 0.8)(7), pal_lancet(alpha = 0.8)(7))
show_col(palette)
palette <- c("#D43F3ACC","#EEA236CC","#AD002ACC","#46B8DACC","#357EBDCC","#9632B8CC","#B8B8B8CC","#00468BCC","#ED0000CC","#42B540CC","#0099B4CC","#925E9FCC")

# Input OTU Table ---------------------------------------------------------
tableV4095 <- read.csv(file = "dataPANAM/PANAM2/V4-result-095/OTU_distribution_tax.txt", sep = "\t")
# Prepare data inf --------------------------------------------------------
infdataini <- read.table(file = "rawdata/data-inf.txt", sep = "\t", header = FALSE)
infdataini <- infdataini[,-1]
infdataini$Rep <- rep("_1", each = 114)
infdataini$variable <- paste(infdataini$V3,infdataini$Rep, sep = "")
infdataini <- infdataini[,-5]
infdataini <- infdataini[,-2]
infdatainif <- separate(infdataini, variable, c("Conditions","Dates","Replicats"), sep = "_")
x <- 0
for (i in infdatainif$Replicats) { 
  x <- x+1
  if (i == "01") infdatainif[x,6] <- "1"
  if (i == "02") infdatainif[x,6] <- "2"
}
infdatainif$OSTA <- rep("OSTA", each = 114)
infdatainif$variable <- paste(infdatainif$V2,infdatainif$OSTA, sep = "")
infdatainif <- infdatainif[,-1]
infdatainif <- infdatainif[,-6]
infdatainif <- separate(infdatainif, variable, c("Cin","variable"), sep = "_")
infdatainif <- infdatainif[,-6]
infdataini<-filter(infdatainif, V5 == "V4")
infdataini[infdataini == "DJOG"] <- "D;J;O;G"
infdataini[infdataini == "DJOP"] <- "D;J;O;P"
infdataini[infdataini == "DJAG"] <- "D;J;A;G"
infdataini[infdataini == "DJAP"] <- "D;J;A;P"
infdataini[infdataini == "DNOG"] <- "D;N;O;G"
infdataini[infdataini == "DNOP"] <- "D;N;O;P"
infdataini[infdataini == "DNAG"] <- "D;N;A;G"
infdataini[infdataini == "DNAP"] <- "D;N;A;P"
infdataini<-as.data.frame(infdataini)
infdataini <- separate(infdataini, Conditions, c("ADN","Jour_Nuit","Oxique_Anoxique","Grande_Petite"),sep = ";")
pattern <- c("variable","V4","V5","ADN","Jour_Nuit","Oxique_Anoxique","Grande_Petite","Dates","Replicats")
samples_df <- infdataini[,pattern]
samples_df$Conditions <- paste(samples_df$ADN,samples_df$Jour_Nuit,samples_df$Oxique_Anoxique,samples_df$Grande_Petite, sep = "")
colnames(samples_df) <- c("sample","Technologie","Regions","ADN","Cycles","Fraction-Oxygène","Fraction-Taille","Dates","Replicats","Conditions")
w <- 0
for (i in samples_df$Cycles) { w <- w +1
if (i == "J") { samples_df[w,"Cycles"] <- "Jour"}
if (i == "N") { samples_df[w,"Cycles"] <- "Nuit"}}
w <- 0
for (i in samples_df$`Fraction-Oxygène`) { w <- w +1
if (i == "O") { samples_df[w,"Fraction-Oxygène"] <- "Oxique"}
if (i == "A") { samples_df[w,"Fraction-Oxygène"] <- "Anoxique"}}
w <- 0
for (i in samples_df$`Fraction-Taille`) { w <- w +1
if (i == "G") { samples_df[w,"Fraction-Taille"] <- "Grande"}
if (i == "P") { samples_df[w,"Fraction-Taille"] <- "Petite"}}
# Remove FR sample (V9 in reality and not V4)
FR <- grep(samples_df[,"sample"], pattern = "FROSTA")
samples_df <- samples_df[-FR,]

# Prepare  Object -------------------------------------------------
  # OTU+SUM ---------------------------------------------------------------------
pattern <- c(grep(pattern = "OSTA", colnames(tableV4095), value = FALSE, fixed = FALSE))
otu_mat <- tableV4095 %>% select(OTU_Id,pattern)
row.names(otu_mat)<-otu_mat$OTU_Id
pattern <- c(grep(pattern = "OSTA", colnames(otu_mat), value = FALSE, fixed = FALSE))
otu_mat[,pattern] <- lapply(otu_mat[,pattern], as.numeric)

# Prep table Corres
tblc <- samples_df[,-2:-7]
tblc$Paste <- paste(tblc$Conditions,tblc$Dates, sep = "_")
tblc <- tblc[,c(-4,-2)]
tblc1 <- tblc %>% filter(tblc$Rep == 1)
tblc2 <- tblc %>% filter(tblc$Rep == 2)
tblcx <- merge(x = tblc1,y = tblc2, by = "Paste")
norep <- c(grep(pattern = "DMOSTA", tblc1$sample),
           grep(pattern = "DVOSTA", tblc1$sample),
           grep(pattern = "DWOSTA", tblc1$sample),
           grep(pattern = "DXOSTA", tblc1$sample),
           grep(pattern = "FOOSTA", tblc1$sample),
           grep(pattern = "FHOSTA", tblc1$sample),
           grep(pattern = "ESOSTA", tblc1$sample))
tblcnorep <- tblc1[norep,]
tblcnorep <- merge(tblcnorep,tblcnorep, by = "Paste")
tblcx <- rbind(tblcx, tblcnorep)
# Pool
w <- 0
data_pool <- as.data.frame(tableV4095[,1])
colnames(data_pool) <- "OTU_Id"
row.names(data_pool) <- data_pool$OTU_Id
for ( h in tblcx[,"sample.x"] ) { w <- w+1
f <- tblcx %>% filter(sample.x == h) %>% select(sample.y)
f <- f$sample.y
print (h)
print (f)
for (r in row.names(otu_mat)) {
  if (otu_mat[r,h] * otu_mat[r,f] == 0) { data_pool[r,h] <- round((otu_mat[r,h]+otu_mat[r,f]))}
  if (otu_mat[r,h] * otu_mat[r,f] != 0) { data_pool[r,h] <- round((otu_mat[r,h]+otu_mat[r,f])/2)}
}}
  # Rarefy ---------------------------------------------------------------------
data_pool <- data_pool %>% select(-OTU_Id)
data_poolx <- t(data_pool)
# Ref : Jun Chen et al. (2012). Associating microbiome composition with environmental covariates using generalized UniFrac distances. 28(16): 2106–2113.
otu_mat_rare <- as.data.frame(t(Rarefy(data_poolx)$otu.tab.rff))
  # PA-AB ------------------------------------------------------------------
#Sequence AB
raw_Sequence_Norm <- otu_mat_rare
raw_Sequence_Norm$OTU_Id <- row.names(raw_Sequence_Norm)
#OTU PA
raw_OTU_Norm <- otu_mat_rare
pattern <- c(grep(pattern = "OSTA", colnames(raw_OTU_Norm), value = FALSE, fixed = FALSE))
raw_OTU_Norm[,pattern][raw_OTU_Norm[,pattern] != 0] <- 1
raw_OTU_Norm$OTU_Id <- row.names(raw_OTU_Norm)
#Stat Rarefy
  #Sequence
avRarefyS <- as.data.frame(colSums(data_pool))
apRarefyS <- as.data.frame(colSums(raw_Sequence_Norm[,pattern]))
statRarefy <- cbind(avRarefyS,apRarefyS)
  #OTU
otu_mat_uniq <- data_pool[,pattern]
otu_mat_uniq[otu_mat_uniq != 0] <- 1
avRarefyO <- as.data.frame(colSums(otu_mat_uniq))
apRarefyO <- as.data.frame(colSums(raw_OTU_Norm[,pattern]))
statRarefy <- cbind(statRarefy,avRarefyO,apRarefyO)
statRarefy["Total",]<-colSums(statRarefy)
colnames(statRarefy) <- c("avRarefy-Sequence","apRarefy-Sequence","avRarefy-OTU","apRarefy-OTU")

OTUav <- otu_mat_uniq
OTUav$OTU_Id <- row.names(OTUav)
pattern <- c(grep(pattern = "OSTA", colnames(OTUav), value = FALSE, fixed = FALSE))
for ( i in row.names(OTUav)) { if (sum(OTUav[i,pattern]) == 0) { OTUav[i,"OTU_Id"] <- "Uniq"}}

OTUap <- raw_OTU_Norm
OTUap$OTU_Id <- row.names(OTUap)
pattern <- c(grep(pattern = "OSTA", colnames(OTUap), value = FALSE, fixed = FALSE))
for ( i in row.names(OTUap)) { if (sum(OTUap[i,pattern]) == 0) { OTUap[i,"OTU_Id"] <- "Uniq"}}
OTUap <- OTUap %>% filter(OTU_Id != "Uniq")

statRarefy["Total","avRarefy-OTU"] <- nrow(OTUav)
statRarefy["Total","apRarefy-OTU"] <- nrow(OTUap)
write.table(statRarefy, file = "Analyse-Composition-Rarefy/TableOnly/StatRarefy_withoutDuplicat.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

#Figure

statRarefyOTU <- statRarefy %>% select("avRarefy-OTU","apRarefy-OTU")
colnames(statRarefyOTU) <- c("NON","OUI")
statRarefyOTU$Echantillons <- row.names(statRarefyOTU)
statRarefymelt <- melt(statRarefyOTU[pattern,])
my_comp <- list(c("OUI","NON"))

svglite("Analyse-Composition-Rarefy/Figure-Sum/Stat-Rarefy-V4_095.svg",width = 3.00,height = 4.00)
R <- ggplot(statRarefymelt, aes(y = value, x = variable)) + geom_boxplot() + geom_point(aes(color = Echantillons), size = 2.5) + 
  stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
  stat_compare_means(method="wilcox.test", label.y = max(statRarefymelt$value)+0.12*max(statRarefymelt$value)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  labs(x="Raréfaction",y="Nombre d'OTUs",x) #+ guides(color = FALSE)
R <- R + theme(legend.position="none")
print(R)
dev.off() 

  # Sum files ---------------------------------------------------------------
pattern <- c(grep(pattern = "OSTA", colnames(raw_Sequence_Norm), value = TRUE, fixed = FALSE)) 
Sum <- as.data.frame(colSums((raw_Sequence_Norm[,pattern])))
colnames(Sum)<- "Séquences"
Sum$OTUs <- colSums((raw_OTU_Norm[,pattern]))
Fusion <- c(grep(pattern = "DMOSTA", row.names(Sum)),
            grep(pattern = "DVOSTA", row.names(Sum)),
            grep(pattern = "DWOSTA", row.names(Sum)),
            grep(pattern = "DXOSTA", row.names(Sum)),
            grep(pattern = "FOOSTA", row.names(Sum)),
            grep(pattern = "FHOSTA", row.names(Sum)),
            grep(pattern = "ESOSTA", row.names(Sum)))
Sum[,"Fusion"] <- "OUI"
Sum[Fusion,"Fusion"] <- "NON"
Sum$Echantillons <- rownames(Sum)
#Plot
#Après rarefaction
# OTU
my_comp <- list(c("OUI","NON"))
svglite("Analyse-Composition-Rarefy/Figure-Sum/Analyse-SumApRare-OTU-V4_095.svg",width = 3.00,height = 4.00)
J <- ggplot(Sum, aes(y = OTUs, x = Fusion)) + geom_boxplot() + geom_point(aes(color = Echantillons), size = 2.5) + 
  stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
  stat_compare_means(method="wilcox.test", label.y = max(Sum$OTUs)+0.12*max(Sum$OTUs)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  labs(x="Regroupement",y="Nombre d'OTUs",x) #+ guides(color = FALSE)
J <- J + theme(legend.position="none")
print(J)
dev.off()  
# Avant rarefaction
  #Séquence
SequenceavRare <- statRarefy 
SequenceavRare$Echantillons <- rownames(SequenceavRare)
SequenceavRare <- SequenceavRare %>% filter(Echantillons != "Total") %>% select(Echantillons,`avRarefy-Sequence`)
Sum <- merge(Sum,SequenceavRare, by = "Echantillons")
my_comp <- list(c("OUI","NON"))
svglite("Analyse-Composition-Rarefy/Figure-Sum/Analyse-SumAvRare-Séquences-V4_095.svg",width = 3.00,height = 4.00)
G <- ggplot(Sum, aes(y = `avRarefy-Sequence`, x = Fusion)) + geom_boxplot() + geom_point(aes(color = Echantillons), size = 2.5) + 
  stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
  stat_compare_means(method="wilcox.test", label.y = max(Sum$`avRarefy-Sequence`)+0.12*max(Sum$`avRarefy-Sequence`)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  labs(x="Regroupement",y="Nombre de séquences",x) #+ guides(color = FALSE)
G <- G + theme(legend.position="none")
print(G)
dev.off()  
#OTU
OTUavRare <- statRarefy 
OTUavRare$Echantillons <- rownames(OTUavRare)
OTUavRare <- OTUavRare %>% filter(Echantillons != "Total") %>% select(Echantillons,`avRarefy-OTU`)
Sum <- merge(Sum,OTUavRare, by = "Echantillons")
my_comp <- list(c("OUI","NON"))
svglite("Analyse-Composition-Rarefy/Figure-Sum/Analyse-SumAvRare-OTU-V4_095.svg",width = 3.00,height = 4.00)
H <- ggplot(Sum, aes(y = `avRarefy-OTU`, x = Fusion)) + geom_boxplot() + geom_point(aes(color = Echantillons), size = 2.5) + 
  stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
  stat_compare_means(method="wilcox.test", label.y = max(Sum$`avRarefy-OTU`)+0.12*max(Sum$`avRarefy-OTU`)) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  labs(x="Regroupement",y="Nombre d'OTUs",x) #+ guides(color = FALSE)
H <- H + theme(legend.position="none")
print(H)
dev.off()  
# Create sort Condition pattern -------------------------------------------
#Cycles
patternCyclesJour <- samples_df %>% filter(Cycles == "Jour") %>% filter(Replicats == 1)
patternCyclesJour <- patternCyclesJour$sample
patternCyclesNuit <- samples_df %>% filter(Cycles == "Nuit") %>% filter(Replicats == 1)
patternCyclesNuit <- patternCyclesNuit$sample
#Fraction oxygène
patternFractionOxique <- samples_df %>% filter(`Fraction-Oxygène` == "Oxique") %>% filter(Replicats == 1)
patternFractionOxique <- patternFractionOxique$sample
patternFractionAnoxique <- samples_df %>% filter(`Fraction-Oxygène` == "Anoxique") %>% filter(Replicats == 1)
patternFractionAnoxique <- patternFractionAnoxique$sample
#Fraction taille
patternFractionPetite <- samples_df %>% filter(`Fraction-Taille` == "Petite") %>% filter(Replicats == 1)
patternFractionPetite <- patternFractionPetite$sample
patternFractionGrande <- samples_df %>% filter(`Fraction-Taille` == "Grande") %>% filter(Replicats == 1)
patternFractionGrande <- patternFractionGrande$sample
#Date
patternDates04 <- samples_df %>% filter(Dates == "04") %>% filter(Replicats == 1)
patternDates04 <- patternDates04$sample
patternDates06 <- samples_df %>% filter(Dates == "06") %>% filter(Replicats == 1)
patternDates06 <- patternDates06$sample
patternDates09 <- samples_df %>% filter(Dates == "09") %>% filter(Replicats == 1)
patternDates09 <- patternDates09$sample
patternDates11 <- samples_df %>% filter(Dates == "11") %>% filter(Replicats == 1)
patternDates11 <- patternDates11$sample



# Create dataframe Condition --------------------------------------------
  # Sequence ----------------------------------------------------------------
    # Condition Cycles ---------------------------------------------------------
# Sequence
#Jour
Jour_Sequence_Norm <- raw_Sequence_Norm %>% select(OTU_Id,patternCyclesJour)
Jour_Sequence_Norm$TotalJour <- rep(0,each = nrow(Jour_Sequence_Norm))
for (i in row.names(Jour_Sequence_Norm)) {Jour_Sequence_Norm[i,"TotalJour"] <- sum(raw_Sequence_Norm[i,patternCyclesJour])}
Jour_Sequence_Norm <- Jour_Sequence_Norm %>% select(-patternCyclesJour)
#Nuit
Nuit_Sequence_Norm <- raw_Sequence_Norm %>% select(OTU_Id,patternCyclesNuit)
Nuit_Sequence_Norm$TotalNuit <- rep(0,each = nrow(Nuit_Sequence_Norm))
for (i in row.names(Nuit_Sequence_Norm)) {Nuit_Sequence_Norm[i,"TotalNuit"] <- sum(raw_Sequence_Norm[i,patternCyclesNuit])}
Nuit_Sequence_Norm <- Nuit_Sequence_Norm %>% select(-patternCyclesNuit)
#Cycles
Cycles_Sequence_Norm <- merge(x = Jour_Sequence_Norm,y = Nuit_Sequence_Norm, by = "OTU_Id")

    # Condition Fraction Oxygène ---------------------------------------------------------
# Sequence
#Oxique
Oxique_Sequence_Norm <- raw_Sequence_Norm %>% select(OTU_Id,patternFractionOxique)
Oxique_Sequence_Norm$TotalOxique <- rep(0,each = nrow(Oxique_Sequence_Norm))
for (i in row.names(Oxique_Sequence_Norm)) {Oxique_Sequence_Norm[i,"TotalOxique"] <- sum(raw_Sequence_Norm[i,patternFractionOxique])}
Oxique_Sequence_Norm <- Oxique_Sequence_Norm %>% select(-patternFractionOxique)
#Anoxique
Anoxique_Sequence_Norm <- raw_Sequence_Norm %>% select(OTU_Id,patternFractionAnoxique)
Anoxique_Sequence_Norm$TotalAnoxique <- rep(0,each = nrow(Anoxique_Sequence_Norm))
for (i in row.names(Anoxique_Sequence_Norm)) {Anoxique_Sequence_Norm[i,"TotalAnoxique"] <- sum(raw_Sequence_Norm[i,patternFractionAnoxique])}
Anoxique_Sequence_Norm <- Anoxique_Sequence_Norm %>% select(-patternFractionAnoxique)
#FractionO
FractionO_Sequence_Norm <- merge(x = Oxique_Sequence_Norm,y = Anoxique_Sequence_Norm, by = "OTU_Id")


    # Condition Fraction Taille ---------------------------------------------------------
# Sequence
#Petite
Petite_Sequence_Norm <- raw_Sequence_Norm %>% select(OTU_Id,patternFractionPetite)
Petite_Sequence_Norm$TotalPetite <- rep(0,each = nrow(Petite_Sequence_Norm))
for (i in row.names(Petite_Sequence_Norm)) {Petite_Sequence_Norm[i,"TotalPetite"] <- sum(raw_Sequence_Norm[i,patternFractionPetite])}
Petite_Sequence_Norm <- Petite_Sequence_Norm %>% select(-patternFractionPetite)
#Grande
Grande_Sequence_Norm <- raw_Sequence_Norm %>% select(OTU_Id,patternFractionGrande)
Grande_Sequence_Norm$TotalGrande <- rep(0,each = nrow(Grande_Sequence_Norm))
for (i in row.names(Grande_Sequence_Norm)) {Grande_Sequence_Norm[i,"TotalGrande"] <- sum(raw_Sequence_Norm[i,patternFractionGrande])}
Grande_Sequence_Norm <- Grande_Sequence_Norm %>% select(-patternFractionGrande)
#FractionT
FractionT_Sequence_Norm <- merge(x = Petite_Sequence_Norm,y = Grande_Sequence_Norm, by = "OTU_Id")



    # Condition Dates ---------------------------------------------------------
# Sequence
#04
`04_Sequence_Norm` <- raw_Sequence_Norm %>% select(OTU_Id,patternDates04)
`04_Sequence_Norm`$Total04 <- rep(0,each = nrow(`04_Sequence_Norm`))
for (i in row.names(`04_Sequence_Norm`)) {`04_Sequence_Norm`[i,"Total04"] <- sum(raw_Sequence_Norm[i,patternDates04])}
`04_Sequence_Norm` <- `04_Sequence_Norm` %>% select(-patternDates04)
#06
`06_Sequence_Norm` <- raw_Sequence_Norm %>% select(OTU_Id,patternDates06)
`06_Sequence_Norm`$Total06 <- rep(0,each = nrow(`06_Sequence_Norm`))
for (i in row.names(`06_Sequence_Norm`)) {`06_Sequence_Norm`[i,"Total06"] <- sum(raw_Sequence_Norm[i,patternDates06])}
`06_Sequence_Norm` <- `06_Sequence_Norm` %>% select(-patternDates06)
#09
`09_Sequence_Norm` <- raw_Sequence_Norm %>% select(OTU_Id,patternDates09)
`09_Sequence_Norm`$Total09 <- rep(0,each = nrow(`09_Sequence_Norm`))
for (i in row.names(`09_Sequence_Norm`)) {`09_Sequence_Norm`[i,"Total09"] <- sum(raw_Sequence_Norm[i,patternDates09])}
`09_Sequence_Norm` <- `09_Sequence_Norm` %>% select(-patternDates09)
#11
`11_Sequence_Norm` <- raw_Sequence_Norm %>% select(OTU_Id,patternDates11)
`11_Sequence_Norm`$Total11 <- rep(0,each = nrow(`11_Sequence_Norm`))
for (i in row.names(`11_Sequence_Norm`)) {`11_Sequence_Norm`[i,"Total11"] <- sum(raw_Sequence_Norm[i,patternDates11])}
`11_Sequence_Norm` <- `11_Sequence_Norm` %>% select(-patternDates11)
#FractionT
Dates_Sequence_Norm <- merge(x = `04_Sequence_Norm`,y = `06_Sequence_Norm`, by = "OTU_Id")
Dates_Sequence_Norm <- merge(x = `Dates_Sequence_Norm`,y = `09_Sequence_Norm`, by = "OTU_Id")
Dates_Sequence_Norm <- merge(x = `Dates_Sequence_Norm`,y = `11_Sequence_Norm`, by = "OTU_Id")
  # OTU ----------------------------------------------------------------
    # Condition Cycles ---------------------------------------------------------
# OTU
#Jour
Jour_OTU_Norm <- raw_OTU_Norm %>% select(OTU_Id,patternCyclesJour)
Jour_OTU_Norm$TotalJour <- rep(0,each = nrow(Jour_OTU_Norm))
for (i in row.names(Jour_OTU_Norm)) {if (sum(raw_OTU_Norm[i,patternCyclesJour]) > 0) {Jour_OTU_Norm[i,"TotalJour"] <- 1}}
Jour_OTU_Norm <- Jour_OTU_Norm %>% select(-patternCyclesJour)
#Nuit
Nuit_OTU_Norm <- raw_OTU_Norm %>% select(OTU_Id,patternCyclesNuit)
Nuit_OTU_Norm$TotalNuit <- rep(0,each = nrow(Nuit_OTU_Norm))
for (i in row.names(Nuit_OTU_Norm)) {if (sum(raw_OTU_Norm[i,patternCyclesNuit]) > 0) {Nuit_OTU_Norm[i,"TotalNuit"] <- 1}}
Nuit_OTU_Norm <- Nuit_OTU_Norm %>% select(-patternCyclesNuit)
#Cycles
Cycles_OTU_Norm <- merge(x = Jour_OTU_Norm,y = Nuit_OTU_Norm, by = "OTU_Id")

    # Condition Fraction Oxygène ---------------------------------------------------------
# OTU
#Oxique
Oxique_OTU_Norm <- raw_OTU_Norm %>% select(OTU_Id,patternFractionOxique)
Oxique_OTU_Norm$TotalOxique <- rep(0,each = nrow(Oxique_OTU_Norm))
for (i in row.names(Oxique_OTU_Norm)) {if (sum(raw_OTU_Norm[i,patternFractionOxique]) > 0) {Oxique_OTU_Norm[i,"TotalOxique"] <- 1}}
Oxique_OTU_Norm <- Oxique_OTU_Norm %>% select(-patternFractionOxique)
#Anoxique
Anoxique_OTU_Norm <- raw_OTU_Norm %>% select(OTU_Id,patternFractionAnoxique)
Anoxique_OTU_Norm$TotalAnoxique <- rep(0,each = nrow(Anoxique_OTU_Norm))
for (i in row.names(Anoxique_OTU_Norm)) {if (sum(raw_OTU_Norm[i,patternFractionAnoxique]) > 0) {Anoxique_OTU_Norm[i,"TotalAnoxique"] <- 1}}
Anoxique_OTU_Norm <- Anoxique_OTU_Norm %>% select(-patternFractionAnoxique)
#FractionO
FractionO_OTU_Norm <- merge(x = Oxique_OTU_Norm,y = Anoxique_OTU_Norm, by = "OTU_Id")


    # Condition Fraction Taille ---------------------------------------------------------
# OTU
#Petite
Petite_OTU_Norm <- raw_OTU_Norm %>% select(OTU_Id,patternFractionPetite)
Petite_OTU_Norm$TotalPetite <- rep(0,each = nrow(Petite_OTU_Norm))
for (i in row.names(Petite_OTU_Norm)) {if (sum(raw_OTU_Norm[i,patternFractionPetite]) > 0) {Petite_OTU_Norm[i,"TotalPetite"] <- 1}}
Petite_OTU_Norm <- Petite_OTU_Norm %>% select(-patternFractionPetite)
#Grande
Grande_OTU_Norm <- raw_OTU_Norm %>% select(OTU_Id,patternFractionGrande)
Grande_OTU_Norm$TotalGrande <- rep(0,each = nrow(Grande_OTU_Norm))
for (i in row.names(Grande_OTU_Norm)) {if (sum(raw_OTU_Norm[i,patternFractionGrande]) > 0) {Grande_OTU_Norm[i,"TotalGrande"] <- 1}}
Grande_OTU_Norm <- Grande_OTU_Norm %>% select(-patternFractionGrande)
#FractionT
FractionT_OTU_Norm <- merge(x = Petite_OTU_Norm,y = Grande_OTU_Norm, by = "OTU_Id")



    # Condition Dates ---------------------------------------------------------
# OTU
#04
`04_OTU_Norm` <- raw_OTU_Norm %>% select(OTU_Id,patternDates04)
`04_OTU_Norm`$Total04 <- rep(0,each = nrow(`04_OTU_Norm`))
for (i in row.names(`04_OTU_Norm`)) {if (sum(raw_OTU_Norm[i,patternDates04]) > 0) {`04_OTU_Norm`[i,"Total04"] <- 1}}
`04_OTU_Norm` <- `04_OTU_Norm` %>% select(-patternDates04)
#06
`06_OTU_Norm` <- raw_OTU_Norm %>% select(OTU_Id,patternDates06)
`06_OTU_Norm`$Total06 <- rep(0,each = nrow(`06_OTU_Norm`))
for (i in row.names(`06_OTU_Norm`)) {if (sum(raw_OTU_Norm[i,patternDates06]) > 0) {`06_OTU_Norm`[i,"Total06"] <- 1}}
`06_OTU_Norm` <- `06_OTU_Norm` %>% select(-patternDates06)
#09
`09_OTU_Norm` <- raw_OTU_Norm %>% select(OTU_Id,patternDates09)
`09_OTU_Norm`$Total09 <- rep(0,each = nrow(`09_OTU_Norm`))
for (i in row.names(`09_OTU_Norm`)) {if (sum(raw_OTU_Norm[i,patternDates09]) > 0) {`09_OTU_Norm`[i,"Total09"] <- 1}}
`09_OTU_Norm` <- `09_OTU_Norm` %>% select(-patternDates09)
#11
`11_OTU_Norm` <- raw_OTU_Norm %>% select(OTU_Id,patternDates11)
`11_OTU_Norm`$Total11 <- rep(0,each = nrow(`11_OTU_Norm`))
for (i in row.names(`11_OTU_Norm`)) {if (sum(raw_OTU_Norm[i,patternDates11]) > 0) {`11_OTU_Norm`[i,"Total11"] <- 1}}
`11_OTU_Norm` <- `11_OTU_Norm` %>% select(-patternDates11)
#FractionT
Dates_OTU_Norm <- merge(x = `04_OTU_Norm`,y = `06_OTU_Norm`, by = "OTU_Id")
Dates_OTU_Norm <- merge(x = `Dates_OTU_Norm`,y = `09_OTU_Norm`, by = "OTU_Id")
Dates_OTU_Norm <- merge(x = `Dates_OTU_Norm`,y = `11_OTU_Norm`, by = "OTU_Id")




# AFC Plot General ----------------------------------------------------------------
  # Sequence ----------------------------------------------------------------
dt <- as.data.frame(raw_Sequence_Norm)
dt <- dt %>% select(-OTU_Id)
#CA
res.ca <- CA(dt, graph = FALSE,ncp = 2 )
fviz_ca_row(res.ca, repel = FALSE, label = "none")
p <- get_ca_row(res.ca)
coord <- p$coord
coord <- as.data.frame(coord)
coord$OTU_Id <- row.names(coord)
dataSequence <- merge(x = coord, y = Cycles_Sequence_Norm, by = "OTU_Id")
dataSequence$Jour <- rep(0,each = nrow(dataSequence))
dataSequence$Nuit <- rep(0,each = nrow(dataSequence))
for (i in rownames(dataSequence)) { if ( dataSequence[i,"TotalJour"] != 0) { dataSequence[i,"Jour"] <- 1}}
for (i in rownames(dataSequence)) { if ( dataSequence[i,"TotalNuit"] != 0) { dataSequence[i,"Nuit"] <- 2}}
dataSequence$Cycles <- rep(0,each = nrow(dataSequence))
for (i in rownames(dataSequence)) { dataSequence[i,"Cycles"] <- dataSequence[i,"Jour"] + dataSequence[i,"Nuit"] 
if ( dataSequence[i,"Cycles"] == 1) { dataSequence[i,"Cycles"] <- "Jour"}
if ( dataSequence[i,"Cycles"] == 2) { dataSequence[i,"Cycles"] <- "Nuit"}
if ( dataSequence[i,"Cycles"] == 3) { dataSequence[i,"Cycles"] <- "Communs"}}
#FractionO
dataSequence <- merge(x = dataSequence, y = FractionO_Sequence_Norm, by = "OTU_Id")
dataSequence$Oxique <- rep(0,each = nrow(dataSequence))
dataSequence$Anoxique <- rep(0,each = nrow(dataSequence))
for (i in rownames(dataSequence)) { if ( dataSequence[i,"TotalOxique"] != 0) { dataSequence[i,"Oxique"] <- 1}}
for (i in rownames(dataSequence)) { if ( dataSequence[i,"TotalAnoxique"] != 0) { dataSequence[i,"Anoxique"] <- 2}}
dataSequence$FractionO <- rep(0,each = nrow(dataSequence))
for (i in rownames(dataSequence)) { dataSequence[i,"FractionO"] <- dataSequence[i,"Oxique"] + dataSequence[i,"Anoxique"] 
if ( dataSequence[i,"FractionO"] == 1) { dataSequence[i,"FractionO"] <- "Oxique"}
if ( dataSequence[i,"FractionO"] == 2) { dataSequence[i,"FractionO"] <- "Anoxique"}
if ( dataSequence[i,"FractionO"] == 3) { dataSequence[i,"FractionO"] <- "Communs"}}
#FractionT
dataSequence <- merge(x = dataSequence, y = FractionT_Sequence_Norm, by = "OTU_Id")
dataSequence$Petite <- rep(0,each = nrow(dataSequence))
dataSequence$Grande <- rep(0,each = nrow(dataSequence))
for (i in rownames(dataSequence)) { if ( dataSequence[i,"TotalPetite"] != 0) { dataSequence[i,"Petite"] <- 1}}
for (i in rownames(dataSequence)) { if ( dataSequence[i,"TotalGrande"] != 0) { dataSequence[i,"Grande"] <- 2}}
dataSequence$FractionT <- rep(0,each = nrow(dataSequence))
for (i in rownames(dataSequence)) { dataSequence[i,"FractionT"] <- dataSequence[i,"Petite"] + dataSequence[i,"Grande"] 
if ( dataSequence[i,"FractionT"] == 1) { dataSequence[i,"FractionT"] <- "Petite"}
if ( dataSequence[i,"FractionT"] == 2) { dataSequence[i,"FractionT"] <- "Grande"}
if ( dataSequence[i,"FractionT"] == 3) { dataSequence[i,"FractionT"] <- "Communs"}}
#Dates
dataSequence <- merge(x = dataSequence, y = Dates_Sequence_Norm, by = "OTU_Id")
dataSequence04 <- dataSequence %>% filter(Total04 != 0) %>% select(-Total06,-Total09,-Total11)
dataSequence06 <- dataSequence %>% filter(Total06 != 0) %>% select(-Total04,-Total09,-Total11)
dataSequence09 <- dataSequence %>% filter(Total09 != 0) %>% select(-Total04,-Total06,-Total11)
dataSequence11 <- dataSequence %>% filter(Total11 != 0) %>% select(-Total04,-Total06,-Total09)

Xseq <- fviz_screeplot (res.ca, addlabels = TRUE, ylim = c(0, 50))
Xseq <- Xseq$data
Dim1Seq <- paste("Dim 1 [",round(Xseq %>% filter(dim == 1) %>% select(eig),1),"%]", sep = "")
Dim2Seq <- paste("Dim 2 [",round(Xseq %>% filter(dim == 2) %>% select(eig),1),"%]", sep = "")

  # OTU ----------------------------------------------------------------
dt <- as.data.frame(raw_OTU_Norm)
dt <- dt %>% select(-OTU_Id)
#CA
res.ca <- CA(dt, graph = FALSE,ncp = 2 )
fviz_ca_row(res.ca, repel = FALSE, label = "none")
p <- get_ca_row(res.ca)
coord <- p$coord
coord <- as.data.frame(coord)
coord$OTU_Id <- row.names(coord)
#Cycles
dataOTU <- merge(x = coord, y = Cycles_OTU_Norm, by = "OTU_Id")
dataOTU$Jour <- rep(0,each = nrow(dataOTU))
dataOTU$Nuit <- rep(0,each = nrow(dataOTU))
for (i in rownames(dataOTU)) { if ( dataOTU[i,"TotalJour"] != 0) { dataOTU[i,"Jour"] <- 1}}
for (i in rownames(dataOTU)) { if ( dataOTU[i,"TotalNuit"] != 0) { dataOTU[i,"Nuit"] <- 2}}
dataOTU$Cycles <- rep(0,each = nrow(dataOTU))
for (i in rownames(dataOTU)) { dataOTU[i,"Cycles"] <- dataOTU[i,"Jour"] + dataOTU[i,"Nuit"] 
if ( dataOTU[i,"Cycles"] == 1) { dataOTU[i,"Cycles"] <- "Jour"}
if ( dataOTU[i,"Cycles"] == 2) { dataOTU[i,"Cycles"] <- "Nuit"}
if ( dataOTU[i,"Cycles"] == 3) { dataOTU[i,"Cycles"] <- "Communs"}}
#FractionO
dataOTU <- merge(x = dataOTU, y = FractionO_OTU_Norm, by = "OTU_Id")
dataOTU$Oxique <- rep(0,each = nrow(dataOTU))
dataOTU$Anoxique <- rep(0,each = nrow(dataOTU))
for (i in rownames(dataOTU)) { if ( dataOTU[i,"TotalOxique"] != 0) { dataOTU[i,"Oxique"] <- 1}}
for (i in rownames(dataOTU)) { if ( dataOTU[i,"TotalAnoxique"] != 0) { dataOTU[i,"Anoxique"] <- 2}}
dataOTU$FractionO <- rep(0,each = nrow(dataOTU))
for (i in rownames(dataOTU)) { dataOTU[i,"FractionO"] <- dataOTU[i,"Oxique"] + dataOTU[i,"Anoxique"] 
if ( dataOTU[i,"FractionO"] == 1) { dataOTU[i,"FractionO"] <- "Oxique"}
if ( dataOTU[i,"FractionO"] == 2) { dataOTU[i,"FractionO"] <- "Anoxique"}
if ( dataOTU[i,"FractionO"] == 3) { dataOTU[i,"FractionO"] <- "Communs"}}
#FractionT
dataOTU <- merge(x = dataOTU, y = FractionT_OTU_Norm, by = "OTU_Id")
dataOTU$Petite <- rep(0,each = nrow(dataOTU))
dataOTU$Grande <- rep(0,each = nrow(dataOTU))
for (i in rownames(dataOTU)) { if ( dataOTU[i,"TotalPetite"] != 0) { dataOTU[i,"Petite"] <- 1}}
for (i in rownames(dataOTU)) { if ( dataOTU[i,"TotalGrande"] != 0) { dataOTU[i,"Grande"] <- 2}}
dataOTU$FractionT <- rep(0,each = nrow(dataOTU))
for (i in rownames(dataOTU)) { dataOTU[i,"FractionT"] <- dataOTU[i,"Petite"] + dataOTU[i,"Grande"] 
if ( dataOTU[i,"FractionT"] == 1) { dataOTU[i,"FractionT"] <- "Petite"}
if ( dataOTU[i,"FractionT"] == 2) { dataOTU[i,"FractionT"] <- "Grande"}
if ( dataOTU[i,"FractionT"] == 3) { dataOTU[i,"FractionT"] <- "Communs"}}
#Dates
dataOTU <- merge(x = dataOTU, y = Dates_OTU_Norm, by = "OTU_Id")
dataOTU04 <- dataOTU %>% filter(Total04 != 0) %>% select(-Total06,-Total09,-Total11)
dataOTU06 <- dataOTU %>% filter(Total06 != 0) %>% select(-Total04,-Total09,-Total11)
dataOTU09 <- dataOTU %>% filter(Total09 != 0) %>% select(-Total04,-Total06,-Total11)
dataOTU11 <- dataOTU %>% filter(Total11 != 0) %>% select(-Total04,-Total06,-Total09)

XOTU <- fviz_screeplot (res.ca, addlabels = TRUE, ylim = c(0, 50))
XOTU <- XOTU$data
Dim1OTU <- paste("Dim 1 [",round(XOTU %>% filter(dim == 1) %>% select(eig),1),"%]", sep = "")
Dim2OTU <- paste("Dim 2 [",round(XOTU %>% filter(dim == 2) %>% select(eig),1),"%]", sep = "")


# AFC Taxonomy ------------------------------------------------------------
  # AFC Sequence ----------------------------------------------------------------
dataSequenceTaxo <- dataSequence
taxo <- tableV4095 %>% select(OTU_Id,NN.taxonomy)
dataSequenceTaxo <- merge(x = dataSequenceTaxo, y =  taxo, by = "OTU_Id")
dataSequenceTaxo <- separate(dataSequenceTaxo, NN.taxonomy, c('Domain', 'Supergroup','Division'), sep = ";")
Parasitemin <- c("Apicomplexa","Perkinsea","Chytridiomycota","Cryptomycota","Ichthyophonida","Cercozoa","Oomycetes","Pirsonia","Labyrinthulida")
dataSequenceTaxo$Parasitism <- rep("Autres",each = nrow(dataSequenceTaxo))
w <- 0
for (j in dataSequenceTaxo$Division) { w <- w+1
for (i in Parasitemin) { if (i == j) { dataSequenceTaxo[w,"Parasitism"] <- "Potentiels parasites" }}}
#dataSequenceTaxo$Cyclesx <- paste(dataSequenceTaxo$Cycles,dataSequenceTaxo$Parasitism, sep = " - ")
dataSequenceTaxo <- dataSequenceTaxo %>% filter(Parasitism == "Potentiels parasites")
    # Figure 3 dates ----------------------------------------------------------
#Figure Cycles
a <- ggplot(dataSequenceTaxo, aes(y = `Dim 2`, x = `Dim 1`, color = Cycles)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = Cycles, alpha = Cycles),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1Seq,y=Dim2Seq,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(a)

#Figure FractionO
b <- ggplot(dataSequenceTaxo, aes(y = `Dim 2`, x = `Dim 1`, color = FractionO)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = FractionO, alpha = FractionO),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1Seq,y=Dim2Seq,color = "Zones", alpha = "Zones", linetype = "Zones")
print(b)
#Figure FractionT
c <- ggplot(dataSequenceTaxo, aes(y = `Dim 2`, x = `Dim 1`, color = FractionT)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = FractionT, alpha = FractionT),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1Seq,y=Dim2Seq,color = "Fractions", alpha = "Fractions", linetype = "Fractions")
print(c)
#Coplot
svglite("Analyse-Composition-Rarefy/AFC-Parasite/Analyse-AFC-Sequence-Parasite-V4_095.svg",width = 10.00,height = 6.00)
Sum_plot <- plot_grid(a, c, b, ncol = 2, nrow = 2, rel_widths = c(3,3), rel_heights = c(3,3))
print(Sum_plot)
dev.off()  
    # Toute Condition Avril ---------------------------------------------------

taxo <- tableV4095 %>% select(OTU_Id,NN.taxonomy)
dataSequenceTaxo04 <- merge(x = dataSequence04, y =  taxo, by = "OTU_Id")
dataSequenceTaxo04 <- separate(dataSequenceTaxo04, NN.taxonomy, c('Domain', 'Supergroup','Division'), sep = ";")
Parasitemin <- c("Apicomplexa","Perkinsea","Chytridiomycota","Cryptomycota","Ichthyophonida","Cercozoa","Oomycetes","Pirsonia","Labyrinthulida")
dataSequenceTaxo04$Parasitism <- rep("Autres",each = nrow(dataSequenceTaxo04))
w <- 0
for (j in dataSequenceTaxo04$Division) { w <- w+1
for (i in Parasitemin) { if (i == j) { dataSequenceTaxo04[w,"Parasitism"] <- "Potentiels parasites" }}}

dataSequenceTaxo04 <-  dataSequenceTaxo04 %>% filter(Parasitism == "Potentiels parasites")



aw <- ggplot(dataSequenceTaxo04, aes(y = `Dim 2`, x = `Dim 1`, color = Cycles)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = Cycles, alpha = Cycles),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1Seq,y=Dim2Seq,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(aw)

#Figure FractionO
bw <- ggplot(dataSequenceTaxo04, aes(y = `Dim 2`, x = `Dim 1`, color = FractionO)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = FractionO, alpha = FractionO),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1Seq,y=Dim2Seq,color = "Zones", alpha = "Zones", linetype = "Zones")
print(bw)

#Figure FractionT
cw <- ggplot(dataSequenceTaxo04, aes(y = `Dim 2`, x = `Dim 1`, color = FractionT)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = FractionT, alpha = FractionT),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1Seq,y=Dim2Seq,color = "Fraction", alpha = "Fraction", linetype = "Fraction")
print(cw)

svglite("Analyse-Composition-Rarefy/AFC-Parasite/Analyse-AFC-Sequence-Parasite-Avril-V4_095.svg",width = 10.00,height = 6.00)
b_plot <- plot_grid(aw, cw, bw, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(b_plot)
dev.off()  
    # Toute Condition Juin ---------------------------------------------------
taxo <- tableV4095 %>% select(OTU_Id,NN.taxonomy)
dataSequenceTaxo06 <- merge(x = dataSequence06, y =  taxo, by = "OTU_Id")
dataSequenceTaxo06 <- separate(dataSequenceTaxo06, NN.taxonomy, c('Domain', 'Supergroup','Division'), sep = ";")
Parasitemin <- c("Apicomplexa","Perkinsea","Chytridiomycota","Cryptomycota","Ichthyophonida","Cercozoa","Oomycetes","Pirsonia","Labyrinthulida")
dataSequenceTaxo06$Parasitism <- rep("Autres",each = nrow(dataSequenceTaxo06))
w <- 0
for (j in dataSequenceTaxo06$Division) { w <- w+1
for (i in Parasitemin) { if (i == j) { dataSequenceTaxo06[w,"Parasitism"] <- "Potentiels parasites" }}}
dataSequenceTaxo06 <-  dataSequenceTaxo06 %>% filter(Parasitism == "Potentiels parasites")


aw <- ggplot(dataSequenceTaxo06, aes(y = `Dim 2`, x = `Dim 1`, color = Cycles)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = Cycles, alpha = Cycles),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1Seq,y=Dim2Seq,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(aw)

#Figure FractionO
bw <- ggplot(dataSequenceTaxo06, aes(y = `Dim 2`, x = `Dim 1`, color = FractionO)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = FractionO, alpha = FractionO),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1Seq,y=Dim2Seq,color = "Zones", alpha = "Zones", linetype = "Zones")
print(bw)

#Figure FractionT
cw <- ggplot(dataSequenceTaxo06, aes(y = `Dim 2`, x = `Dim 1`, color = FractionT)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = FractionT, alpha = FractionT),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1Seq,y=Dim2Seq,color = "Fraction", alpha = "Fraction", linetype = "Fraction")
print(cw)

#Potentiel parasite


svglite("Analyse-Composition-Rarefy/AFC-Parasite/Analyse-AFC-Sequence-Parasite-Juin-V4_095.svg",width = 10.00,height = 6.00)
b_plot <- plot_grid(aw, cw, bw, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(b_plot)
dev.off()  
    # Toute Condition Septembre ---------------------------------------------------
taxo <- tableV4095 %>% select(OTU_Id,NN.taxonomy)
dataSequenceTaxo09 <- merge(x = dataSequence09, y =  taxo, by = "OTU_Id")
dataSequenceTaxo09 <- separate(dataSequenceTaxo09, NN.taxonomy, c('Domain', 'Supergroup','Division'), sep = ";")
Parasitemin <- c("Apicomplexa","Perkinsea","Chytridiomycota","Cryptomycota","Ichthyophonida","Cercozoa","Oomycetes","Pirsonia","Labyrinthulida")
dataSequenceTaxo09$Parasitism <- rep("Autres",each = nrow(dataSequenceTaxo09))
w <- 0
for (j in dataSequenceTaxo09$Division) { w <- w+1
for (i in Parasitemin) { if (i == j) { dataSequenceTaxo09[w,"Parasitism"] <- "Potentiels parasites" }}}
dataSequenceTaxo09 <-  dataSequenceTaxo09 %>% filter(Parasitism == "Potentiels parasites")



aw <- ggplot(dataSequenceTaxo09, aes(y = `Dim 2`, x = `Dim 1`, color = Cycles)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = Cycles, alpha = Cycles),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1Seq,y=Dim2Seq,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(aw)

#Figure FractionO
bw <- ggplot(dataSequenceTaxo09, aes(y = `Dim 2`, x = `Dim 1`, color = FractionO)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = FractionO, alpha = FractionO),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1Seq,y=Dim2Seq,color = "Zones", alpha = "Zones", linetype = "Zones")
print(bw)

#Figure FractionT
cw <- ggplot(dataSequenceTaxo09, aes(y = `Dim 2`, x = `Dim 1`, color = FractionT)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = FractionT, alpha = FractionT),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1Seq,y=Dim2Seq,color = "Fraction", alpha = "Fraction", linetype = "Fraction")
print(cw)

#Potentiel parasite

svglite("Analyse-Composition-Rarefy/AFC-Parasite/Analyse-AFC-Sequence-Parasite-Septembre-V4_095.svg",width = 10.00,height = 6.00)
b_plot <- plot_grid(aw, cw, bw, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(b_plot)
dev.off()  
    # Toute Condition Novembre ---------------------------------------------------
taxo <- tableV4095 %>% select(OTU_Id,NN.taxonomy)
dataSequenceTaxo11 <- merge(x = dataSequence11, y =  taxo, by = "OTU_Id")
dataSequenceTaxo11 <- separate(dataSequenceTaxo11, NN.taxonomy, c('Domain', 'Supergroup','Division'), sep = ";")
Parasitemin <- c("Apicomplexa","Perkinsea","Chytridiomycota","Cryptomycota","Ichthyophonida","Cercozoa","Oomycetes","Pirsonia","Labyrinthulida")
dataSequenceTaxo11$Parasitism <- rep("Autres",each = nrow(dataSequenceTaxo11))
w <- 0
for (j in dataSequenceTaxo11$Division) { w <- w+1
for (i in Parasitemin) { if (i == j) { dataSequenceTaxo11[w,"Parasitism"] <- "Potentiels parasites" }}}
dataSequenceTaxo11 <-  dataSequenceTaxo11 %>% filter(Parasitism == "Potentiels parasites")

aw <- ggplot(dataSequenceTaxo11, aes(y = `Dim 2`, x = `Dim 1`, color = Cycles)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = Cycles, alpha = Cycles),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1Seq,y=Dim2Seq,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(aw)

#Figure FractionO
bw <- ggplot(dataSequenceTaxo11, aes(y = `Dim 2`, x = `Dim 1`, color = FractionO)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = FractionO, alpha = FractionO),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1Seq,y=Dim2Seq,color = "Zones", alpha = "Zones", linetype = "Zones")
print(bw)

#Figure FractionT
cw <- ggplot(dataSequenceTaxo11, aes(y = `Dim 2`, x = `Dim 1`, color = FractionT)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = FractionT, alpha = FractionT),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1Seq,y=Dim2Seq,color = "Fraction", alpha = "Fraction", linetype = "Fraction")
print(cw)

#Potentiel parasite

svglite("Analyse-Composition-Rarefy/AFC-Parasite/Analyse-AFC-Sequence-Parasite-Novembre-V4_095.svg",width = 10.00,height = 6.00)
b_plot <- plot_grid(aw, cw, bw, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(b_plot)
dev.off()  

    # Dates ---------------------------------------------------
dataSequenceTaxoDate <- dataSequenceTaxo
for (i in row.names(dataSequenceTaxoDate)) {dataSequenceTaxoDate[i,"Domain"] <- sum(dataSequenceTaxoDate[i,"Total04"], dataSequenceTaxoDate[i,"Total06"], dataSequenceTaxoDate[i,"Total09"], dataSequenceTaxoDate[i,"Total11"])}
for (i in row.names(dataSequenceTaxoDate)) { if (dataSequenceTaxoDate[i,"Total04"] == 1) { dataSequenceTaxoDate[i,"Total04"] <- "Avril"}}
for (i in row.names(dataSequenceTaxoDate)) { if (dataSequenceTaxoDate[i,"Total06"] == 1) { dataSequenceTaxoDate[i,"Total06"] <- "Juin"}}
for (i in row.names(dataSequenceTaxoDate)) { if (dataSequenceTaxoDate[i,"Total09"] == 1) { dataSequenceTaxoDate[i,"Total09"] <- "Septembre"}}
for (i in row.names(dataSequenceTaxoDate)) { if (dataSequenceTaxoDate[i,"Total11"] == 1) { dataSequenceTaxoDate[i,"Total11"] <- "Novembre"}}
for (i in row.names(dataSequenceTaxoDate)) { if (dataSequenceTaxoDate[i,"Domain"] > 1) { dataSequenceTaxoDate[i,"Domain"] <- "Communs"}}
for (i in row.names(dataSequenceTaxoDate)) { if (dataSequenceTaxoDate[i,"Domain"] == 1) { dataSequenceTaxoDate[i,"Domain"] <- c(dataSequenceTaxoDate[i,"Total04"],dataSequenceTaxoDate[i,"Total06"],dataSequenceTaxoDate[i,"Total09"],dataSequenceTaxoDate[i,"Total11"])[c(dataSequenceTaxoDate[i,"Total04"],dataSequenceTaxoDate[i,"Total06"],dataSequenceTaxoDate[i,"Total09"],dataSequenceTaxoDate[i,"Total11"])!=0]}}
#Figure Cycles
as <- ggplot(dataSequenceTaxoDate, aes(y = `Dim 2`, x = `Dim 1`, color = Domain)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = Domain, alpha = Domain),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2,2,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1,0.1,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF","lightgreen","black")) +
  labs(x=Dim1Seq,y=Dim2Seq,color = "Dates", alpha = "Dates", linetype = "Dates")
print(as)


  # AFC OTU ---------------------------------------------------------------------
dataOTUTaxo <- dataOTU
taxo <- tableV4095 %>% select(OTU_Id,NN.taxonomy)
dataOTUTaxo <- merge(x = dataOTUTaxo, y =  taxo, by = "OTU_Id")
dataOTUTaxo <- separate(dataOTUTaxo, NN.taxonomy, c('Domain', 'Supergroup','Division'), sep = ";")
Parasitemin <- c("Apicomplexa","Perkinsea","Chytridiomycota","Cryptomycota","Ichthyophonida","Cercozoa","Oomycetes","Pirsonia","Labyrinthulida")
dataOTUTaxo$Parasitism <- rep("Autres",each = nrow(dataOTUTaxo))
w <- 0
for (j in dataOTUTaxo$Division) { w <- w+1
for (i in Parasitemin) { if (i == j) { dataOTUTaxo[w,"Parasitism"] <- "Potentiels parasites" }}}
dataOTUTaxo <- dataOTUTaxo %>% filter(Parasitism == "Potentiels parasites")

    # Figure 3 dates ----------------------------------------------------------
#Figure Cycles
as <- ggplot(dataOTUTaxo, aes(y = `Dim 2`, x = `Dim 1`, color = Cycles)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = Cycles, alpha = Cycles),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1OTU,y=Dim2OTU,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(as)

#Figure FractionO
bs <- ggplot(dataOTUTaxo, aes(y = `Dim 2`, x = `Dim 1`, color = FractionO)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = FractionO, alpha = FractionO),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1OTU,y=Dim2OTU,color = "Zones", alpha = "Zones", linetype = "Zones")
print(bs)
#Figure FractionT
cs <- ggplot(dataOTUTaxo, aes(y = `Dim 2`, x = `Dim 1`, color = FractionT)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = FractionT, alpha = FractionT),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1OTU,y=Dim2OTU,color = "Fractions", alpha = "Fractions", linetype = "Fractions")
print(cs)
#Coplot
svglite("Analyse-Composition-Rarefy/AFC-Parasite/Analyse-AFC-OTU-Parasite-V4_095.svg",width = 10.00,height = 6.00)
Sum_plot <- plot_grid(as, cs, bs, ncol = 2, nrow = 2, rel_widths = c(3,3), rel_heights = c(3,3))
print(Sum_plot)
dev.off()  



    # Toute Condition Avril ---------------------------------------------------

taxo <- tableV4095 %>% select(OTU_Id,NN.taxonomy)
dataOTUTaxo04 <- merge(x = dataOTU04, y =  taxo, by = "OTU_Id")
dataOTUTaxo04 <- separate(dataOTUTaxo04, NN.taxonomy, c('Domain', 'Supergroup','Division'), sep = ";")
Parasitemin <- c("Apicomplexa","Perkinsea","Chytridiomycota","Cryptomycota","Ichthyophonida","Cercozoa","Oomycetes","Pirsonia","Labyrinthulida")
dataOTUTaxo04$Parasitism <- rep("Autres",each = nrow(dataOTUTaxo04))
w <- 0
for (j in dataOTUTaxo04$Division) { w <- w+1
for (i in Parasitemin) { if (i == j) { dataOTUTaxo04[w,"Parasitism"] <- "Potentiels parasites" }}}

dataOTUTaxo04 <-  dataOTUTaxo04 %>% filter(Parasitism == "Potentiels parasites")



aw <- ggplot(dataOTUTaxo04, aes(y = `Dim 2`, x = `Dim 1`, color = Cycles)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = Cycles, alpha = Cycles),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1OTU,y=Dim2OTU,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(aw)

#Figure FractionO
bw <- ggplot(dataOTUTaxo04, aes(y = `Dim 2`, x = `Dim 1`, color = FractionO)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = FractionO, alpha = FractionO),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1OTU,y=Dim2OTU,color = "Zones", alpha = "Zones", linetype = "Zones")
print(bw)

#Figure FractionT
cw <- ggplot(dataOTUTaxo04, aes(y = `Dim 2`, x = `Dim 1`, color = FractionT)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = FractionT, alpha = FractionT),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1OTU,y=Dim2OTU,color = "Fraction", alpha = "Fraction", linetype = "Fraction")
print(cw)

svglite("Analyse-Composition-Rarefy/AFC-Parasite/Analyse-AFC-OTU-Parasite-Avril-V4_095.svg",width = 10.00,height = 6.00)
b_plot <- plot_grid(aw, cw, bw, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(b_plot)
dev.off()  
    # Toute Condition Juin ---------------------------------------------------
taxo <- tableV4095 %>% select(OTU_Id,NN.taxonomy)
dataOTUTaxo06 <- merge(x = dataOTU06, y =  taxo, by = "OTU_Id")
dataOTUTaxo06 <- separate(dataOTUTaxo06, NN.taxonomy, c('Domain', 'Supergroup','Division'), sep = ";")
Parasitemin <- c("Apicomplexa","Perkinsea","Chytridiomycota","Cryptomycota","Ichthyophonida","Cercozoa","Oomycetes","Pirsonia","Labyrinthulida")
dataOTUTaxo06$Parasitism <- rep("Autres",each = nrow(dataOTUTaxo06))
w <- 0
for (j in dataOTUTaxo06$Division) { w <- w+1
for (i in Parasitemin) { if (i == j) { dataOTUTaxo06[w,"Parasitism"] <- "Potentiels parasites" }}}
dataOTUTaxo06 <-  dataOTUTaxo06 %>% filter(Parasitism == "Potentiels parasites")


aw <- ggplot(dataOTUTaxo06, aes(y = `Dim 2`, x = `Dim 1`, color = Cycles)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = Cycles, alpha = Cycles),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1OTU,y=Dim2OTU,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(aw)

#Figure FractionO
bw <- ggplot(dataOTUTaxo06, aes(y = `Dim 2`, x = `Dim 1`, color = FractionO)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = FractionO, alpha = FractionO),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1OTU,y=Dim2OTU,color = "Zones", alpha = "Zones", linetype = "Zones")
print(bw)

#Figure FractionT
cw <- ggplot(dataOTUTaxo06, aes(y = `Dim 2`, x = `Dim 1`, color = FractionT)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = FractionT, alpha = FractionT),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1OTU,y=Dim2OTU,color = "Fraction", alpha = "Fraction", linetype = "Fraction")
print(cw)

#Potentiel parasite


svglite("Analyse-Composition-Rarefy/AFC-Parasite/Analyse-AFC-OTU-Parasite-Juin-V4_095.svg",width = 10.00,height = 6.00)
b_plot <- plot_grid(aw, cw, bw, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(b_plot)
dev.off()  
    # Toute Condition Septembre ---------------------------------------------------
taxo <- tableV4095 %>% select(OTU_Id,NN.taxonomy)
dataOTUTaxo09 <- merge(x = dataOTU09, y =  taxo, by = "OTU_Id")
dataOTUTaxo09 <- separate(dataOTUTaxo09, NN.taxonomy, c('Domain', 'Supergroup','Division'), sep = ";")
Parasitemin <- c("Apicomplexa","Perkinsea","Chytridiomycota","Cryptomycota","Ichthyophonida","Cercozoa","Oomycetes","Pirsonia","Labyrinthulida")
dataOTUTaxo09$Parasitism <- rep("Autres",each = nrow(dataOTUTaxo09))
w <- 0
for (j in dataOTUTaxo09$Division) { w <- w+1
for (i in Parasitemin) { if (i == j) { dataOTUTaxo09[w,"Parasitism"] <- "Potentiels parasites" }}}
dataOTUTaxo09 <-  dataOTUTaxo09 %>% filter(Parasitism == "Potentiels parasites")



aw <- ggplot(dataOTUTaxo09, aes(y = `Dim 2`, x = `Dim 1`, color = Cycles)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = Cycles, alpha = Cycles),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1OTU,y=Dim2OTU,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(aw)

#Figure FractionO
bw <- ggplot(dataOTUTaxo09, aes(y = `Dim 2`, x = `Dim 1`, color = FractionO)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = FractionO, alpha = FractionO),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1OTU,y=Dim2OTU,color = "Zones", alpha = "Zones", linetype = "Zones")
print(bw)

#Figure FractionT
cw <- ggplot(dataOTUTaxo09, aes(y = `Dim 2`, x = `Dim 1`, color = FractionT)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = FractionT, alpha = FractionT),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1OTU,y=Dim2OTU,color = "Fraction", alpha = "Fraction", linetype = "Fraction")
print(cw)

#Potentiel parasite

svglite("Analyse-Composition-Rarefy/AFC-Parasite/Analyse-AFC-OTU-Parasite-Septembre-V4_095.svg",width = 10.00,height = 6.00)
b_plot <- plot_grid(aw, cw, bw, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(b_plot)
dev.off()  
    # Toute Condition Novembre ---------------------------------------------------
taxo <- tableV4095 %>% select(OTU_Id,NN.taxonomy)
dataOTUTaxo11 <- merge(x = dataOTU11, y =  taxo, by = "OTU_Id")
dataOTUTaxo11 <- separate(dataOTUTaxo11, NN.taxonomy, c('Domain', 'Supergroup','Division'), sep = ";")
Parasitemin <- c("Apicomplexa","Perkinsea","Chytridiomycota","Cryptomycota","Ichthyophonida","Cercozoa","Oomycetes","Pirsonia","Labyrinthulida")
dataOTUTaxo11$Parasitism <- rep("Autres",each = nrow(dataOTUTaxo11))
w <- 0
for (j in dataOTUTaxo11$Division) { w <- w+1
for (i in Parasitemin) { if (i == j) { dataOTUTaxo11[w,"Parasitism"] <- "Potentiels parasites" }}}
dataOTUTaxo11 <-  dataOTUTaxo11 %>% filter(Parasitism == "Potentiels parasites")

aw <- ggplot(dataOTUTaxo11, aes(y = `Dim 2`, x = `Dim 1`, color = Cycles)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = Cycles, alpha = Cycles),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1OTU,y=Dim2OTU,color = "Cycles", alpha = "Cycles", linetype = "Cycles")
print(aw)

#Figure FractionO
bw <- ggplot(dataOTUTaxo11, aes(y = `Dim 2`, x = `Dim 1`, color = FractionO)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = FractionO, alpha = FractionO),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(2,0,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0.1,0,0.1)) +
  scale_color_manual(values = c("#F8766D","lightgrey","#00A5FF")) +
  labs(x=Dim1OTU,y=Dim2OTU,color = "Zones", alpha = "Zones", linetype = "Zones")
print(bw)

#Figure FractionT
cw <- ggplot(dataOTUTaxo11, aes(y = `Dim 2`, x = `Dim 1`, color = FractionT)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
  geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
  geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
  #facet_grid(~`Dates`, switch = "x") + theme_bw() +
  #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
  stat_ellipse(aes(linetype = FractionT, alpha = FractionT),geom = "polygon",type = "norm") +
  scale_linetype_manual(values = c(0,2,2)) +
  #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
  #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
  #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
  theme(axis.title = element_text(face="bold", size=12), 
        axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
        title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10)) +
  scale_alpha_manual(values = c(0,0.1,0.1)) +
  scale_color_manual(values = c("lightgrey","#F8766D","#00A5FF")) +
  labs(x=Dim1OTU,y=Dim2OTU,color = "Fraction", alpha = "Fraction", linetype = "Fraction")
print(cw)

#Potentiel parasite

svglite("Analyse-Composition-Rarefy/AFC-Parasite/Analyse-AFC-OTU-Parasite-Novembre-V4_095.svg",width = 10.00,height = 6.00)
b_plot <- plot_grid(aw, cw, bw, ncol = 2, nrow = 2, rel_widths = c(3,3),rel_heights = c(3,3))
print(b_plot)
dev.off()  


# HIST  --------------------------------------------------------------------
  # HIST OTU ONLY---------------------------------------------------------------------
    # Cycles only -------------------------------------------------------------------
  onlyJourOTU <- dataOTUTaxo %>% filter(Cycles == "Jour") %>% select(OTU_Id,TotalJour,Division)
  row.names(onlyJourOTU)<-onlyJourOTU$OTU_Id
  onlyJourOTU <- onlyJourOTU %>% select(-OTU_Id)
  onlyJourOTU <- onlyJourOTU %>% group_by(Division) %>% summarise_all(funs(sum))
  onlyJourOTU$TotalJour <- onlyJourOTU$TotalJour*100/sum(onlyJourOTU$TotalJour)
  onlyJourOTU <- onlyJourOTU %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(TotalJour) - 0.5*TotalJour)
  onlyJourOTU$label <- paste(round(onlyJourOTU$TotalJour,1), "%", sep = "")
  for (i in rownames(onlyJourOTU)) {
    if (onlyJourOTU[i,"label"] == "0%") { onlyJourOTU[i,"label"] <- NA}}
  for (i in rownames(onlyJourOTU)) {
    if (is.na(onlyJourOTU[i,"label"]) == FALSE) { onlyJourOTU[i,"label"] <- paste(onlyJourOTU[i,"Division"]," : ",onlyJourOTU[i,"label"], sep = "")}}
  onlyJourOTU$Cycles<- rep("Jour", each = nrow(onlyJourOTU))
  ##
  onlyNuitOTU <- dataOTUTaxo %>% filter(Cycles == "Nuit") %>% select(OTU_Id,TotalNuit,Division)
  row.names(onlyNuitOTU)<-onlyNuitOTU$OTU_Id
  onlyNuitOTU <- onlyNuitOTU %>% select(-OTU_Id)
  onlyNuitOTU <- onlyNuitOTU %>% group_by(Division) %>% summarise_all(funs(sum))
  onlyNuitOTU$TotalNuit <- onlyNuitOTU$TotalNuit*100/sum(onlyNuitOTU$TotalNuit)
  onlyNuitOTU <- onlyNuitOTU %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(TotalNuit) - 0.5*TotalNuit)
  onlyNuitOTU$label <- paste(round(onlyNuitOTU$TotalNuit,1), "%", sep = "")
  for (i in rownames(onlyNuitOTU)) {
    if (onlyNuitOTU[i,"label"] == "0%") { onlyNuitOTU[i,"label"] <- NA}}
  for (i in rownames(onlyNuitOTU)) {
    if (is.na(onlyNuitOTU[i,"label"]) == FALSE) { onlyNuitOTU[i,"label"] <- paste(onlyNuitOTU[i,"Division"]," : ",onlyNuitOTU[i,"label"], sep = "")}}
  onlyNuitOTU$Cycles<- rep("Nuit", each = nrow(onlyNuitOTU))
  ##
  colnames(onlyJourOTU)[2]  <- "value"
  onlyJourOTU$Sum <- rep(0, each = nrow(onlyJourOTU))
  onlyJourOTU$Count <- rep(0, each = nrow(onlyJourOTU))
  for (i in onlyJourOTU$Division) { onlyJourOTU$Count[which(onlyJourOTU$Division == i)] <- nrow(dataOTUTaxo %>% filter(Cycles == "Jour") %>% filter(Division == i))}
  onlyJourOTU$Sum <- sum(onlyJourOTU$Count)
  colnames(onlyNuitOTU)[2]  <- "value"
  onlyNuitOTU$Sum <- rep(0, each = nrow(onlyNuitOTU))
  onlyNuitOTU$Count <- rep(0, each = nrow(onlyNuitOTU))
  for (i in onlyNuitOTU$Division) { onlyNuitOTU$Count[which(onlyNuitOTU$Division == i)] <- nrow(dataOTUTaxo %>% filter(Cycles == "Nuit") %>% filter(Division == i))}
  onlyNuitOTU$Sum <- sum(onlyNuitOTU$Count)
  onlyCyclesOTU <- rbind(onlyJourOTU,onlyNuitOTU)
  #Figure
  az <- ggplot(onlyCyclesOTU, mapping = aes(y= value, x = Cycles, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
    geom_label(aes(y = 104,label = paste(Sum," OTUs",sep ="")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
    scale_fill_manual(values = palette)
    #scale_fill_locuszoom()
  legendOTU <- get_legend(az)
  az <- az + theme(legend.position = "none") + labs(x="Cycles",y="OTUs (%)")
  print(az)
  
  
    # FractionO only -------------------------------------------------------------------
  onlyOxiqueOTU <- dataOTUTaxo %>% filter(FractionO == "Oxique") %>% select(OTU_Id,TotalOxique,Division)
  row.names(onlyOxiqueOTU)<-onlyOxiqueOTU$OTU_Id
  onlyOxiqueOTU <- onlyOxiqueOTU %>% select(-OTU_Id)
  onlyOxiqueOTU <- onlyOxiqueOTU %>% group_by(Division) %>% summarise_all(funs(sum))
  onlyOxiqueOTU$TotalOxique <- onlyOxiqueOTU$TotalOxique*100/sum(onlyOxiqueOTU$TotalOxique)
  onlyOxiqueOTU <- onlyOxiqueOTU %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(TotalOxique) - 0.5*TotalOxique)
  onlyOxiqueOTU$label <- paste(round(onlyOxiqueOTU$TotalOxique,1), "%", sep = "")
  for (i in rownames(onlyOxiqueOTU)) {
    if (onlyOxiqueOTU[i,"label"] == "0%") { onlyOxiqueOTU[i,"label"] <- NA}}
  for (i in rownames(onlyOxiqueOTU)) {
    if (is.na(onlyOxiqueOTU[i,"label"]) == FALSE) { onlyOxiqueOTU[i,"label"] <- paste(onlyOxiqueOTU[i,"Division"]," : ",onlyOxiqueOTU[i,"label"], sep = "")}}
  onlyOxiqueOTU$FractionO<- rep("Oxique", each = nrow(onlyOxiqueOTU))
  ##
  onlyAnoxiqueOTU <- dataOTUTaxo %>% filter(FractionO == "Anoxique") %>% select(OTU_Id,TotalAnoxique,Division)
  row.names(onlyAnoxiqueOTU)<-onlyAnoxiqueOTU$OTU_Id
  onlyAnoxiqueOTU <- onlyAnoxiqueOTU %>% select(-OTU_Id)
  onlyAnoxiqueOTU <- onlyAnoxiqueOTU %>% group_by(Division) %>% summarise_all(funs(sum))
  onlyAnoxiqueOTU$TotalAnoxique <- onlyAnoxiqueOTU$TotalAnoxique*100/sum(onlyAnoxiqueOTU$TotalAnoxique)
  onlyAnoxiqueOTU <- onlyAnoxiqueOTU %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(TotalAnoxique) - 0.5*TotalAnoxique)
  onlyAnoxiqueOTU$label <- paste(round(onlyAnoxiqueOTU$TotalAnoxique,1), "%", sep = "")
  for (i in rownames(onlyAnoxiqueOTU)) {
    if (onlyAnoxiqueOTU[i,"label"] == "0%") { onlyAnoxiqueOTU[i,"label"] <- NA}}
  for (i in rownames(onlyAnoxiqueOTU)) {
    if (is.na(onlyAnoxiqueOTU[i,"label"]) == FALSE) { onlyAnoxiqueOTU[i,"label"] <- paste(onlyAnoxiqueOTU[i,"Division"]," : ",onlyAnoxiqueOTU[i,"label"], sep = "")}}
  onlyAnoxiqueOTU$FractionO<- rep("Anoxique", each = nrow(onlyAnoxiqueOTU))
  ##
  colnames(onlyOxiqueOTU)[2]  <- "value"
  onlyOxiqueOTU$Sum <- rep(0, each = nrow(onlyOxiqueOTU))
  onlyOxiqueOTU$Count <- rep(0, each = nrow(onlyOxiqueOTU))
  for (i in onlyOxiqueOTU$Division) { onlyOxiqueOTU$Count[which(onlyOxiqueOTU$Division == i)] <- nrow(dataOTUTaxo %>% filter(FractionO == "Oxique") %>% filter(Division == i))}
  onlyOxiqueOTU$Sum <- sum(onlyOxiqueOTU$Count)
  colnames(onlyAnoxiqueOTU)[2]  <- "value"
  onlyAnoxiqueOTU$Sum <- rep(0, each = nrow(onlyAnoxiqueOTU))
  onlyAnoxiqueOTU$Count <- rep(0, each = nrow(onlyAnoxiqueOTU))
  for (i in onlyAnoxiqueOTU$Division) { onlyAnoxiqueOTU$Count[which(onlyAnoxiqueOTU$Division == i)] <- nrow(dataOTUTaxo %>% filter(FractionO == "Anoxique") %>% filter(Division == i))}
  onlyAnoxiqueOTU$Sum <- sum(onlyAnoxiqueOTU$Count)
  onlyFractionOOTU <- rbind(onlyOxiqueOTU,onlyAnoxiqueOTU)
  #Figure
  ax <- ggplot(onlyFractionOOTU, mapping = aes(y= value, x = FractionO, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
    geom_label(aes(y = 104,label = paste(Sum," OTUs",sep ="")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
    scale_fill_manual(values = palette)
  ax <- ax + theme(legend.position = "none") + labs(x="Zones",y="OTUs (%)")
  print(ax)    
    # FractionT only -------------------------------------------------------------------
  onlyPetiteOTU <- dataOTUTaxo %>% filter(FractionT == "Petite") %>% select(OTU_Id,TotalPetite,Division)
  row.names(onlyPetiteOTU)<-onlyPetiteOTU$OTU_Id
  onlyPetiteOTU <- onlyPetiteOTU %>% select(-OTU_Id)
  onlyPetiteOTU <- onlyPetiteOTU %>% group_by(Division) %>% summarise_all(funs(sum))
  onlyPetiteOTU$TotalPetite <- onlyPetiteOTU$TotalPetite*100/sum(onlyPetiteOTU$TotalPetite)
  onlyPetiteOTU <- onlyPetiteOTU %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(TotalPetite) - 0.5*TotalPetite)
  onlyPetiteOTU$label <- paste(round(onlyPetiteOTU$TotalPetite,1), "%", sep = "")
  for (i in rownames(onlyPetiteOTU)) {
    if (onlyPetiteOTU[i,"label"] == "0%") { onlyPetiteOTU[i,"label"] <- NA}}
  for (i in rownames(onlyPetiteOTU)) {
    if (is.na(onlyPetiteOTU[i,"label"]) == FALSE) { onlyPetiteOTU[i,"label"] <- paste(onlyPetiteOTU[i,"Division"]," : ",onlyPetiteOTU[i,"label"], sep = "")}}
  onlyPetiteOTU$FractionT<- rep("Petite", each = nrow(onlyPetiteOTU))
  ##
  onlyGrandeOTU <- dataOTUTaxo %>% filter(FractionT == "Grande") %>% select(OTU_Id,TotalGrande,Division)
  row.names(onlyGrandeOTU)<-onlyGrandeOTU$OTU_Id
  onlyGrandeOTU <- onlyGrandeOTU %>% select(-OTU_Id)
  onlyGrandeOTU <- onlyGrandeOTU %>% group_by(Division) %>% summarise_all(funs(sum))
  onlyGrandeOTU$TotalGrande <- onlyGrandeOTU$TotalGrande*100/sum(onlyGrandeOTU$TotalGrande)
  onlyGrandeOTU <- onlyGrandeOTU %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(TotalGrande) - 0.5*TotalGrande)
  onlyGrandeOTU$label <- paste(round(onlyGrandeOTU$TotalGrande,1), "%", sep = "")
  for (i in rownames(onlyGrandeOTU)) {
    if (onlyGrandeOTU[i,"label"] == "0%") { onlyGrandeOTU[i,"label"] <- NA}}
  for (i in rownames(onlyGrandeOTU)) {
    if (is.na(onlyGrandeOTU[i,"label"]) == FALSE) { onlyGrandeOTU[i,"label"] <- paste(onlyGrandeOTU[i,"Division"]," : ",onlyGrandeOTU[i,"label"], sep = "")}}
  onlyGrandeOTU$FractionT<- rep("Grande", each = nrow(onlyGrandeOTU))
  ##
  colnames(onlyPetiteOTU)[2]  <- "value"
  onlyPetiteOTU$Sum <- rep(0, each = nrow(onlyPetiteOTU))
  onlyPetiteOTU$Count <- rep(0, each = nrow(onlyPetiteOTU))
  for (i in onlyPetiteOTU$Division) { onlyPetiteOTU$Count[which(onlyPetiteOTU$Division == i)] <- nrow(dataOTUTaxo %>% filter(FractionT == "Petite") %>% filter(Division == i))}
  onlyPetiteOTU$Sum <- sum(onlyPetiteOTU$Count)
  colnames(onlyGrandeOTU)[2]  <- "value"
  onlyGrandeOTU$Sum <- rep(0, each = nrow(onlyGrandeOTU))
  onlyGrandeOTU$Count <- rep(0, each = nrow(onlyGrandeOTU))
  for (i in onlyGrandeOTU$Division) { onlyGrandeOTU$Count[which(onlyGrandeOTU$Division == i)] <- nrow(dataOTUTaxo %>% filter(FractionT == "Grande") %>% filter(Division == i))}
  onlyGrandeOTU$Sum <- sum(onlyGrandeOTU$Count)
  onlyFractionTOTU <- rbind(onlyPetiteOTU,onlyGrandeOTU)
  #Figure
  ay <- ggplot(onlyFractionTOTU, mapping = aes(y= value, x = FractionT, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
    geom_label(aes(y = 104,label = paste(Sum," OTUs",sep ="")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
    scale_fill_manual(values = palette)
  ay <- ay + theme(legend.position = "none") + labs(x="Fractions",y="OTUs (%)")
  print(ay) 
    # Coplot -------------------------------------------------------------------
  svglite("Analyse-Composition-Rarefy/HistOnly/Analyse-OTU-Parasite-only-V4_095.svg",width = 8.00,height = 6.00)
  b_plot <- plot_grid(az,ax,ay,legendOTU,  ncol = 4, nrow = 1, rel_widths = c(3,3,3,2),rel_heights = c(3))
  print(b_plot)
  dev.off()  
    # Coplot X AFC -------------------------------------------------------------------
  f <- ggplot() + theme_void()
  print(f)
  
  svglite("Analyse-Composition-Rarefy/Biplot/Analyse-OTU-Parasite-onlyX-V4_095.svg",width = 10.00,height = 13.00)
  b_plot <- plot_grid(as,az,f,bs,ax,legendOTU,cs,ay,f, ncol = 3, nrow = 3, rel_widths = c(3,1.5,1),rel_heights = c(3,3,3))
  print(b_plot)
  dev.off()  
  svglite("Analyse-Composition-Rarefy/Biplot/Analyse-SequenceXOTU-Parasite-onlyX-V4_095.svg",width = 10.00,height = 13.00)
  b_plot <- plot_grid(a,az,f,b,ax,legendOTU,c,ay,f, ncol = 3, nrow = 3, rel_widths = c(3,1.5,1),rel_heights = c(3,3,3))
  print(b_plot)
  dev.off()  
  
  
  
  # HIST OTU TOTAL---------------------------------------------------------------------
    # Cycles -------------------------------------------------------------------
  totalJourOTU <- dataOTUTaxo %>% select(OTU_Id,TotalJour,Division)
  row.names(totalJourOTU)<-totalJourOTU$OTU_Id
  totalJourOTU <- totalJourOTU %>% select(-OTU_Id)
  totalJourOTU <- totalJourOTU %>% group_by(Division) %>% summarise_all(funs(sum))
  totalJourOTU$TotalJour <- totalJourOTU$TotalJour*100/sum(totalJourOTU$TotalJour)
  totalJourOTU <- totalJourOTU %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(TotalJour) - 0.5*TotalJour)
  totalJourOTU$label <- paste(round(totalJourOTU$TotalJour,1), "%", sep = "")
  for (i in rownames(totalJourOTU)) {
    if (totalJourOTU[i,"label"] == "0%") { totalJourOTU[i,"label"] <- NA}}
  for (i in rownames(totalJourOTU)) {
    if (is.na(totalJourOTU[i,"label"]) == FALSE) { totalJourOTU[i,"label"] <- paste(totalJourOTU[i,"Division"]," : ",totalJourOTU[i,"label"], sep = "")}}
  totalJourOTU$Cycles<- rep("Jour", each = nrow(totalJourOTU))
  ##
  totalNuitOTU <- dataOTUTaxo %>% select(OTU_Id,TotalNuit,Division)
  row.names(totalNuitOTU)<-totalNuitOTU$OTU_Id
  totalNuitOTU <- totalNuitOTU %>% select(-OTU_Id)
  totalNuitOTU <- totalNuitOTU %>% group_by(Division) %>% summarise_all(funs(sum))
  totalNuitOTU$TotalNuit <- totalNuitOTU$TotalNuit*100/sum(totalNuitOTU$TotalNuit)
  totalNuitOTU <- totalNuitOTU %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(TotalNuit) - 0.5*TotalNuit)
  totalNuitOTU$label <- paste(round(totalNuitOTU$TotalNuit,1), "%", sep = "")
  for (i in rownames(totalNuitOTU)) {
    if (totalNuitOTU[i,"label"] == "0%") { totalNuitOTU[i,"label"] <- NA}}
  for (i in rownames(totalNuitOTU)) {
    if (is.na(totalNuitOTU[i,"label"]) == FALSE) { totalNuitOTU[i,"label"] <- paste(totalNuitOTU[i,"Division"]," : ",totalNuitOTU[i,"label"], sep = "")}}
  totalNuitOTU$Cycles<- rep("Nuit", each = nrow(totalNuitOTU))
  ##
  colnames(totalJourOTU)[2]  <- "value"
  totalJourOTU$Sum <- rep(0, each = nrow(totalJourOTU))
  totalJourOTU$Count <- rep(0, each = nrow(totalJourOTU))
  for (i in totalJourOTU$Division) { totalJourOTU$Count[which(totalJourOTU$Division == i)] <- sum(dataOTUTaxo  %>% filter(Division == i) %>% select(TotalJour))}
  totalJourOTU$Sum <- sum(totalJourOTU$Count)
  colnames(totalNuitOTU)[2]  <- "value"
  totalNuitOTU$Sum <- rep(0, each = nrow(totalNuitOTU))
  totalNuitOTU$Count <- rep(0, each = nrow(totalNuitOTU))
  for (i in totalNuitOTU$Division) { totalNuitOTU$Count[which(totalNuitOTU$Division == i)] <- sum(dataOTUTaxo  %>% filter(Division == i) %>% select(TotalNuit))}
  totalNuitOTU$Sum <- sum(totalNuitOTU$Count)
  totalCyclesOTU <- rbind(totalJourOTU,totalNuitOTU)

  #Figure
  bz <- ggplot(totalCyclesOTU, mapping = aes(y= value, x = Cycles, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + scale_fill_manual(values = palette) + 
    geom_label(aes(y = 104,label = paste(Sum," OTUs",sep ="")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
  legendOTU <- get_legend(bz)
  bz <- bz + theme(legend.position = "none") + labs(x="Cycles",y="OTUs (%)")
  print(bz)
  
    # FractionO -------------------------------------------------------------------
  totalOxiqueOTU <- dataOTUTaxo %>% select(OTU_Id,TotalOxique,Division)
  row.names(totalOxiqueOTU)<-totalOxiqueOTU$OTU_Id
  totalOxiqueOTU <- totalOxiqueOTU %>% select(-OTU_Id)
  totalOxiqueOTU <- totalOxiqueOTU %>% group_by(Division) %>% summarise_all(funs(sum))
  totalOxiqueOTU$TotalOxique <- totalOxiqueOTU$TotalOxique*100/sum(totalOxiqueOTU$TotalOxique)
  totalOxiqueOTU <- totalOxiqueOTU %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(TotalOxique) - 0.5*TotalOxique)
  totalOxiqueOTU$label <- paste(round(totalOxiqueOTU$TotalOxique,1), "%", sep = "")
  for (i in rownames(totalOxiqueOTU)) {
    if (totalOxiqueOTU[i,"label"] == "0%") { totalOxiqueOTU[i,"label"] <- NA}}
  for (i in rownames(totalOxiqueOTU)) {
    if (is.na(totalOxiqueOTU[i,"label"]) == FALSE) { totalOxiqueOTU[i,"label"] <- paste(totalOxiqueOTU[i,"Division"]," : ",totalOxiqueOTU[i,"label"], sep = "")}}
  totalOxiqueOTU$FractionO<- rep("Oxique", each = nrow(totalOxiqueOTU))
  ##
  totalAnoxiqueOTU <- dataOTUTaxo %>% select(OTU_Id,TotalAnoxique,Division)
  row.names(totalAnoxiqueOTU)<-totalAnoxiqueOTU$OTU_Id
  totalAnoxiqueOTU <- totalAnoxiqueOTU %>% select(-OTU_Id)
  totalAnoxiqueOTU <- totalAnoxiqueOTU %>% group_by(Division) %>% summarise_all(funs(sum))
  totalAnoxiqueOTU$TotalAnoxique <- totalAnoxiqueOTU$TotalAnoxique*100/sum(totalAnoxiqueOTU$TotalAnoxique)
  totalAnoxiqueOTU <- totalAnoxiqueOTU %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(TotalAnoxique) - 0.5*TotalAnoxique)
  totalAnoxiqueOTU$label <- paste(round(totalAnoxiqueOTU$TotalAnoxique,1), "%", sep = "")
  for (i in rownames(totalAnoxiqueOTU)) {
    if (totalAnoxiqueOTU[i,"label"] == "0%") { totalAnoxiqueOTU[i,"label"] <- NA}}
  for (i in rownames(totalAnoxiqueOTU)) {
    if (is.na(totalAnoxiqueOTU[i,"label"]) == FALSE) { totalAnoxiqueOTU[i,"label"] <- paste(totalAnoxiqueOTU[i,"Division"]," : ",totalAnoxiqueOTU[i,"label"], sep = "")}}
  totalAnoxiqueOTU$FractionO<- rep("Anoxique", each = nrow(totalAnoxiqueOTU))
  ##
  colnames(totalOxiqueOTU)[2]  <- "value"
  totalOxiqueOTU$Sum <- rep(0, each = nrow(totalOxiqueOTU))
  totalOxiqueOTU$Count <- rep(0, each = nrow(totalOxiqueOTU))
  for (i in totalOxiqueOTU$Division) { totalOxiqueOTU$Count[which(totalOxiqueOTU$Division == i)] <- sum(dataOTUTaxo  %>% filter(Division == i) %>% select(TotalOxique))}
  totalOxiqueOTU$Sum <- sum(totalOxiqueOTU$Count)
  colnames(totalAnoxiqueOTU)[2]  <- "value"
  totalAnoxiqueOTU$Sum <- rep(0, each = nrow(totalAnoxiqueOTU))
  totalAnoxiqueOTU$Count <- rep(0, each = nrow(totalAnoxiqueOTU))
  for (i in totalAnoxiqueOTU$Division) { totalAnoxiqueOTU$Count[which(totalAnoxiqueOTU$Division == i)] <- sum(dataOTUTaxo  %>% filter(Division == i) %>% select(TotalAnoxique))}
  totalAnoxiqueOTU$Sum <- sum(totalAnoxiqueOTU$Count)
  totalFractionOOTU <- rbind(totalOxiqueOTU,totalAnoxiqueOTU)
  #Figure
  bx <- ggplot(totalFractionOOTU, mapping = aes(y= value, x = FractionO, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + scale_fill_manual(values = palette) + 
    geom_label(aes(y = 104,label = paste(Sum," OTUs",sep ="")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
  bx <- bx + theme(legend.position = "none") + labs(x="Zones",y="OTUs (%)")
  print(bx)    
    # FractionT -------------------------------------------------------------------
  totalPetiteOTU <- dataOTUTaxo %>% select(OTU_Id,TotalPetite,Division)
  row.names(totalPetiteOTU)<-totalPetiteOTU$OTU_Id
  totalPetiteOTU <- totalPetiteOTU %>% select(-OTU_Id)
  totalPetiteOTU <- totalPetiteOTU %>% group_by(Division) %>% summarise_all(funs(sum))
  totalPetiteOTU$TotalPetite <- totalPetiteOTU$TotalPetite*100/sum(totalPetiteOTU$TotalPetite)
  totalPetiteOTU <- totalPetiteOTU %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(TotalPetite) - 0.5*TotalPetite)
  totalPetiteOTU$label <- paste(round(totalPetiteOTU$TotalPetite,1), "%", sep = "")
  for (i in rownames(totalPetiteOTU)) {
    if (totalPetiteOTU[i,"label"] == "0%") { totalPetiteOTU[i,"label"] <- NA}}
  for (i in rownames(totalPetiteOTU)) {
    if (is.na(totalPetiteOTU[i,"label"]) == FALSE) { totalPetiteOTU[i,"label"] <- paste(totalPetiteOTU[i,"Division"]," : ",totalPetiteOTU[i,"label"], sep = "")}}
  totalPetiteOTU$FractionT<- rep("Petite", each = nrow(totalPetiteOTU))
  ##
  totalGrandeOTU <- dataOTUTaxo %>% select(OTU_Id,TotalGrande,Division)
  row.names(totalGrandeOTU)<-totalGrandeOTU$OTU_Id
  totalGrandeOTU <- totalGrandeOTU %>% select(-OTU_Id)
  totalGrandeOTU <- totalGrandeOTU %>% group_by(Division) %>% summarise_all(funs(sum))
  totalGrandeOTU$TotalGrande <- totalGrandeOTU$TotalGrande*100/sum(totalGrandeOTU$TotalGrande)
  totalGrandeOTU <- totalGrandeOTU %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(TotalGrande) - 0.5*TotalGrande)
  totalGrandeOTU$label <- paste(round(totalGrandeOTU$TotalGrande,1), "%", sep = "")
  for (i in rownames(totalGrandeOTU)) {
    if (totalGrandeOTU[i,"label"] == "0%") { totalGrandeOTU[i,"label"] <- NA}}
  for (i in rownames(totalGrandeOTU)) {
    if (is.na(totalGrandeOTU[i,"label"]) == FALSE) { totalGrandeOTU[i,"label"] <- paste(totalGrandeOTU[i,"Division"]," : ",totalGrandeOTU[i,"label"], sep = "")}}
  totalGrandeOTU$FractionT<- rep("Grande", each = nrow(totalGrandeOTU))
  ##
  colnames(totalPetiteOTU)[2]  <- "value"
  totalPetiteOTU$Sum <- rep(0, each = nrow(totalPetiteOTU))
  totalPetiteOTU$Count <- rep(0, each = nrow(totalPetiteOTU))
  for (i in totalPetiteOTU$Division) { totalPetiteOTU$Count[which(totalPetiteOTU$Division == i)] <- sum(dataOTUTaxo  %>% filter(Division == i) %>% select(TotalPetite))}
  totalPetiteOTU$Sum <- sum(totalPetiteOTU$Count)
  colnames(totalGrandeOTU)[2]  <- "value"
  totalGrandeOTU$Sum <- rep(0, each = nrow(totalGrandeOTU))
  totalGrandeOTU$Count <- rep(0, each = nrow(totalGrandeOTU))
  for (i in totalGrandeOTU$Division) { totalGrandeOTU$Count[which(totalGrandeOTU$Division == i)] <- sum(dataOTUTaxo  %>% filter(Division == i) %>% select(TotalGrande))}
  totalGrandeOTU$Sum <- sum(totalGrandeOTU$Count)
  totalFractionTOTU <- rbind(totalPetiteOTU,totalGrandeOTU)
  #Figure
  by <- ggplot(totalFractionTOTU, mapping = aes(y= value, x = FractionT, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + scale_fill_manual(values = palette)+ 
    geom_label(aes(y = 104,label = paste(Sum," OTUs",sep ="")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
  by <- by + theme(legend.position = "none") + labs(x="Fractions",y="OTUs (%)")
  print(by) 
    # Dates -------------------------------------------------------------------
  total04OTU <- dataOTUTaxo %>% select(OTU_Id,Total04,Division)
  row.names(total04OTU)<-total04OTU$OTU_Id
  total04OTU <- total04OTU %>% select(-OTU_Id)
  total04OTU <- total04OTU %>% group_by(Division) %>% summarise_all(funs(sum))
  total04OTU$Total04 <- total04OTU$Total04*100/sum(total04OTU$Total04)
  total04OTU <- total04OTU %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(Total04) - 0.5*Total04)
  total04OTU$label <- paste(round(total04OTU$Total04,1), "%", sep = "")
  for (i in rownames(total04OTU)) {
    if (total04OTU[i,"label"] == "0%") { total04OTU[i,"label"] <- NA}}
  for (i in rownames(total04OTU)) {
    if (is.na(total04OTU[i,"label"]) == FALSE) { total04OTU[i,"label"] <- paste(total04OTU[i,"Division"]," : ",total04OTU[i,"label"], sep = "")}}
  total04OTU$Dates<- rep("04", each = nrow(total04OTU))
  ##
  total06OTU <- dataOTUTaxo %>% select(OTU_Id,Total06,Division)
  row.names(total06OTU)<-total06OTU$OTU_Id
  total06OTU <- total06OTU %>% select(-OTU_Id)
  total06OTU <- total06OTU %>% group_by(Division) %>% summarise_all(funs(sum))
  total06OTU$Total06 <- total06OTU$Total06*100/sum(total06OTU$Total06)
  total06OTU <- total06OTU %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(Total06) - 0.5*Total06)
  total06OTU$label <- paste(round(total06OTU$Total06,1), "%", sep = "")
  for (i in rownames(total06OTU)) {
    if (total06OTU[i,"label"] == "0%") { total06OTU[i,"label"] <- NA}}
  for (i in rownames(total06OTU)) {
    if (is.na(total06OTU[i,"label"]) == FALSE) { total06OTU[i,"label"] <- paste(total06OTU[i,"Division"]," : ",total06OTU[i,"label"], sep = "")}}
  total06OTU$Dates<- rep("06", each = nrow(total06OTU))
  ##
  total09OTU <- dataOTUTaxo %>% select(OTU_Id,Total09,Division)
  row.names(total09OTU)<-total09OTU$OTU_Id
  total09OTU <- total09OTU %>% select(-OTU_Id)
  total09OTU <- total09OTU %>% group_by(Division) %>% summarise_all(funs(sum))
  total09OTU$Total09 <- total09OTU$Total09*100/sum(total09OTU$Total09)
  total09OTU <- total09OTU %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(Total09) - 0.5*Total09)
  total09OTU$label <- paste(round(total09OTU$Total09,1), "%", sep = "")
  for (i in rownames(total09OTU)) {
    if (total09OTU[i,"label"] == "0%") { total09OTU[i,"label"] <- NA}}
  for (i in rownames(total09OTU)) {
    if (is.na(total09OTU[i,"label"]) == FALSE) { total09OTU[i,"label"] <- paste(total09OTU[i,"Division"]," : ",total09OTU[i,"label"], sep = "")}}
  total09OTU$Dates<- rep("09", each = nrow(total09OTU))
  ##
  total11OTU <- dataOTUTaxo %>% select(OTU_Id,Total11,Division)
  row.names(total11OTU)<-total11OTU$OTU_Id
  total11OTU <- total11OTU %>% select(-OTU_Id)
  total11OTU <- total11OTU %>% group_by(Division) %>% summarise_all(funs(sum))
  total11OTU$Total11 <- total11OTU$Total11*100/sum(total11OTU$Total11)
  total11OTU <- total11OTU %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(Total11) - 0.5*Total11)
  total11OTU$label <- paste(round(total11OTU$Total11,1), "%", sep = "")
  for (i in rownames(total11OTU)) {
    if (total11OTU[i,"label"] == "0%") { total11OTU[i,"label"] <- NA}}
  for (i in rownames(total11OTU)) {
    if (is.na(total11OTU[i,"label"]) == FALSE) { total11OTU[i,"label"] <- paste(total11OTU[i,"Division"]," : ",total11OTU[i,"label"], sep = "")}}
  total11OTU$Dates<- rep("11", each = nrow(total11OTU))
  ##
  colnames(total04OTU)[2]  <- "value"
  total04OTU$Sum <- rep(0, each = nrow(total04OTU))
  total04OTU$Count <- rep(0, each = nrow(total04OTU))
  for (i in total04OTU$Division) { total04OTU$Count[which(total04OTU$Division == i)] <- sum(dataOTUTaxo  %>% filter(Division == i) %>% select(Total04))}
  total04OTU$Sum <- sum(total04OTU$Count)
  colnames(total06OTU)[2]  <- "value"
  total06OTU$Sum <- rep(0, each = nrow(total06OTU))
  total06OTU$Count <- rep(0, each = nrow(total06OTU))
  for (i in total06OTU$Division) { total06OTU$Count[which(total06OTU$Division == i)] <- sum(dataOTUTaxo  %>% filter(Division == i) %>% select(Total06))}
  total06OTU$Sum <- sum(total06OTU$Count)
  colnames(total09OTU)[2]  <- "value"
  total09OTU$Sum <- rep(0, each = nrow(total09OTU))
  total09OTU$Count <- rep(0, each = nrow(total09OTU))
  for (i in total09OTU$Division) { total09OTU$Count[which(total09OTU$Division == i)] <- sum(dataOTUTaxo  %>% filter(Division == i) %>% select(Total09))}
  total09OTU$Sum <- sum(total09OTU$Count)
  colnames(total11OTU)[2]  <- "value"
  total11OTU$Sum <- rep(0, each = nrow(total11OTU))
  total11OTU$Count <- rep(0, each = nrow(total11OTU))
  for (i in total11OTU$Division) { total11OTU$Count[which(total11OTU$Division == i)] <- sum(dataOTUTaxo  %>% filter(Division == i) %>% select(Total11))}
  total11OTU$Sum <- sum(total11OTU$Count)
  totalDatesTOTU <- rbind(total04OTU,total06OTU,total09OTU,total11OTU)
  #Figure
  bw <- ggplot(totalDatesTOTU, mapping = aes(y= value, x = Dates, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + scale_fill_manual(values = palette)+ 
    geom_label(aes(y = 104,label = paste(Sum," OTUs",sep ="")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
  bw <- bw + theme(legend.position = "none") + labs(x="Dates",y="OTUs (%)")
  print(bw) 
    # Coplot -------------------------------------------------------------------
  svglite("Analyse-Composition-Rarefy/HistOnly/Analyse-OTU-Parasite-Total-V4_095.svg",width = 12.00,height = 6.00)
  b_plot <- plot_grid(bz,bx,by,bw,legendOTU, ncol = 5, nrow = 1, rel_widths = c(3,3,3,5,2),rel_heights = c(3))
  print(b_plot)
  dev.off()
  # HIST Sequence ONLY---------------------------------------------------------------------
    # Cycles only -------------------------------------------------------------------
  onlyJourSequence <- dataSequenceTaxo %>% filter(Cycles == "Jour") %>% select(OTU_Id,TotalJour,Division)
  row.names(onlyJourSequence)<-onlyJourSequence$OTU_Id
  onlyJourSequence <- onlyJourSequence %>% select(-OTU_Id)
  onlyJourSequence <- onlyJourSequence %>% group_by(Division) %>% summarise_all(funs(sum))
  onlyJourSequence$TotalJour <- onlyJourSequence$TotalJour*100/sum(onlyJourSequence$TotalJour)
  onlyJourSequence <- onlyJourSequence %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(TotalJour) - 0.5*TotalJour)
  onlyJourSequence$label <- paste(round(onlyJourSequence$TotalJour,1), "%", sep = "")
  for (i in rownames(onlyJourSequence)) {
    if (onlyJourSequence[i,"label"] == "0%") { onlyJourSequence[i,"label"] <- NA}}
  for (i in rownames(onlyJourSequence)) {
    if (is.na(onlyJourSequence[i,"label"]) == FALSE) { onlyJourSequence[i,"label"] <- paste(onlyJourSequence[i,"Division"]," : ",onlyJourSequence[i,"label"], sep = "")}}
  onlyJourSequence$Cycles<- rep("Jour", each = nrow(onlyJourSequence))
  ##
  onlyNuitSequence <- dataSequenceTaxo %>% filter(Cycles == "Nuit") %>% select(OTU_Id,TotalNuit,Division)
  row.names(onlyNuitSequence)<-onlyNuitSequence$OTU_Id
  onlyNuitSequence <- onlyNuitSequence %>% select(-OTU_Id)
  onlyNuitSequence <- onlyNuitSequence %>% group_by(Division) %>% summarise_all(funs(sum))
  onlyNuitSequence$TotalNuit <- onlyNuitSequence$TotalNuit*100/sum(onlyNuitSequence$TotalNuit)
  onlyNuitSequence <- onlyNuitSequence %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(TotalNuit) - 0.5*TotalNuit)
  onlyNuitSequence$label <- paste(round(onlyNuitSequence$TotalNuit,1), "%", sep = "")
  for (i in rownames(onlyNuitSequence)) {
    if (onlyNuitSequence[i,"label"] == "0%") { onlyNuitSequence[i,"label"] <- NA}}
  for (i in rownames(onlyNuitSequence)) {
    if (is.na(onlyNuitSequence[i,"label"]) == FALSE) { onlyNuitSequence[i,"label"] <- paste(onlyNuitSequence[i,"Division"]," : ",onlyNuitSequence[i,"label"], sep = "")}}
  onlyNuitSequence$Cycles<- rep("Nuit", each = nrow(onlyNuitSequence))
  ##
  colnames(onlyJourSequence)[2]  <- "value"
  onlyJourSequence$Sum <- rep(0, each = nrow(onlyJourSequence))
  onlyJourSequence$Count <- rep(0, each = nrow(onlyJourSequence))
  for (i in onlyJourSequence$Division) { onlyJourSequence$Count[which(onlyJourSequence$Division == i)] <- sum(dataSequenceTaxo %>% filter(Cycles == "Jour") %>% filter(Division == i) %>% select(TotalJour))}
  onlyJourSequence$Sum <- sum(onlyJourSequence$Count)
  colnames(onlyNuitSequence)[2]  <- "value"
  onlyNuitSequence$Sum <- rep(0, each = nrow(onlyNuitSequence))
  onlyNuitSequence$Count <- rep(0, each = nrow(onlyNuitSequence))
  for (i in onlyNuitSequence$Division) { onlyNuitSequence$Count[which(onlyNuitSequence$Division == i)] <- sum(dataSequenceTaxo %>% filter(Cycles == "Nuit") %>% filter(Division == i) %>% select(TotalNuit))}
  onlyNuitSequence$Sum <- sum(onlyNuitSequence$Count)
  onlyCyclesSequence <- rbind(onlyJourSequence,onlyNuitSequence)
  #Figure
  iz <- ggplot(onlyCyclesSequence, mapping = aes(y= value, x = Cycles, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
    geom_label(aes(y = 106,label = paste(Sum,"séquences",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
    scale_fill_manual(values = palette)
  legendSequence <- get_legend(iz)
  iz <- iz + theme(legend.position = "none") + labs(x="Cycles",y="Séquences (%)")
  print(iz)
  
    # FractionO only -------------------------------------------------------------------
  onlyOxiqueSequence <- dataSequenceTaxo %>% filter(FractionO == "Oxique") %>% select(OTU_Id,TotalOxique,Division)
  row.names(onlyOxiqueSequence)<-onlyOxiqueSequence$OTU_Id
  onlyOxiqueSequence <- onlyOxiqueSequence %>% select(-OTU_Id)
  onlyOxiqueSequence <- onlyOxiqueSequence %>% group_by(Division) %>% summarise_all(funs(sum))
  onlyOxiqueSequence$TotalOxique <- onlyOxiqueSequence$TotalOxique*100/sum(onlyOxiqueSequence$TotalOxique)
  onlyOxiqueSequence <- onlyOxiqueSequence %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(TotalOxique) - 0.5*TotalOxique)
  onlyOxiqueSequence$label <- paste(round(onlyOxiqueSequence$TotalOxique,1), "%", sep = "")
  for (i in rownames(onlyOxiqueSequence)) {
    if (onlyOxiqueSequence[i,"label"] == "0%") { onlyOxiqueSequence[i,"label"] <- NA}}
  for (i in rownames(onlyOxiqueSequence)) {
    if (is.na(onlyOxiqueSequence[i,"label"]) == FALSE) { onlyOxiqueSequence[i,"label"] <- paste(onlyOxiqueSequence[i,"Division"]," : ",onlyOxiqueSequence[i,"label"], sep = "")}}
  onlyOxiqueSequence$FractionO<- rep("Oxique", each = nrow(onlyOxiqueSequence))
  ##
  onlyAnoxiqueSequence <- dataSequenceTaxo %>% filter(FractionO == "Anoxique") %>% select(OTU_Id,TotalAnoxique,Division)
  row.names(onlyAnoxiqueSequence)<-onlyAnoxiqueSequence$OTU_Id
  onlyAnoxiqueSequence <- onlyAnoxiqueSequence %>% select(-OTU_Id)
  onlyAnoxiqueSequence <- onlyAnoxiqueSequence %>% group_by(Division) %>% summarise_all(funs(sum))
  onlyAnoxiqueSequence$TotalAnoxique <- onlyAnoxiqueSequence$TotalAnoxique*100/sum(onlyAnoxiqueSequence$TotalAnoxique)
  onlyAnoxiqueSequence <- onlyAnoxiqueSequence %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(TotalAnoxique) - 0.5*TotalAnoxique)
  onlyAnoxiqueSequence$label <- paste(round(onlyAnoxiqueSequence$TotalAnoxique,1), "%", sep = "")
  for (i in rownames(onlyAnoxiqueSequence)) {
    if (onlyAnoxiqueSequence[i,"label"] == "0%") { onlyAnoxiqueSequence[i,"label"] <- NA}}
  for (i in rownames(onlyAnoxiqueSequence)) {
    if (is.na(onlyAnoxiqueSequence[i,"label"]) == FALSE) { onlyAnoxiqueSequence[i,"label"] <- paste(onlyAnoxiqueSequence[i,"Division"]," : ",onlyAnoxiqueSequence[i,"label"], sep = "")}}
  onlyAnoxiqueSequence$FractionO<- rep("Anoxique", each = nrow(onlyAnoxiqueSequence))
  ##
  colnames(onlyOxiqueSequence)[2]  <- "value"
  onlyOxiqueSequence$Sum <- rep(0, each = nrow(onlyOxiqueSequence))
  onlyOxiqueSequence$Count <- rep(0, each = nrow(onlyOxiqueSequence))
  for (i in onlyOxiqueSequence$Division) { onlyOxiqueSequence$Count[which(onlyOxiqueSequence$Division == i)] <- sum(dataSequenceTaxo %>% filter(FractionO == "Oxique") %>% filter(Division == i) %>% select(TotalOxique))}
  onlyOxiqueSequence$Sum <- sum(onlyOxiqueSequence$Count)
  colnames(onlyAnoxiqueSequence)[2]  <- "value"
  onlyAnoxiqueSequence$Sum <- rep(0, each = nrow(onlyAnoxiqueSequence))
  onlyAnoxiqueSequence$Count <- rep(0, each = nrow(onlyAnoxiqueSequence))
  for (i in onlyAnoxiqueSequence$Division) { onlyAnoxiqueSequence$Count[which(onlyAnoxiqueSequence$Division == i)] <- sum(dataSequenceTaxo %>% filter(FractionO == "Anoxique") %>% filter(Division == i) %>% select(TotalAnoxique))}
  onlyAnoxiqueSequence$Sum <- sum(onlyAnoxiqueSequence$Count)
  onlyFractionOSequence <- rbind(onlyOxiqueSequence,onlyAnoxiqueSequence)
  #Figure
  ix <- ggplot(onlyFractionOSequence, mapping = aes(y= value, x = FractionO, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
    geom_label(aes(y = 106,label = paste(Sum,"séquences",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
    scale_fill_manual(values = palette)
  ix <- ix + theme(legend.position = "none") + labs(x="Zones",y="Séquences (%)")
  print(ix)    
    # FractionT only -------------------------------------------------------------------
  onlyPetiteSequence <- dataSequenceTaxo %>% filter(FractionT == "Petite") %>% select(OTU_Id,TotalPetite,Division)
  row.names(onlyPetiteSequence)<-onlyPetiteSequence$OTU_Id
  onlyPetiteSequence <- onlyPetiteSequence %>% select(-OTU_Id)
  onlyPetiteSequence <- onlyPetiteSequence %>% group_by(Division) %>% summarise_all(funs(sum))
  onlyPetiteSequence$TotalPetite <- onlyPetiteSequence$TotalPetite*100/sum(onlyPetiteSequence$TotalPetite)
  onlyPetiteSequence <- onlyPetiteSequence %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(TotalPetite) - 0.5*TotalPetite)
  onlyPetiteSequence$label <- paste(round(onlyPetiteSequence$TotalPetite,1), "%", sep = "")
  for (i in rownames(onlyPetiteSequence)) {
    if (onlyPetiteSequence[i,"label"] == "0%") { onlyPetiteSequence[i,"label"] <- NA}}
  for (i in rownames(onlyPetiteSequence)) {
    if (is.na(onlyPetiteSequence[i,"label"]) == FALSE) { onlyPetiteSequence[i,"label"] <- paste(onlyPetiteSequence[i,"Division"]," : ",onlyPetiteSequence[i,"label"], sep = "")}}
  onlyPetiteSequence$FractionT<- rep("Petite", each = nrow(onlyPetiteSequence))
  ##
  onlyGrandeSequence <- dataSequenceTaxo %>% filter(FractionT == "Grande") %>% select(OTU_Id,TotalGrande,Division)
  row.names(onlyGrandeSequence)<-onlyGrandeSequence$OTU_Id
  onlyGrandeSequence <- onlyGrandeSequence %>% select(-OTU_Id)
  onlyGrandeSequence <- onlyGrandeSequence %>% group_by(Division) %>% summarise_all(funs(sum))
  onlyGrandeSequence$TotalGrande <- onlyGrandeSequence$TotalGrande*100/sum(onlyGrandeSequence$TotalGrande)
  onlyGrandeSequence <- onlyGrandeSequence %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(TotalGrande) - 0.5*TotalGrande)
  onlyGrandeSequence$label <- paste(round(onlyGrandeSequence$TotalGrande,1), "%", sep = "")
  for (i in rownames(onlyGrandeSequence)) {
    if (onlyGrandeSequence[i,"label"] == "0%") { onlyGrandeSequence[i,"label"] <- NA}}
  for (i in rownames(onlyGrandeSequence)) {
    if (is.na(onlyGrandeSequence[i,"label"]) == FALSE) { onlyGrandeSequence[i,"label"] <- paste(onlyGrandeSequence[i,"Division"]," : ",onlyGrandeSequence[i,"label"], sep = "")}}
  onlyGrandeSequence$FractionT<- rep("Grande", each = nrow(onlyGrandeSequence))
  ##
  colnames(onlyPetiteSequence)[2]  <- "value"
  onlyPetiteSequence$Sum <- rep(0, each = nrow(onlyPetiteSequence))
  onlyPetiteSequence$Count <- rep(0, each = nrow(onlyPetiteSequence))
  for (i in onlyPetiteSequence$Division) { onlyPetiteSequence$Count[which(onlyPetiteSequence$Division == i)] <- sum(dataSequenceTaxo %>% filter(FractionT == "Petite") %>% filter(Division == i) %>% select(TotalPetite))}
  onlyPetiteSequence$Sum <- sum(onlyPetiteSequence$Count)
  colnames(onlyGrandeSequence)[2]  <- "value"
  onlyGrandeSequence$Sum <- rep(0, each = nrow(onlyGrandeSequence))
  onlyGrandeSequence$Count <- rep(0, each = nrow(onlyGrandeSequence))
  for (i in onlyGrandeSequence$Division) { onlyGrandeSequence$Count[which(onlyGrandeSequence$Division == i)] <- sum(dataSequenceTaxo %>% filter(FractionT == "Grande") %>% filter(Division == i) %>% select(TotalGrande))}
  onlyGrandeSequence$Sum <- sum(onlyGrandeSequence$Count)
  onlyFractionTSequence <- rbind(onlyPetiteSequence,onlyGrandeSequence)
  #Figure
  iy <- ggplot(onlyFractionTSequence, mapping = aes(y= value, x = FractionT, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
    geom_label(aes(y = 106,label = paste(Sum,"séquences",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
    scale_fill_manual(values = palette)
  iy <- iy + theme(legend.position = "none") + labs(x="Fractions",y="Séquences (%)")
  print(iy) 
    # Coplot -------------------------------------------------------------------
  svglite("Analyse-Composition-Rarefy/HistOnly/Analyse-Sequence-Parasite-only-V4_095.svg",width = 8.00,height = 6.00)
  b_plot <- plot_grid(iz,ix,iy,legendSequence,  ncol = 4, nrow = 1, rel_widths = c(3,3,3,2),rel_heights = c(3))
  print(b_plot)
  dev.off()  

    # Coplot X AFC -------------------------------------------------------------------
  f <- ggplot() + theme_void()
  print(f)
  
  svglite("Analyse-Composition-Rarefy/Biplot/Analyse-Sequence-Parasite-onlyX-V4_095.svg",width = 10.00,height = 13.00)
  b_plot <- plot_grid(a,iz,f,b,ix,legendSequence,c,iy,f, ncol = 3, nrow = 3, rel_widths = c(3,1.5,1),rel_heights = c(3,3,3))
  print(b_plot)
  dev.off()  
  
  svglite("Analyse-Composition-Rarefy/Biplot/Analyse-OTUXSequence-Parasite-onlyX-V4_095.svg",width = 10.00,height = 13.00)
  b_plot <- plot_grid(as,iz,f,bs,ix,legendSequence,cs,iy,f, ncol = 3, nrow = 3, rel_widths = c(3,1.5,1),rel_heights = c(3,3,3))
  print(b_plot)
  dev.off()  
  
  
  # HIST Sequence TOTAL---------------------------------------------------------------------
    # Cycles -------------------------------------------------------------------
  totalJourSequence <- dataSequenceTaxo %>% select(OTU_Id,TotalJour,Division)
  row.names(totalJourSequence)<-totalJourSequence$OTU_Id
  totalJourSequence <- totalJourSequence %>% select(-OTU_Id)
  totalJourSequence <- totalJourSequence %>% group_by(Division) %>% summarise_all(funs(sum))
  totalJourSequence$TotalJour <- totalJourSequence$TotalJour*100/sum(totalJourSequence$TotalJour)
  totalJourSequence <- totalJourSequence %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(TotalJour) - 0.5*TotalJour)
  totalJourSequence$label <- paste(round(totalJourSequence$TotalJour,1), "%", sep = "")
  for (i in rownames(totalJourSequence)) {
    if (totalJourSequence[i,"label"] == "0%") { totalJourSequence[i,"label"] <- NA}}
  for (i in rownames(totalJourSequence)) {
    if (is.na(totalJourSequence[i,"label"]) == FALSE) { totalJourSequence[i,"label"] <- paste(totalJourSequence[i,"Division"]," : ",totalJourSequence[i,"label"], sep = "")}}
  totalJourSequence$Cycles<- rep("Jour", each = nrow(totalJourSequence))
  ##
  totalNuitSequence <- dataSequenceTaxo %>% select(OTU_Id,TotalNuit,Division)
  row.names(totalNuitSequence)<-totalNuitSequence$OTU_Id
  totalNuitSequence <- totalNuitSequence %>% select(-OTU_Id)
  totalNuitSequence <- totalNuitSequence %>% group_by(Division) %>% summarise_all(funs(sum))
  totalNuitSequence$TotalNuit <- totalNuitSequence$TotalNuit*100/sum(totalNuitSequence$TotalNuit)
  totalNuitSequence <- totalNuitSequence %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(TotalNuit) - 0.5*TotalNuit)
  totalNuitSequence$label <- paste(round(totalNuitSequence$TotalNuit,1), "%", sep = "")
  for (i in rownames(totalNuitSequence)) {
    if (totalNuitSequence[i,"label"] == "0%") { totalNuitSequence[i,"label"] <- NA}}
  for (i in rownames(totalNuitSequence)) {
    if (is.na(totalNuitSequence[i,"label"]) == FALSE) { totalNuitSequence[i,"label"] <- paste(totalNuitSequence[i,"Division"]," : ",totalNuitSequence[i,"label"], sep = "")}}
  totalNuitSequence$Cycles<- rep("Nuit", each = nrow(totalNuitSequence))
  ##
  colnames(totalJourSequence)[2]  <- "value"
  totalJourSequence$Sum <- rep(0, each = nrow(totalJourSequence))
  totalJourSequence$Count <- rep(0, each = nrow(totalJourSequence))
  for (i in totalJourSequence$Division) { totalJourSequence$Count[which(totalJourSequence$Division == i)] <- sum(dataSequenceTaxo  %>% filter(Division == i) %>% select(TotalJour))}
  totalJourSequence$Sum <- sum(totalJourSequence$Count)
  colnames(totalNuitSequence)[2]  <- "value"
  totalNuitSequence$Sum <- rep(0, each = nrow(totalNuitSequence))
  totalNuitSequence$Count <- rep(0, each = nrow(totalNuitSequence))
  for (i in totalNuitSequence$Division) { totalNuitSequence$Count[which(totalNuitSequence$Division == i)] <- sum(dataSequenceTaxo  %>% filter(Division == i) %>% select(TotalNuit))}
  totalNuitSequence$Sum <- sum(totalNuitSequence$Count)
  totalCyclesSequence <- rbind(totalJourSequence,totalNuitSequence)
  totalCyclesSequence$percent <- paste("(",round(totalCyclesSequence$Sum*100/(colSums(Sum %>% select(Séquences))/2),1)," %)",sep ="")
  
  #Figure
  jz <- ggplot(totalCyclesSequence, mapping = aes(y= value, x = Cycles, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
    scale_fill_manual(values = palette) + 
    geom_label(aes(y = 106,label = paste(Sum,"séquences",percent,sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
  legendSequence <- get_legend(jz)
  jz <- jz + theme(legend.position = "none") + labs(x="Cycles",y="Séquences (%)")
  print(jz)
  
    # FractionO -------------------------------------------------------------------
  totalOxiqueSequence <- dataSequenceTaxo %>% select(OTU_Id,TotalOxique,Division)
  row.names(totalOxiqueSequence)<-totalOxiqueSequence$OTU_Id
  totalOxiqueSequence <- totalOxiqueSequence %>% select(-OTU_Id)
  totalOxiqueSequence <- totalOxiqueSequence %>% group_by(Division) %>% summarise_all(funs(sum))
  totalOxiqueSequence$TotalOxique <- totalOxiqueSequence$TotalOxique*100/sum(totalOxiqueSequence$TotalOxique)
  totalOxiqueSequence <- totalOxiqueSequence %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(TotalOxique) - 0.5*TotalOxique)
  totalOxiqueSequence$label <- paste(round(totalOxiqueSequence$TotalOxique,1), "%", sep = "")
  for (i in rownames(totalOxiqueSequence)) {
    if (totalOxiqueSequence[i,"label"] == "0%") { totalOxiqueSequence[i,"label"] <- NA}}
  for (i in rownames(totalOxiqueSequence)) {
    if (is.na(totalOxiqueSequence[i,"label"]) == FALSE) { totalOxiqueSequence[i,"label"] <- paste(totalOxiqueSequence[i,"Division"]," : ",totalOxiqueSequence[i,"label"], sep = "")}}
  totalOxiqueSequence$FractionO<- rep("Oxique", each = nrow(totalOxiqueSequence))
  ##
  totalAnoxiqueSequence <- dataSequenceTaxo %>% select(OTU_Id,TotalAnoxique,Division)
  row.names(totalAnoxiqueSequence)<-totalAnoxiqueSequence$OTU_Id
  totalAnoxiqueSequence <- totalAnoxiqueSequence %>% select(-OTU_Id)
  totalAnoxiqueSequence <- totalAnoxiqueSequence %>% group_by(Division) %>% summarise_all(funs(sum))
  totalAnoxiqueSequence$TotalAnoxique <- totalAnoxiqueSequence$TotalAnoxique*100/sum(totalAnoxiqueSequence$TotalAnoxique)
  totalAnoxiqueSequence <- totalAnoxiqueSequence %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(TotalAnoxique) - 0.5*TotalAnoxique)
  totalAnoxiqueSequence$label <- paste(round(totalAnoxiqueSequence$TotalAnoxique,1), "%", sep = "")
  for (i in rownames(totalAnoxiqueSequence)) {
    if (totalAnoxiqueSequence[i,"label"] == "0%") { totalAnoxiqueSequence[i,"label"] <- NA}}
  for (i in rownames(totalAnoxiqueSequence)) {
    if (is.na(totalAnoxiqueSequence[i,"label"]) == FALSE) { totalAnoxiqueSequence[i,"label"] <- paste(totalAnoxiqueSequence[i,"Division"]," : ",totalAnoxiqueSequence[i,"label"], sep = "")}}
  totalAnoxiqueSequence$FractionO<- rep("Anoxique", each = nrow(totalAnoxiqueSequence))
  ##
  colnames(totalOxiqueSequence)[2]  <- "value"
  totalOxiqueSequence$Sum <- rep(0, each = nrow(totalOxiqueSequence))
  totalOxiqueSequence$Count <- rep(0, each = nrow(totalOxiqueSequence))
  for (i in totalOxiqueSequence$Division) { totalOxiqueSequence$Count[which(totalOxiqueSequence$Division == i)] <- sum(dataSequenceTaxo  %>% filter(Division == i) %>% select(TotalOxique))}
  totalOxiqueSequence$Sum <- sum(totalOxiqueSequence$Count)
  colnames(totalAnoxiqueSequence)[2]  <- "value"
  totalAnoxiqueSequence$Sum <- rep(0, each = nrow(totalAnoxiqueSequence))
  totalAnoxiqueSequence$Count <- rep(0, each = nrow(totalAnoxiqueSequence))
  for (i in totalAnoxiqueSequence$Division) { totalAnoxiqueSequence$Count[which(totalAnoxiqueSequence$Division == i)] <- sum(dataSequenceTaxo  %>% filter(Division == i) %>% select(TotalAnoxique))}
  totalAnoxiqueSequence$Sum <- sum(totalAnoxiqueSequence$Count)
  totalFractionOSequence <- rbind(totalOxiqueSequence,totalAnoxiqueSequence)
  totalFractionOSequence$percent <- paste("(",round(totalFractionOSequence$Sum*100/(colSums(Sum %>% select(Séquences))/2),1)," %)",sep ="")
  
  #Figure
  jx <- ggplot(totalFractionOSequence, mapping = aes(y= value, x = FractionO, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
    scale_fill_manual(values = palette) + 
    geom_label(aes(y = 106,label = paste(Sum,"séquences",percent,sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
  jx <- jx + theme(legend.position = "none") + labs(x="Zones",y="Séquences (%)")
  print(jx)    
    # FractionT -------------------------------------------------------------------
  totalPetiteSequence <- dataSequenceTaxo %>% select(OTU_Id,TotalPetite,Division)
  row.names(totalPetiteSequence)<-totalPetiteSequence$OTU_Id
  totalPetiteSequence <- totalPetiteSequence %>% select(-OTU_Id)
  totalPetiteSequence <- totalPetiteSequence %>% group_by(Division) %>% summarise_all(funs(sum))
  totalPetiteSequence$TotalPetite <- totalPetiteSequence$TotalPetite*100/sum(totalPetiteSequence$TotalPetite)
  totalPetiteSequence <- totalPetiteSequence %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(TotalPetite) - 0.5*TotalPetite)
  totalPetiteSequence$label <- paste(round(totalPetiteSequence$TotalPetite,1), "%", sep = "")
  for (i in rownames(totalPetiteSequence)) {
    if (totalPetiteSequence[i,"label"] == "0%") { totalPetiteSequence[i,"label"] <- NA}}
  for (i in rownames(totalPetiteSequence)) {
    if (is.na(totalPetiteSequence[i,"label"]) == FALSE) { totalPetiteSequence[i,"label"] <- paste(totalPetiteSequence[i,"Division"]," : ",totalPetiteSequence[i,"label"], sep = "")}}
  totalPetiteSequence$FractionT<- rep("Petite", each = nrow(totalPetiteSequence))
  ##
  totalGrandeSequence <- dataSequenceTaxo %>% select(OTU_Id,TotalGrande,Division)
  row.names(totalGrandeSequence)<-totalGrandeSequence$OTU_Id
  totalGrandeSequence <- totalGrandeSequence %>% select(-OTU_Id)
  totalGrandeSequence <- totalGrandeSequence %>% group_by(Division) %>% summarise_all(funs(sum))
  totalGrandeSequence$TotalGrande <- totalGrandeSequence$TotalGrande*100/sum(totalGrandeSequence$TotalGrande)
  totalGrandeSequence <- totalGrandeSequence %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(TotalGrande) - 0.5*TotalGrande)
  totalGrandeSequence$label <- paste(round(totalGrandeSequence$TotalGrande,1), "%", sep = "")
  for (i in rownames(totalGrandeSequence)) {
    if (totalGrandeSequence[i,"label"] == "0%") { totalGrandeSequence[i,"label"] <- NA}}
  for (i in rownames(totalGrandeSequence)) {
    if (is.na(totalGrandeSequence[i,"label"]) == FALSE) { totalGrandeSequence[i,"label"] <- paste(totalGrandeSequence[i,"Division"]," : ",totalGrandeSequence[i,"label"], sep = "")}}
  totalGrandeSequence$FractionT<- rep("Grande", each = nrow(totalGrandeSequence))
  ##
  colnames(totalPetiteSequence)[2]  <- "value"
  totalPetiteSequence$Sum <- rep(0, each = nrow(totalPetiteSequence))
  totalPetiteSequence$Count <- rep(0, each = nrow(totalPetiteSequence))
  for (i in totalPetiteSequence$Division) { totalPetiteSequence$Count[which(totalPetiteSequence$Division == i)] <- sum(dataSequenceTaxo  %>% filter(Division == i) %>% select(TotalPetite))}
  totalPetiteSequence$Sum <- sum(totalPetiteSequence$Count)
  colnames(totalGrandeSequence)[2]  <- "value"
  totalGrandeSequence$Sum <- rep(0, each = nrow(totalGrandeSequence))
  totalGrandeSequence$Count <- rep(0, each = nrow(totalGrandeSequence))
  for (i in totalGrandeSequence$Division) { totalGrandeSequence$Count[which(totalGrandeSequence$Division == i)] <- sum(dataSequenceTaxo  %>% filter(Division == i) %>% select(TotalGrande))}
  totalGrandeSequence$Sum <- sum(totalGrandeSequence$Count)
  totalFractionTSequence <- rbind(totalPetiteSequence,totalGrandeSequence)
  totalFractionTSequence$percent <- paste("(",round(totalFractionTSequence$Sum*100/(colSums(Sum %>% select(Séquences))/2),1)," %)",sep ="")
  
  #Figure
  jy <- ggplot(totalFractionTSequence, mapping = aes(y= value, x = FractionT, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
    scale_fill_manual(values = palette) + 
    geom_label(aes(y = 106,label = paste(Sum,"séquences",percent,sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
  jy <- jy + theme(legend.position = "none") + labs(x="Fractions",y="Séquences (%)")
  print(jy) 
    # Dates -------------------------------------------------------------------
  total04Sequence <- dataSequenceTaxo %>% select(OTU_Id,Total04,Division)
  row.names(total04Sequence)<-total04Sequence$OTU_Id
  total04Sequence <- total04Sequence %>% select(-OTU_Id)
  total04Sequence <- total04Sequence %>% group_by(Division) %>% summarise_all(funs(sum))
  total04Sequence$Total04 <- total04Sequence$Total04*100/sum(total04Sequence$Total04)
  total04Sequence <- total04Sequence %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(Total04) - 0.5*Total04)
  total04Sequence$label <- paste(round(total04Sequence$Total04,1), "%", sep = "")
  for (i in rownames(total04Sequence)) {
    if (total04Sequence[i,"label"] == "0%") { total04Sequence[i,"label"] <- NA}}
  for (i in rownames(total04Sequence)) {
    if (is.na(total04Sequence[i,"label"]) == FALSE) { total04Sequence[i,"label"] <- paste(total04Sequence[i,"Division"]," : ",total04Sequence[i,"label"], sep = "")}}
  total04Sequence$Dates<- rep("04", each = nrow(total04Sequence))
  ##
  total06Sequence <- dataSequenceTaxo %>% select(OTU_Id,Total06,Division)
  row.names(total06Sequence)<-total06Sequence$OTU_Id
  total06Sequence <- total06Sequence %>% select(-OTU_Id)
  total06Sequence <- total06Sequence %>% group_by(Division) %>% summarise_all(funs(sum))
  total06Sequence$Total06 <- total06Sequence$Total06*100/sum(total06Sequence$Total06)
  total06Sequence <- total06Sequence %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(Total06) - 0.5*Total06)
  total06Sequence$label <- paste(round(total06Sequence$Total06,1), "%", sep = "")
  for (i in rownames(total06Sequence)) {
    if (total06Sequence[i,"label"] == "0%") { total06Sequence[i,"label"] <- NA}}
  for (i in rownames(total06Sequence)) {
    if (is.na(total06Sequence[i,"label"]) == FALSE) { total06Sequence[i,"label"] <- paste(total06Sequence[i,"Division"]," : ",total06Sequence[i,"label"], sep = "")}}
  total06Sequence$Dates<- rep("06", each = nrow(total06Sequence))
  ##
  total09Sequence <- dataSequenceTaxo %>% select(OTU_Id,Total09,Division)
  row.names(total09Sequence)<-total09Sequence$OTU_Id
  total09Sequence <- total09Sequence %>% select(-OTU_Id)
  total09Sequence <- total09Sequence %>% group_by(Division) %>% summarise_all(funs(sum))
  total09Sequence$Total09 <- total09Sequence$Total09*100/sum(total09Sequence$Total09)
  total09Sequence <- total09Sequence %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(Total09) - 0.5*Total09)
  total09Sequence$label <- paste(round(total09Sequence$Total09,1), "%", sep = "")
  for (i in rownames(total09Sequence)) {
    if (total09Sequence[i,"label"] == "0%") { total09Sequence[i,"label"] <- NA}}
  for (i in rownames(total09Sequence)) {
    if (is.na(total09Sequence[i,"label"]) == FALSE) { total09Sequence[i,"label"] <- paste(total09Sequence[i,"Division"]," : ",total09Sequence[i,"label"], sep = "")}}
  total09Sequence$Dates<- rep("09", each = nrow(total09Sequence))
  ##
  total11Sequence <- dataSequenceTaxo %>% select(OTU_Id,Total11,Division)
  row.names(total11Sequence)<-total11Sequence$OTU_Id
  total11Sequence <- total11Sequence %>% select(-OTU_Id)
  total11Sequence <- total11Sequence %>% group_by(Division) %>% summarise_all(funs(sum))
  total11Sequence$Total11 <- total11Sequence$Total11*100/sum(total11Sequence$Total11)
  total11Sequence <- total11Sequence %>%
    arrange(desc(Division)) %>%
    mutate(lab.ypos = cumsum(Total11) - 0.5*Total11)
  total11Sequence$label <- paste(round(total11Sequence$Total11,1), "%", sep = "")
  for (i in rownames(total11Sequence)) {
    if (total11Sequence[i,"label"] == "0%") { total11Sequence[i,"label"] <- NA}}
  for (i in rownames(total11Sequence)) {
    if (is.na(total11Sequence[i,"label"]) == FALSE) { total11Sequence[i,"label"] <- paste(total11Sequence[i,"Division"]," : ",total11Sequence[i,"label"], sep = "")}}
  total11Sequence$Dates<- rep("11", each = nrow(total11Sequence))
  ##
  colnames(total04Sequence)[2]  <- "value"
  total04Sequence$Sum <- rep(0, each = nrow(total04Sequence))
  total04Sequence$Count <- rep(0, each = nrow(total04Sequence))
  for (i in total04Sequence$Division) { total04Sequence$Count[which(total04Sequence$Division == i)] <- sum(dataSequenceTaxo  %>% filter(Division == i) %>% select(Total04))}
  total04Sequence$Sum <- sum(total04Sequence$Count)
  colnames(total06Sequence)[2]  <- "value"
  total06Sequence$Sum <- rep(0, each = nrow(total06Sequence))
  total06Sequence$Count <- rep(0, each = nrow(total06Sequence))
  for (i in total06Sequence$Division) { total06Sequence$Count[which(total06Sequence$Division == i)] <- sum(dataSequenceTaxo  %>% filter(Division == i) %>% select(Total06))}
  total06Sequence$Sum <- sum(total06Sequence$Count)
  colnames(total09Sequence)[2]  <- "value"
  total09Sequence$Sum <- rep(0, each = nrow(total09Sequence))
  total09Sequence$Count <- rep(0, each = nrow(total09Sequence))
  for (i in total09Sequence$Division) { total09Sequence$Count[which(total09Sequence$Division == i)] <- sum(dataSequenceTaxo  %>% filter(Division == i) %>% select(Total09))}
  total09Sequence$Sum <- sum(total09Sequence$Count)
  colnames(total11Sequence)[2]  <- "value"
  total11Sequence$Sum <- rep(0, each = nrow(total11Sequence))
  total11Sequence$Count <- rep(0, each = nrow(total11Sequence))
  for (i in total11Sequence$Division) { total11Sequence$Count[which(total11Sequence$Division == i)] <- sum(dataSequenceTaxo  %>% filter(Division == i) %>% select(Total11))}
  total11Sequence$Sum <- sum(total11Sequence$Count)
  totalDatesTSequence <- rbind(total04Sequence,total06Sequence,total09Sequence,total11Sequence)
  totalDatesTSequence$percent <- paste("(",round(totalDatesTSequence$Sum*100/(colSums(Sum %>% select(Séquences))/4),1)," %)",sep ="")
  
  #Figure
  jw <- ggplot(totalDatesTSequence, mapping = aes(y= value, x = Dates, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
    scale_fill_manual(values = palette) + 
    geom_label(aes(y = 106,label = paste(Sum,"séquences",percent,sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
  jw <- jw + theme(legend.position = "none") + labs(x="Dates",y="Séquences (%)")
  print(jw) 
    # Coplot -------------------------------------------------------------------
  svglite("Analyse-Composition-Rarefy/HistOnly/Analyse-Sequence-Parasite-Total-V4_095.svg",width = 12.00,height = 7.00)
  b_plot <- plot_grid(jz,jx,jy,jw, legendSequence, ncol = 5, nrow = 1, rel_widths = c(3,3,3,5,2),rel_heights = c(3))
  print(b_plot)
  dev.off()
# Table OTUs majoritaires -------------------------------------------------
  # Table ONLY---------------------------------------------------------------------
  totaltaxo <- tableV4095 %>% select(OTU_Id,NN.taxonomy,LCA.taxonomy,Best_hit_identity,Identity....)
    # Cycles only -------------------------------------------------------------------
  onlyJourTable <- dataSequenceTaxo %>% filter(Cycles == "Jour") %>% select(OTU_Id,TotalJour,Division,Cycles)
  x <- sum(onlyJourTable$TotalJour)
  onlyJourTable <- onlyJourTable %>% filter(TotalJour > 0.01*x)
  colnames(onlyJourTable)[2]  <- "value"
  ##
  onlyNuitTable <- dataSequenceTaxo %>% filter(Cycles == "Nuit") %>% select(OTU_Id,TotalNuit,Division,Cycles)
  y <- sum(onlyNuitTable$TotalNuit)
  onlyNuitTable <- onlyNuitTable %>% filter(TotalNuit > 0.01*y)
  colnames(onlyNuitTable)[2]  <- "value"
  onlyCyclesTable <- rbind(onlyJourTable,onlyNuitTable)
  onlyCyclesTable <- merge(onlyCyclesTable,totaltaxo, by = "OTU_Id")
  ##
  write.table(onlyCyclesTable, file = "Analyse-Composition-Rarefy/TableOnly/TableOnlyCycles.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(onlyCyclesOTU, file = "Analyse-Composition-Rarefy/TableOnly/OnlyCyclesOTU.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(onlyCyclesSequence, file = "Analyse-Composition-Rarefy/TableOnly/OnlyCyclesSequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
    # FractionO only -------------------------------------------------------------------
  onlyOxiqueTable <- dataSequenceTaxo %>% filter(FractionO == "Oxique") %>% select(OTU_Id,TotalOxique,Division,FractionO)
  x <- sum(onlyOxiqueTable$TotalOxique)
  onlyOxiqueTable <- onlyOxiqueTable %>% filter(TotalOxique > 0.01*x)
  colnames(onlyOxiqueTable)[2]  <- "value"
  ##
  onlyAnoxiqueTable <- dataSequenceTaxo %>% filter(FractionO == "Anoxique") %>% select(OTU_Id,TotalAnoxique,Division,FractionO)
  y <- sum(onlyAnoxiqueTable$TotalAnoxique)
  onlyAnoxiqueTable <- onlyAnoxiqueTable %>% filter(TotalAnoxique > 0.01*y)
  colnames(onlyAnoxiqueTable)[2]  <- "value"
  onlyFractionOTable <- rbind(onlyOxiqueTable,onlyAnoxiqueTable)
  onlyFractionOTable <- merge(onlyFractionOTable,totaltaxo, by = "OTU_Id")
  ##
  write.table(onlyFractionOTable, file = "Analyse-Composition-Rarefy/TableOnly/TableOnlyFractionO.txt", sep = "\t", col.names = TRUE, row.names = FALSE,quote = FALSE)
  write.table(onlyFractionOOTU, file = "Analyse-Composition-Rarefy/TableOnly/OnlyFractionOOTU.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(onlyFractionOSequence, file = "Analyse-Composition-Rarefy/TableOnly/OnlyFractionOSequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
    # FractionT only -------------------------------------------------------------------
  onlyPetiteTable <- dataSequenceTaxo %>% filter(FractionT == "Petite") %>% select(OTU_Id,TotalPetite,Division,FractionT)
  x <- sum(onlyPetiteTable$TotalPetite)
  onlyPetiteTable <- onlyPetiteTable %>% filter(TotalPetite > 0.01*x)
  colnames(onlyPetiteTable)[2]  <- "value"
  ##
  onlyGrandeTable <- dataSequenceTaxo %>% filter(FractionT == "Grande") %>% select(OTU_Id,TotalGrande,Division,FractionT)
  y <- sum(onlyGrandeTable$TotalGrande)
  onlyGrandeTable <- onlyGrandeTable %>% filter(TotalGrande > 0.01*y)
  colnames(onlyGrandeTable)[2]  <- "value"
  onlyFractionTTable <- rbind(onlyPetiteTable,onlyGrandeTable)
  onlyFractionTTable <- merge(onlyFractionTTable,totaltaxo, by = "OTU_Id")
  ##
  write.table(onlyFractionTTable, file = "Analyse-Composition-Rarefy/TableOnly/TableOnlyFractionT.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(onlyFractionTOTU, file = "Analyse-Composition-Rarefy/TableOnly/OnlyFractionTOTU.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(onlyFractionTSequence, file = "Analyse-Composition-Rarefy/TableOnly/OnlyFractionTSequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  # Table Total---------------------------------------------------------------------
  totaltaxo <- tableV4095 %>% select(OTU_Id,NN.taxonomy,LCA.taxonomy,Best_hit_identity,Identity....)
    # Cycles Total -------------------------------------------------------------------
  totalJourTable <- dataSequenceTaxo %>% select(OTU_Id,TotalJour,Division,Cycles)
  totalJourTable$Cycles <- rep("Jour",each = nrow(totalJourTable))
  x <- sum(totalJourTable$TotalJour)
  totalJourTable <- totalJourTable %>% filter(TotalJour > 0.05*x)
  colnames(totalJourTable)[2]  <- "value"
  ##
  totalNuitTable <- dataSequenceTaxo %>% select(OTU_Id,TotalNuit,Division,Cycles)
  totalNuitTable$Cycles <- rep("Nuit",each = nrow(totalNuitTable))
  y <- sum(totalNuitTable$TotalNuit)
  totalNuitTable <- totalNuitTable %>% filter(TotalNuit > 0.05*y)
  colnames(totalNuitTable)[2]  <- "value"
  totalCyclesTable <- rbind(totalJourTable,totalNuitTable)
  totalCyclesTable <- merge(totalCyclesTable,totaltaxo, by = "OTU_Id")
  ##
  write.table(totalCyclesTable, file = "Analyse-Composition-Rarefy/TableOnly/TableTotalCycles.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(totalCyclesOTU, file = "Analyse-Composition-Rarefy/TableOnly/TotalCyclesOTU.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(totalCyclesSequence, file = "Analyse-Composition-Rarefy/TableOnly/TotalCyclesSequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
    # FractionO Total -------------------------------------------------------------------
  totalOxiqueTable <- dataSequenceTaxo %>% select(OTU_Id,TotalOxique,Division,FractionO)
  totalOxiqueTable$FractionO <- rep("Oxique",each = nrow(totalOxiqueTable))
  x <- sum(totalOxiqueTable$TotalOxique)
  totalOxiqueTable <- totalOxiqueTable %>% filter(TotalOxique > 0.05*x)
  colnames(totalOxiqueTable)[2]  <- "value"
  ##
  totalAnoxiqueTable <- dataSequenceTaxo %>% select(OTU_Id,TotalAnoxique,Division,FractionO)
  totalAnoxiqueTable$FractionO <- rep("Anoxique",each = nrow(totalAnoxiqueTable))
  y <- sum(totalAnoxiqueTable$TotalAnoxique)
  totalAnoxiqueTable <- totalAnoxiqueTable %>% filter(TotalAnoxique > 0.05*y)
  colnames(totalAnoxiqueTable)[2]  <- "value"
  totalFractionOTable <- rbind(totalOxiqueTable,totalAnoxiqueTable)
  totalFractionOTable <- merge(totalFractionOTable,totaltaxo, by = "OTU_Id")
  ##
  write.table(totalFractionOTable, file = "Analyse-Composition-Rarefy/TableOnly/TableTotalFractionO.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(totalFractionOOTU, file = "Analyse-Composition-Rarefy/TableOnly/TotalFractionOOTU.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(totalFractionOSequence, file = "Analyse-Composition-Rarefy/TableOnly/TotalFractionOSequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
    # FractionT Total -------------------------------------------------------------------
  totalPetiteTable <- dataSequenceTaxo %>% select(OTU_Id,TotalPetite,Division,FractionT)
  totalPetiteTable$FractionT <- rep("Petite",each = nrow(totalPetiteTable))
  x <- sum(totalPetiteTable$TotalPetite)
  totalPetiteTable <- totalPetiteTable %>% filter(TotalPetite > 0.05*x)
  colnames(totalPetiteTable)[2]  <- "value"
  ##
  totalGrandeTable <- dataSequenceTaxo %>% select(OTU_Id,TotalGrande,Division,FractionT)
  totalGrandeTable$FractionT <- rep("Grande",each = nrow(totalGrandeTable))
  y <- sum(totalGrandeTable$TotalGrande)
  totalGrandeTable <- totalGrandeTable %>% filter(TotalGrande > 0.05*y)
  colnames(totalGrandeTable)[2]  <- "value"
  totalFractionTTable <- rbind(totalPetiteTable,totalGrandeTable)
  totalFractionTTable <- merge(totalFractionTTable,totaltaxo, by = "OTU_Id")
  ##
  write.table(totalFractionTTable, file = "Analyse-Composition-Rarefy/TableOnly/TableTotalFractionT.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(totalFractionTOTU, file = "Analyse-Composition-Rarefy/TableOnly/TotalFractionTOTU.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(totalFractionTSequence, file = "Analyse-Composition-Rarefy/TableOnly/TotalFractionTSequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
    # Dates Total -------------------------------------------------------------------
  total04Table <- dataSequenceTaxo %>% select(OTU_Id,Total04,Division)
  total04Table$Dates <- rep(04,each = nrow(total04Table))
  x <- sum(total04Table$Total04)
  total04Table <- total04Table %>% filter(Total04 > 0.05*x)
  colnames(total04Table)[2]  <- "value"
  ##
  total06Table <- dataSequenceTaxo %>% select(OTU_Id,Total06,Division)
  total06Table$Dates <- rep(06,each = nrow(total06Table))
  y <- sum(total06Table$Total06)
  total06Table <- total06Table %>% filter(Total06 > 0.05*y)
  colnames(total06Table)[2]  <- "value"
  ##
  total09Table <- dataSequenceTaxo %>% select(OTU_Id,Total09,Division)
  total09Table$Dates <- rep(09,each = nrow(total09Table))
  z <- sum(total09Table$Total09)
  total09Table <- total09Table %>% filter(Total09 > 0.05*z)
  colnames(total09Table)[2]  <- "value"
  ##
  total11Table <- dataSequenceTaxo %>% select(OTU_Id,Total11,Division)
  total11Table$Dates <- rep(11,each = nrow(total11Table))
  w <- sum(total11Table$Total11)
  total11Table <- total11Table %>% filter(Total11 > 0.05*w)
  colnames(total11Table)[2]  <- "value"
  ##
  totalDatesTable <- rbind(total04Table,total06Table,total09Table,total11Table)
  totalDatesTable <- merge(totalDatesTable,totaltaxo, by = "OTU_Id")
  ##
  write.table(totalDatesTable, file = "Analyse-Composition-Rarefy/TableOnly/TableTotalDates.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(totalDatesTOTU, file = "Analyse-Composition-Rarefy/TableOnly/TotalDatesOTU.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(totalDatesTSequence, file = "Analyse-Composition-Rarefy/TableOnly/TotalDatesSequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  

# Hist OTU majoritaire ----------------------------------------------------
  #Cycle
    #label
  totalCyclesMajor <- totalCyclesTable
  totalCyclesTableJour <- totalCyclesMajor %>% filter(Cycles == "Jour")
  totalCyclesTableJour <- totalCyclesTableJour %>%
    mutate(lab.ypos = cumsum(value) - 0.5*value)
  totalCyclesTableNuit <- totalCyclesMajor %>% filter(Cycles == "Nuit")
  totalCyclesTableNuit <- totalCyclesTableNuit %>%
    mutate(lab.ypos = cumsum(value) - 0.5*value)
  totalCyclesMajor <- rbind(totalCyclesTableJour,totalCyclesTableNuit)
  kl <- ggplot(totalCyclesMajor, mapping = aes(y= value, x = Cycles, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity", color = "black") +
    scale_fill_manual(values = palette) + geom_label(aes(y = lab.ypos,label = OTU_Id),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
  kl <- kl + labs(x="Dates",y="Séquences") + theme(legend.position = "none")
  print(kl) 
  
  #FractionO
    #label
  totalFractionOMajor <- totalFractionOTable
  totalFractionOTableOxique <- totalFractionOMajor %>% filter(FractionO == "Oxique")
  totalFractionOTableOxique <- totalFractionOTableOxique %>%
    mutate(lab.ypos = cumsum(value) - 0.5*value)
  totalFractionOTableAnoxique <- totalFractionOMajor %>% filter(FractionO == "Anoxique")
  totalFractionOTableAnoxique <- totalFractionOTableAnoxique %>%
    mutate(lab.ypos = cumsum(value) - 0.5*value)
  totalFractionOMajor <- rbind(totalFractionOTableOxique,totalFractionOTableAnoxique)
  km <- ggplot(totalFractionOMajor, mapping = aes(y= value, x = FractionO, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity", color = "black") +
    scale_fill_manual(values = palette) + geom_label(aes(y = lab.ypos,label = OTU_Id),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
  km <- km + labs(x="Zones",y="Séquences") + theme(legend.position = "none")
  print(km) 
  
  #FractionT
  #label
  totalFractionTMajor <- totalFractionTTable
  totalFractionTTablePetite <- totalFractionTMajor %>% filter(FractionT == "Petite")
  totalFractionTTablePetite <- totalFractionTTablePetite %>%
    mutate(lab.ypos = cumsum(value) - 0.5*value)
  totalFractionTTableGrande <- totalFractionTMajor %>% filter(FractionT == "Grande")
  totalFractionTTableGrande <- totalFractionTTableGrande %>%
    mutate(lab.ypos = cumsum(value) - 0.5*value)
  totalFractionTMajor <- rbind(totalFractionTTablePetite,totalFractionTTableGrande)
  kn <- ggplot(totalFractionTMajor, mapping = aes(y= value, x = FractionT, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity", color = "black") +
    scale_fill_manual(values = palette) + geom_label(aes(y = lab.ypos,label = OTU_Id),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
  kn <- kn + labs(x="Fraction",y="Séquences") + theme(legend.position = "none")
  print(kn) 
  
  #Dates
  #label
  totalDatesMajor <- totalDatesTable
  totalDatesMajor$Dates <- as.character(totalDatesMajor$Dates)
  totalDatesTable04 <- totalDatesMajor %>% filter(Dates == "4")
  totalDatesTable04 <- totalDatesTable04 %>%
    mutate(lab.ypos = cumsum(value) - 0.5*value)
  totalDatesTable06 <- totalDatesMajor %>% filter(Dates == "6")
  totalDatesTable06 <- totalDatesTable06 %>%
    mutate(lab.ypos = cumsum(value) - 0.5*value)
  totalDatesTable09 <- totalDatesMajor %>% filter(Dates == "9")
  totalDatesTable09 <- totalDatesTable09 %>%
    mutate(lab.ypos = cumsum(value) - 0.5*value)
  totalDatesTable11 <- totalDatesMajor %>% filter(Dates == "11")
  totalDatesTable11 <- totalDatesTable11 %>%
    mutate(lab.ypos = cumsum(value) - 0.5*value)
  totalDatesMajor <- rbind(totalDatesTable04,totalDatesTable06,totalDatesTable09,totalDatesTable11)
  ko <- ggplot(totalDatesMajor, mapping = aes(y= value, x = Dates, fill = Division, group = Division), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity", color = "black") +
    scale_fill_manual(values = palette) + geom_label(aes(y = lab.ypos,label = OTU_Id),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
  ko <- ko + labs(x="Dates",y="Séquences")
  print(ko) 

# Table GENERAL + Parasitisme ------------------------------------------------------
  # Sequence ----------------------------------------------------------------
  Sequence_pool_Parasitism_Norm <- raw_Sequence_Norm
  pattern <- c(grep(pattern = "OSTA", colnames(Sequence_pool_Parasitism_Norm), value = FALSE, fixed = FALSE))
  Sequence_pool_Parasitism_Norm$Total<-rep(0, each = 3402)
  taxo <- tableV4095 %>% select(OTU_Id,NN.taxonomy)
  for (i in row.names(Sequence_pool_Parasitism_Norm)) {Sequence_pool_Parasitism_Norm[i,"Total"] <- sum(Sequence_pool_Parasitism_Norm[i,pattern])}
  w <- 0
  Sequence_pool_Parasitism_Norm <- merge(Sequence_pool_Parasitism_Norm,taxo, by = "OTU_Id")
  colnames(Sequence_pool_Parasitism_Norm)[35] <- "Taxonomy"
  #Total No-norm
  Parasitism_Norm_Sequence <- separate(Sequence_pool_Parasitism_Norm, Taxonomy, c('Domain', 'Supergroup','Division'), sep = ";")
  Parasitism_Norm_Sequence$Taxonomy <- paste(Parasitism_Norm_Sequence$Domain,Parasitism_Norm_Sequence$Supergroup,Parasitism_Norm_Sequence$Division, sep = "/")
  Parasitism_Norm_Sequencex <- Parasitism_Norm_Sequence %>% select(-Division,-Supergroup,-Domain,-OTU_Id)
  #Total Figure
  Parasitism_Norm_Sequencex <- Parasitism_Norm_Sequencex %>% group_by(Taxonomy) %>% summarise_all(funs(sum))
  Parasitism_Norm_Sequencex <- separate(Parasitism_Norm_Sequencex, Taxonomy, c('Domain', 'Supergroup','Division'), sep = "/")
  # Rapport sur 100
  pattern <- c(grep(pattern = c("OSTA"), colnames(Parasitism_Norm_Sequencex), value = TRUE, fixed = FALSE))
  pattern <- c(pattern, "Total")
  for (i in pattern) {Parasitism_Norm_Sequencex[,i]<-Parasitism_Norm_Sequencex[,i]*100/sum(Parasitism_Norm_Sequencex[,i])}
  #Add Column Parasitism
  Parasitism_Norm_Sequencex$Parasitism <- rep("Autres",each = nrow(Parasitism_Norm_Sequencex))
  Parasitemin <- c("Apicomplexa","Perkinsea","Chytridiomycota","Cryptomycota","Ichthyophonida","Cercozoa","Oomycetes","Pirsonia","Labyrinthulida")
  Parasitemax <- c("Apicomplexa","Perkinsea","Hexamitidae","Kinetoplastida","Chytridiomycota","Cryptomycota","Ichthyophonida","Trichomonadida","Cercozoa","Oomycetes","Pirsonia","Dinophyceae","Centramoebida","Tubulinea","Rhizophydium","Zygomycota_1","Florideophyceae","Labyrinthulida")
  w <- 0
  for (j in Parasitism_Norm_Sequencex$Division) { w <- w+1
  for (i in Parasitemin) {
    if (i == j) { Parasitism_Norm_Sequencex[w,"Parasitism"] <- "Parasite" }}}
  #Melting 
  Parasitism_Norm_Sequencex <- as.data.frame(melt(Parasitism_Norm_Sequencex, id = c("Domain","Supergroup","Division","Parasitism")))
  Parasitism_polarNorm_Sequence <- Parasitism_Norm_Sequencex %>% filter(variable == "Total")
  Parasitism_polarNorm_Sequence$Division <- paste(Parasitism_polarNorm_Sequence$Parasitism,Parasitism_polarNorm_Sequence$Division, sep = ";")
  for (i in row.names(Parasitism_polarNorm_Sequence)) { 
    if (Parasitism_polarNorm_Sequence[i,"Parasitism"] == "Autres") { Parasitism_polarNorm_Sequence[i,"Division"] <- "Autres"}}
  Parasitism_polarNorm_Sequence <- Parasitism_polarNorm_Sequence %>% select(-Domain,-variable,-Supergroup,-Parasitism)
  Parasitism_polarNorm_Sequence <- Parasitism_polarNorm_Sequence %>% group_by(Division) %>% summarise_all(funs(sum))
  Parasitism_polarNorm_Sequence <- separate(Parasitism_polarNorm_Sequence, Division, c('Division','Parasitism'))
  for (i in row.names(Parasitism_polarNorm_Sequence)) {
    if (Parasitism_polarNorm_Sequence[i,"Division"] == "Autres") { Parasitism_polarNorm_Sequence[i,"Parasitism"] <- "Autres"}}
  #Label
  Parasitism_polarNorm_Sequence <- Parasitism_polarNorm_Sequence %>%
    arrange(desc(Parasitism)) %>%
    mutate(lab.ypos = cumsum(value) - 0.5*value)
  Parasitism_polarNorm_Sequence$label <- paste(round(Parasitism_polarNorm_Sequence$value,1), "%", sep = "")
  for (i in rownames(Parasitism_polarNorm_Sequence)) {
    if (Parasitism_polarNorm_Sequence[i,"label"] == "0%") { Parasitism_polarNorm_Sequence[i,"label"] <- NA}}
  for (i in rownames(Parasitism_polarNorm_Sequence)) {
    if (is.na(Parasitism_polarNorm_Sequence[i,"label"]) == FALSE) { Parasitism_polarNorm_Sequence[i,"label"] <- paste(Parasitism_polarNorm_Sequence[i,"Parasitism"]," : ",Parasitism_polarNorm_Sequence[i,"label"], sep = "")}}
  #Figure
  pdf("Analyse-Composition-Rarefy/Composition/Analyse-polar-Sequence-Parasitism(min)-Detail-Total-Norm-V4_095.pdf",width = 8.00,height = 8.00)
  ax <- ggplot(Parasitism_polarNorm_Sequence, mapping = aes(y= value, x = 2, fill = Parasitism), Rowv = NA, col = colMain, scale = "column") +
    geom_bar(stat="identity", color = "white", width = 1) + coord_polar("y") + 
    geom_label_repel(aes(y = lab.ypos,label = label),color = "white",segment.color = "black",show.legend = FALSE, nudge_x = 0.5) +
    scale_y_continuous(limits=c(0,sum(Parasitism_polarNorm_Sequence %>% select(value))))+
    xlim(1,2.5) +
    theme_unique_darkbis() + 
    #facet_wrap( ~ variable , nrow = 4) +
    labs(title="Part de Parasitism en Séquences Totale Norm",y = "Total",x="")
  print(ax)
  dev.off()
  
  #palette <- c(pal_locuszoom(alpha = 0.9)(5), pal_lancet(alpha = 0.9)(5))
  pdf("Analyse-Composition-Rarefy/Composition/Analyse-polar-Sequence-Parasitism(min)-Detail-Total-Norm-V4_095.pdf",width = 4.00,height = 4.00)
  ax <- ggplot(Parasitism_polarNorm_Sequence, mapping = aes(y= value, x = 2, fill = Parasitism), Rowv = NA, col = colMain, scale = "column") +
    geom_bar(stat="identity", color = "white", width = 1) + coord_polar("y") + 
    geom_label_repel(aes(y = lab.ypos,label = label),color = "white",size = 3,segment.color = "black",show.legend = FALSE, nudge_x = 0.5) +
    scale_y_continuous(limits=c(0,sum(Parasitism_polarNorm_Sequence %>% select(value))))+
    xlim(1,2.5) +
    theme_unique_darkbis() + 
    #facet_wrap( ~ variable , nrow = 4) +
    labs(y = "Séquences",x="") + scale_fill_manual(values = rev(palette))
  print(ax)
  dev.off()
  
  svglite("Analyse-Composition-Rarefy/Composition/Analyse-polar-Sequence-Parasitism(min)-Detail-Total-Norm-V4_095.svg",width = 4.00,height = 4.00)
  ax <- ggplot(Parasitism_polarNorm_Sequence, mapping = aes(y= value, x = 2, fill = Parasitism), Rowv = NA, col = colMain, scale = "column") +
    geom_bar(stat="identity", color = "white", width = 1) + coord_polar("y") + 
    geom_label_repel(aes(y = lab.ypos,label = label),color = "white",size = 3,segment.color = "black",show.legend = FALSE, nudge_x = 0.5) +
    scale_y_continuous(limits=c(0,sum(Parasitism_polarNorm_Sequence %>% select(value))))+
    xlim(1,2.5) +
    theme_unique_darkbis() + 
    #facet_wrap( ~ variable , nrow = 4) +
    labs(y = "Séquences",x="") + scale_fill_manual(values = rev(palette))
  print(ax)
  dev.off()
  # OTU ----------------------------------------------------------------
  OTU_pool_Parasitism_Norm <- raw_OTU_Norm
  pattern <- c(grep(pattern = "OSTA", colnames(OTU_pool_Parasitism_Norm), value = FALSE, fixed = FALSE))
  OTU_pool_Parasitism_Norm$Total<-rep(0, each = 3402)
  taxo <- tableV4095 %>% select(OTU_Id,NN.taxonomy)
  for (i in row.names(OTU_pool_Parasitism_Norm)) {if (sum(OTU_pool_Parasitism_Norm[i,pattern]) > 0) { OTU_pool_Parasitism_Norm[i,"Total"] <- 1}}
  w <- 0
  OTU_pool_Parasitism_Norm <- merge(OTU_pool_Parasitism_Norm,taxo, by = "OTU_Id")
  colnames(OTU_pool_Parasitism_Norm)[35] <- "Taxonomy"
  #Total No-norm
  Parasitism_Norm_OTU <- separate(OTU_pool_Parasitism_Norm, Taxonomy, c('Domain', 'Supergroup','Division'), sep = ";")
  Parasitism_Norm_OTU$Taxonomy <- paste(Parasitism_Norm_OTU$Domain,Parasitism_Norm_OTU$Supergroup,Parasitism_Norm_OTU$Division, sep = "/")
  Parasitism_Norm_OTUx <- Parasitism_Norm_OTU %>% select(-Division,-Supergroup,-Domain,-OTU_Id)
  #Total Figure
  Parasitism_Norm_OTUx <- Parasitism_Norm_OTUx %>% group_by(Taxonomy) %>% summarise_all(funs(sum))
  Parasitism_Norm_OTUx <- separate(Parasitism_Norm_OTUx, Taxonomy, c('Domain', 'Supergroup','Division'), sep = "/")
  # Rapport sur 100
  pattern <- c(grep(pattern = c("OSTA"), colnames(Parasitism_Norm_OTUx), value = TRUE, fixed = FALSE))
  pattern <- c(pattern, "Total")
  for (i in pattern) {Parasitism_Norm_OTUx[,i]<-Parasitism_Norm_OTUx[,i]*100/sum(Parasitism_Norm_OTUx[,i])}
  #Add Column Parasitism
  Parasitism_Norm_OTUx$Parasitism <- rep("Autres",each = nrow(Parasitism_Norm_OTUx))
  Parasitemin <- c("Apicomplexa","Perkinsea","Chytridiomycota","Cryptomycota","Ichthyophonida","Cercozoa","Oomycetes","Pirsonia","Labyrinthulida")
  Parasitemax <- c("Apicomplexa","Perkinsea","Hexamitidae","Kinetoplastida","Chytridiomycota","Cryptomycota","Ichthyophonida","Trichomonadida","Cercozoa","Oomycetes","Pirsonia","Dinophyceae","Centramoebida","Tubulinea","Rhizophydium","Zygomycota_1","Florideophyceae","Labyrinthulida")
  w <- 0
  for (j in Parasitism_Norm_OTUx$Division) { w <- w+1
  for (i in Parasitemin) {
    if (i == j) { Parasitism_Norm_OTUx[w,"Parasitism"] <- "Parasite" }}}
  #Melting 
  Parasitism_Norm_OTUx <- as.data.frame(melt(Parasitism_Norm_OTUx, id = c("Domain","Supergroup","Division","Parasitism")))
  Parasitism_polarNorm_OTU <- Parasitism_Norm_OTUx %>% filter(variable == "Total")
  Parasitism_polarNorm_OTU$Division <- paste(Parasitism_polarNorm_OTU$Parasitism,Parasitism_polarNorm_OTU$Division, sep = ";")
  for (i in row.names(Parasitism_polarNorm_OTU)) { 
    if (Parasitism_polarNorm_OTU[i,"Parasitism"] == "Autres") { Parasitism_polarNorm_OTU[i,"Division"] <- "Autres"}}
  Parasitism_polarNorm_OTU <- Parasitism_polarNorm_OTU %>% select(-Domain,-variable,-Supergroup,-Parasitism)
  Parasitism_polarNorm_OTU <- Parasitism_polarNorm_OTU %>% group_by(Division) %>% summarise_all(funs(sum))
  Parasitism_polarNorm_OTU <- separate(Parasitism_polarNorm_OTU, Division, c('Division','Parasitism'))
  for (i in row.names(Parasitism_polarNorm_OTU)) {
    if (Parasitism_polarNorm_OTU[i,"Division"] == "Autres") { Parasitism_polarNorm_OTU[i,"Parasitism"] <- "Autres"}}
  #Label
  Parasitism_polarNorm_OTU <- Parasitism_polarNorm_OTU %>%
    arrange(desc(Parasitism)) %>%
    mutate(lab.ypos = cumsum(value) - 0.5*value)
  Parasitism_polarNorm_OTU$label <- paste(round(Parasitism_polarNorm_OTU$value,1), "%", sep = "")
  for (i in rownames(Parasitism_polarNorm_OTU)) {
    if (Parasitism_polarNorm_OTU[i,"label"] == "0%") { Parasitism_polarNorm_OTU[i,"label"] <- NA}}
  for (i in rownames(Parasitism_polarNorm_OTU)) {
    if (is.na(Parasitism_polarNorm_OTU[i,"label"]) == FALSE) { Parasitism_polarNorm_OTU[i,"label"] <- paste(Parasitism_polarNorm_OTU[i,"Parasitism"]," : ",Parasitism_polarNorm_OTU[i,"label"], sep = "")}}
  #Figure
  pdf("Analyse-Composition-Rarefy/Composition/Analyse-polar-OTU-Parasitism(min)-Detail-Total-Norm-V4_095.pdf",width = 8.00,height = 8.00)
  ax <- ggplot(Parasitism_polarNorm_OTU, mapping = aes(y= value, x = 2, fill = Parasitism), Rowv = NA, col = colMain, scale = "column") +
    geom_bar(stat="identity", color = "white", width = 1) + coord_polar("y") + 
    geom_label_repel(aes(y = lab.ypos,label = label),color = "white",segment.color = "black",show.legend = FALSE, nudge_x = 0.5) +
    scale_y_continuous(limits=c(0,sum(Parasitism_polarNorm_OTU %>% select(value))))+
    xlim(1,2.5) +
    theme_unique_darkbis() + 
    #facet_wrap( ~ variable , nrow = 4) +
    labs(title="Part de Parasitism en OTUs Totale Norm",y = "Total",x="")
  print(ax)
  dev.off()
  
  pdf("Analyse-Composition-Rarefy/Composition/Analyse-polar-OTU-Parasitism(min)-Detail-Total-Norm-V4_095.pdf",width = 4.00,height = 4.00)
  ax <- ggplot(Parasitism_polarNorm_OTU, mapping = aes(y= value, x = 2, fill = Parasitism), Rowv = NA, col = colMain, scale = "column") +
    geom_bar(stat="identity", color = "white", width = 1) + coord_polar("y") + 
    geom_label_repel(aes(y = lab.ypos,label = label),size = 3, color = "white",segment.color = "black",show.legend = FALSE, nudge_x = 0.5) +
    scale_y_continuous(limits=c(0,sum(Parasitism_polarNorm_OTU %>% select(value))))+
    xlim(1,2.5) +
    theme_unique_darkbis() + 
    #facet_wrap( ~ variable , nrow = 4) +
    labs(y = "OTUs",x="") + scale_fill_manual(values = rev(palette))
  print(ax)
  dev.off()
  
  svglite("Analyse-Composition-Rarefy/Composition/Analyse-polar-OTU-Parasitism(min)-Detail-Total-Norm-V4_095.svg",width = 4.00,height = 4.00)
  ax <- ggplot(Parasitism_polarNorm_OTU, mapping = aes(y= value, x = 2, fill = Parasitism), Rowv = NA, col = colMain, scale = "column") +
    geom_bar(stat="identity", color = "white", width = 1) + coord_polar("y") + 
    geom_label_repel(aes(y = lab.ypos,label = label),size =3, color = "white",segment.color = "black",show.legend = FALSE, nudge_x = 0.5) +
    scale_y_continuous(limits=c(0,sum(Parasitism_polarNorm_OTU %>% select(value))))+
    xlim(1,2.5) +
    theme_unique_darkbis() + 
    #facet_wrap( ~ variable , nrow = 4) +
    labs(y = "OTUs",x="") + scale_fill_manual(values = rev(palette))
  print(ax)
  dev.off()
  
  #Ad grid : grid.locator(unit="native") 
  #grid.locator(unit="native") 
  #grid.brackets(208,174,133.875,247.5,h = 0.01,lwd=1.5)
  
# Sort Percent no affiliate
  #Best Hit %
  totaltaxofidel <- totaltaxo %>% filter(Identity.... > 95.0)
  nrow(totaltaxofidel)
  #NN.taxonomy
    #Division
     totaltaxosupergroup <- totaltaxo %>% select(OTU_Id,NN.taxonomy)
     totaltaxosupergroup <- separate(totaltaxosupergroup, NN.taxonomy, c('Domain', 'Supergroup','Division','Class','Order'), sep = ";")
     DivPattern <- c("unclassified Stramenopiles","unclassified Fungi","unclassified Choanoflagellida","Fungi incertae sedis","environmental samples")
     nrow(totaltaxosupergroup %>% filter(Division == "unclassified Stramenopiles"))
     nrow(totaltaxosupergroup %>% filter(Division == "unclassified Fungi"))
     nrow(totaltaxosupergroup %>% filter(Division == "unclassified Choanoflagellida"))
     nrow(totaltaxosupergroup %>% filter(Division == "Fungi incertae sedis"))
     nrow(totaltaxosupergroup %>% filter(Division == "environmental samples"))
     #Class
     nrow(totaltaxosupergroup %>% filter(Class == ""))
     nrow(totaltaxosupergroup %>% filter(Class == "Basal fungal lineages"))
     nrow(totaltaxosupergroup %>% filter(Class == "Bicosoecida clade"))
     nrow(totaltaxosupergroup %>% filter(Class == "environmental samples"))
     nrow(totaltaxosupergroup %>% filter(Class == "unclassified Bicosoecida"))
     nrow(totaltaxosupergroup %>% filter(Class == "unclassified Chrysophyceae"))
     nrow(totaltaxosupergroup %>% filter(Class == "unclassified oomycetes"))
     nrow(totaltaxosupergroup %>% filter(Class == "unclassified Perkinsea"))
     #Order
     nrow(totaltaxosupergroup %>% filter(Order == ""))
     nrow(totaltaxosupergroup %>% filter(is.na(Order)))
     nrow(totaltaxosupergroup %>% filter(Order == "Ascomycota incertae sedis"))
     nrow(totaltaxosupergroup %>% filter(Order == "environmental samples"))
     nrow(totaltaxosupergroup %>% filter(Order == "environvironmental samples"))
     nrow(totaltaxosupergroup %>% filter(Order == "unclassified Thraustochytriidae"))

     
     
     

# AFC FractionO avec Dates ------------------------------------------------
     #Figure FractionO
     dataSequenceTaxox <- dataSequenceTaxo
     #Avril
     for (i in row.names(dataSequenceTaxox)) { 
       if (dataSequenceTaxox[i,"FractionO"] != "Communs"){ 
         if (dataSequenceTaxox[i,"Total04"] == 1) { (dataSequenceTaxox[i,"Total04"] <- "Avril")}
         if (dataSequenceTaxox[i,"Total06"] == 1) { (dataSequenceTaxox[i,"Total06"] <- "Juin")}
         if (dataSequenceTaxox[i,"Total09"] == 1) { (dataSequenceTaxox[i,"Total09"] <- "Septembre")}
         if (dataSequenceTaxox[i,"Total11"] == 1) { (dataSequenceTaxox[i,"Total11"] <- "Novembre")}}}
     for (i in row.names(dataSequenceTaxox)) { if (dataSequenceTaxox[i,"Total04"] != "Avril") { dataSequenceTaxox[i,"Total04"] <- "No"}}
     for (i in row.names(dataSequenceTaxox)) { if (dataSequenceTaxox[i,"Total06"] != "Juin") { dataSequenceTaxox[i,"Total06"] <- "No"}}
     for (i in row.names(dataSequenceTaxox)) { if (dataSequenceTaxox[i,"Total09"] != "Septembre") { dataSequenceTaxox[i,"Total09"] <- "No"}}
     for (i in row.names(dataSequenceTaxox)) { if (dataSequenceTaxox[i,"Total11"] != "Novembre") { dataSequenceTaxox[i,"Total11"] <- "No"}}
     
     #Print plot
     bs <- ggplot(dataSequenceTaxox, aes(y = `Dim 2`, x = `Dim 1`, color = FractionO)) + geom_point(size = 2) + #geom_text(aes(label=sample), size = 3, hjust=1.2, vjust=1) +
       geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
       geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
       #facet_grid(~`Dates`, switch = "x") + theme_bw() +
       #stat_ellipse(aes(group = Cycles),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "red") +
       #stat_ellipse(aes(linetype = FractionO, alpha = FractionO),geom = "polygon",type = "norm") +
       stat_ellipse(aes(linetype = Total04,color = Total04, alpha = Total04),geom = "polygon",type = "norm") +
       stat_ellipse(aes(linetype = Total06,color = Total06, alpha = Total06),geom = "polygon",type = "norm") +
       stat_ellipse(aes(linetype = Total09,color = Total09, alpha = Total09),geom = "polygon",type = "norm") +
       stat_ellipse(aes(linetype = Total11,color = Total11, alpha = Total11),geom = "polygon",type = "norm") +
       scale_linetype_manual(values = c(2,2,0,2,2)) +
       #stat_ellipse(geom = "polygon",type = "norm", linetype = 2, alpha = 0.1) +
       #stat_ellipse(aes(group = `Fraction-Oxygène`),geom = "polygon",type = "norm", linetype = 2, alpha = 0.1, color = "grey") +
       #geom_label_repel(aes(label = Cycles), color = 'white',size = 2.5, segment.size = 0.3,segment.color = "black", alpha = 0.8) +
       theme(axis.title = element_text(face="bold", size=12), 
             axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
             title = element_text(face="bold", size=14),
             legend.title = element_text(face="bold"),
             legend.position = "right",
             legend.text = element_text(size=10)) +
       scale_alpha_manual(values = c(0.1,0.1,0,0.1,0.1)) +
       scale_color_manual(values = c("#F8766D","#00A5FF","lightgrey","lightgreen","black","red","blue","green")) +
       labs(x=Dim1Seq,y=Dim2Seq,color = "Dates", alpha = "Dates", linetype = "Dates", shape = "Zones")
     print(bs)
     
     
     