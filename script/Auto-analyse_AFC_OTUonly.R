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
#if (!require("pacman")) {install.packages("pacman", repos="http://cran.rstudio.com/")}
#pacman::p_load("ggplot2", "readxl","dplyr","tidyr","cowplot","FactoMineR","factoextra","reshape2","varhandle","gplots","ggrepel","ggpubr","ggsci","scales","hrbrthemes","GUniFrac","svglite")

pkg <- c("ggplot2", "readxl","dplyr","tidyr","cowplot","FactoMineR","factoextra","reshape2","varhandle","ggrepel","ggpubr","ggsci","scales","hrbrthemes","GUniFrac","svglite")
lapply(pkg, require, character.only = TRUE)

system("mkdir Analyse-Composition-Rarefy")
system("mkdir Analyse-Composition-Rarefy/Figure-Sum")
system("mkdir Analyse-Composition-Rarefy/AFC-Duplicat")
system("mkdir Analyse-Composition-Rarefy/Composition")
system("mkdir Analyse-Composition-Rarefy/HistOnly")
system("mkdir Analyse-Composition-Rarefy/Figure-Polar")
system("mkdir Analyse-Composition-Rarefy/AFC-Parasite")
system("mkdir Analyse-Composition-Rarefy/TableOnly")
system("mkdir Analyse-Composition-Rarefy/Biplot")

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
colnames(samples_df) <- c("sample","Technologie","Regions","ADN","Cycle","Fraction-Oxygène","Fraction-Taille","Dates","Replicats","Conditions")
w <- 0
for (i in samples_df$Cycle) { w <- w +1
if (i == "J") { samples_df[w,"Cycle"] <- "Jour"}
if (i == "N") { samples_df[w,"Cycle"] <- "Nuit"}}
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
  # Rarefy ---------------------------------------------------------------------
pattern <- c(grep(pattern = "OSTA", colnames(tableV4095), value = FALSE, fixed = FALSE))
otu_mat <- tableV4095 %>% select(OTU_Id,pattern)
row.names(otu_mat) <- otu_mat$OTU_Id
otu_mat <- otu_mat %>% select(-OTU_Id)
otu_mat[,] <- lapply(otu_mat[,], as.numeric)
otu_matx <- t(otu_mat)
sequence_mat_rare <- t(Rarefy(otu_matx)$otu.tab.rff)


  # PA-AB ------------------------------------------------------------------
otu_mat_rare<- sequence_mat_rare
pattern <- c(grep(pattern = "OSTA", colnames(otu_mat_rare), value = FALSE, fixed = FALSE))
otu_mat_rare[,pattern][otu_mat_rare[,pattern] != 0] <- 1

#Stat Rarefy
  #Sequence
  avRarefyS <- as.data.frame(colSums(otu_mat))
  apRarefyS <- as.data.frame(colSums(sequence_mat_rare))
  statRarefy <- cbind(avRarefyS,apRarefyS)
  #OTU
  otu_mat_uniq <- otu_mat
  otu_mat_uniq[otu_mat_uniq != 0] <- 1
  avRarefyO <- as.data.frame(colSums(otu_mat_uniq))
  apRarefyO <- as.data.frame(colSums(otu_mat_rare))
  statRarefy <- cbind(statRarefy,avRarefyO,apRarefyO)
  statRarefy["Total",]<-colSums(statRarefy)
  colnames(statRarefy) <- c("avRarefy-Sequence","apRarefy-Sequence","avRarefy-OTU","apRarefy-OTU")
  write.table(statRarefy, file = "Analyse-Composition-Rarefy/TableOnly/StatRarefy_withDuplicat.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
  
  
# OTU ---------------------------------------------------------------------
  # AFC Plot ----------------------------------------------------------------
  dt <- as.data.frame(otu_mat_rare)
  dt <- dt[,-1]  
  res.ca <- CA (dt, graph = FALSE,ncp = 2 )
  pdf("Analyse-Composition-Rarefy/AFC-Duplicat/AFC_OTU-V4_095.pdf",width = 16.00,height = 9.00)
  fviz_ca_col(res.ca, repel = TRUE, col.col = "lightgreen") + theme_unique_dark()
  dev.off()
  p <- get_ca_col(res.ca)
  coord <- p$coord
  coord <- as.data.frame(coord)
  coord$sample <- row.names(coord)
  data <- merge(x = coord, y = samples_df, by= "sample")
  fviz_screeplot (res.ca, addlabels = TRUE, ylim = c(0, 50))
  
  Xseq <- fviz_screeplot (res.ca, addlabels = TRUE, ylim = c(0, 50))
  Xseq <- Xseq$data
  Dim1Seq <- paste("Dim 1 [",round(Xseq %>% filter(dim == 1) %>% select(eig),1),"%]", sep = "")
  Dim2Seq <- paste("Dim 2 [",round(Xseq %>% filter(dim == 2) %>% select(eig),1),"%]", sep = "")
  
  # Plot AFC ----------------------------------------------------------------
    # Overview* -------------------------------------------------------------------
  pdf("Analyse-Composition-Rarefy/AFC-Duplicat/Analyse-OTU-V4_095.pdf",width = 16.00,height = 9.00)
  a <- ggplot(data, aes(y = `Dim 1`, x = `Dim 2`)) + geom_point(aes(color = Conditions, shape = Dates), size = 2) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    facet_grid(~Dates, switch = "x") + theme_bw() +
    geom_label_repel(aes(label = sample,fill = Conditions), color = 'white',size = 2.5, segment.size = 0,segment.color = "black", alpha = 0.8) +
    geom_polygon(aes(color = Conditions)) +
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
  
  svglite("Analyse-Composition-Rarefy/AFC-Duplicat/Analyse-OTU-11-V4_095.svg",width = 8.00,height = 6.00)
  a <- ggplot(data %>% filter(Dates == "11"), aes(y = `Dim 1`, x = `Dim 2`)) + geom_point(aes(color = Conditions, shape = Dates), size = 3) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    geom_label_repel(aes(label = sample,fill = Conditions), color = 'white',size = 4, segment.size = 0.5,segment.color = "black", alpha = 0.8) +
    #geom_polygon(aes(color = Conditions)) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x=Dim1Seq,y=Dim2Seq) + xlim(-1,1) +
    guides(shape = FALSE, color = FALSE) + labs(color = "Conditions")
  print(a)
  dev.off()
  

# Sequence ---------------------------------------------------------------------
  # AFC Plot ----------------------------------------------------------------
  dt <- as.data.frame(sequence_mat_rare)
  dt <- dt[,-1]  
  res.ca <- CA (dt, graph = FALSE,ncp = 2 )
  pdf("Analyse-Composition-Rarefy/AFC-Duplicat/AFC_Seq-V4_095.pdf",width = 16.00,height = 9.00)
  fviz_ca_col(res.ca, repel = TRUE, col.col = "lightgreen") + theme_unique_dark()
  dev.off()
  p <- get_ca_col(res.ca)
  coord <- p$coord
  coord <- as.data.frame(coord)
  coord$sample <- row.names(coord)
  data <- merge(x = coord, y = samples_df, by= "sample")
  fviz_screeplot (res.ca, addlabels = TRUE, ylim = c(0, 50))
  
  Xseq <- fviz_screeplot (res.ca, addlabels = TRUE, ylim = c(0, 50))
  Xseq <- Xseq$data
  Dim1Seq <- paste("Dim 1 [",round(Xseq %>% filter(dim == 1) %>% select(eig),1),"%]", sep = "")
  Dim2Seq <- paste("Dim 2 [",round(Xseq %>% filter(dim == 2) %>% select(eig),1),"%]", sep = "")
  
  # Plot AFC ----------------------------------------------------------------
    # Overview* -------------------------------------------------------------------
  pdf("Analyse-Composition-Rarefy/AFC-Duplicat/Analyse-Sequence-V4_095.pdf",width = 16.00,height = 9.00)
  a <- ggplot(data, aes(y = `Dim 1`, x = `Dim 2`)) + geom_point(aes(color = Conditions, shape = Dates), size = 2) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    facet_grid(~Dates, switch = "x") + theme_bw() +
    geom_label_repel(aes(label = sample,fill = Conditions), color = 'white',size = 2.5, segment.size = 0,segment.color = "black", alpha = 0.8) +
    geom_polygon(aes(color = Conditions)) +
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
  
  svglite("Analyse-Composition-Rarefy/AFC-Duplicat/Analyse-Sequence-11-V4_095.svg",width = 8.00,height = 6.00)
  a <- ggplot(data %>% filter(Dates == "11"), aes(y = `Dim 1`, x = `Dim 2`)) + geom_point(aes(color = Conditions, shape = Dates), size = 3) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    geom_label_repel(aes(label = sample,fill = Conditions), color = 'white',size = 4, segment.size = 0.5,segment.color = "black", alpha = 0.8) +
    #geom_polygon(aes(color = Conditions)) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x=Dim1Seq,y=Dim2Seq) + xlim(-1,1) +
    guides(shape = FALSE, color = FALSE) + labs(color = "Conditions")
  print(a)
  dev.off()
  
  
