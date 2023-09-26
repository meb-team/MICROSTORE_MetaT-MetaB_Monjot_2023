# !/usr/bin/env Rscript
#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 09/06/2022
#
# Script Metatrans_Table analyze
Cstack_info()["size"]
# Set input and output -----------------------------------------------------------
# Detect R or Rstudio
se <- Sys.getenv()
if (length(se[grepl("rstudio",se,ignore.case=TRUE)]) > 0 ) { inputmode <- TRUE }
if (length(se[grepl("rstudio",se,ignore.case=TRUE)]) == 0 ) { inputmode <- FALSE }
#
# Input argument if using Rstudio
if (inputmode == TRUE) {
input <- "main_table.mapping.unique.raw.noHuman.noConta.noMetazoa.annot.tsv"
output <- "V4-unified-correct-paired-out-compo"
tax <- "table_taxonomy.perUnigene.allUnigenes.tsv"
pr2 <- "pr2_version_4.14.0_SSU_dada2.fasta.gz"
}
#
# Input argument if using R
args = commandArgs(trailingOnly=TRUE)
if ( inputmode == FALSE ) {
  input <- args[1]
  output <- args[2]
  tax <- args[3]
  pr2 <- args[4]
}
print(input)
print(tax)
print(output)
print(pr2)
alltranscript_analysis <- "no" # Warning, this option launch DESeq2 on the entire dataset but requires large amount of memory (RAM > 64 Go)
#
# Import package -----------------------------------------------------------
pkg <- c("ggplot2","dplyr","tidyr","cowplot","FactoMineR","factoextra","reshape2","varhandle","ggrepel","ggpubr","ggsci","scales","hrbrthemes","data.table","svglite","treemap", "VennDiagram","stringr","paletteer","elementalist","gtools","vegan","treemapify","gplots","car","scatterplot3d","DESeq2","SARTools","ComplexHeatmap","RColorBrewer","ggh4x")
lapply(pkg, require, character.only = TRUE)
#
paletspe <- rep(c("#e53a35","#f44336","#ef5350","#e57373","#ef9a9a","#f4d4d4",
                  "#7b1fa2","#8e24aa","#8e24aaFF",
                  "#303f9f","#303f9fFF",
                  "#3065c0","#3f88e5","#049be5FF","#64b5f6","#b8d9fd","#64b5f680",
                  "#2e7d32","#689f38",#"#7cb342",
                  "#fbc02c","#fdd835",#"#ffeb3a","#ffee58","#fff176",
                  "#c6c6c6",#"#f49800","#f57c0099","#f4980080",
                  "#616161","#757575","#9e9e9e","#bdbdbd",#"#e0e0e0",
                  "#f57c00",
                  "#1876d2","#1876d295",#"#4cc6da","#80deea","#b2ebf2",
                  "#3897a7","#41acc1","#4cc6da","#80deea","#b2ebf2"),
                4)
#
# Set directory, create file result -----------------------------------------------------------
if (inputmode == TRUE) {
  current <- dirname(rstudioapi::getSourceEditorContext()$path)
  setwd(current)
  output <- paste(current,"../dataDADA2/result",output,"Metatrans", sep = "/")
}
if (inputmode == FALSE) {
  output <- paste("../dataDADA2/result",output,"Metatrans", sep = "/")
}
# Set working directory
if (dir.exists(output) == FALSE) {dir.create(output)}
if (dir.exists(output) == TRUE) { setwd(output) }
# Create result files
dir.create("DESeq2")
dir.create("DESeq2/Summary")
dir.create("DESeq2/Heatmap_Interaction")
dir.create("DESeq2/Heatmap_raw")
dir.create("DESeq2/Quaternaryplot")
dir.create("DESeq2/PCoA")
dir.create("Stat")
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
theme_unique_darktris <- function (base_size = 12, base_family = "") {
  ret <- (theme_bw(base_size = base_size, base_family = base_family) +
            theme(text = element_text(colour = "black"),
                  title = element_text(color = "black", face = "bold", vjust = 0.5, hjust = 0.5,size = 10),
                  axis.ticks = element_blank(),
                  line = element_line(color = "black"),
                  rect = element_rect(fill = "white", color = "black"),
                  axis.title = element_text(color = "black", face = "bold"),
                  axis.text.y = element_blank(),
                  axis.text.x = element_blank(),
                  legend.background = element_rect(fill = NULL, color = NULL),
                  legend.key = element_rect(fill = NA, color = NA, linetype = 0),
                  legend.text = element_text(size = 12,face ="bold"),
                  legend.title = element_text(size = 0, face="bold"),
                  strip.background = element_rect(fill=NA,colour=NA,size=NA,linetype = NULL),
                  strip.text = element_text(color="black",face="bold",vjust=.5,hjust=.5,size = 14),
                  panel.background = element_rect(fill = "white", color = NULL),
                  panel.border = element_blank(),
                  panel.grid.major = element_line(color = "white"),
                  panel.grid.minor = element_line(color = "white"),
                  plot.background = element_rect(fill = "white", colour = "white", linetype = 0)))
  ret
} 
# Input Tables ---------------------------------------------------------
tableVinput <- fread(file = paste("../../../../rawdata",input,sep = "/"), sep = "\t")
# Taxo
Tax_table <- fread(file = paste("../../../../rawdata",tax,sep = "/"), sep = "\t")
# ko2_hierarachy
ko2Met_Table <- fread(file = "../../../../rawdata/ko_to_hierarchy.txt", sep = "\t")
# Taxonomy reference in metabarcoding
PR2_tax <- fread(file = paste("../../../../dataBase",pr2,sep = "/"), sep = "\t" , header = F, col.names = c("Taxonomy"))
PR2_tax <- PR2_tax %>% dplyr::filter(grepl("Eukaryota;",Taxonomy)==T)
# Functional group table 
funct2Last <- fread(file = "../Functional-Analyse/Table/ASV-Table-GrpFunct.csv", sep = "\t")
funct2Last <- as.data.frame(funct2Last)
funct2Last <- funct2Last %>% dplyr::select(Last, nbGroup) %>% distinct()

# Stat 1 --------------------------------------------------------------------
tableVinput$Somme <- rowSums(tableVinput %>% dplyr::select(-unigene, -PfamDom, -ko, -ko_BestHitNotSignif, -GOterms))
nrow(tableVinput)
nrow(tableVinput %>% filter(Somme != 0))
tableVinput0 <- tableVinput %>% dplyr::filter(Somme == 0)
tableVinput <- tableVinput %>% dplyr::select(-Somme)
# Prepare data inf --------------------------------------------------------
infdataini <- read.table(file = "../../../../rawdata/metadata_metaT", sep = ";", header = TRUE)
infdataini$Echantillon <- sapply(strsplit(as.character(infdataini$Echantillon), "_"), function(x) x[[2]])
infdataini$Rep <- "_1"

infdataini$Variable <- paste0(infdataini$Ref..collab,infdataini$Rep)
infdataini <- infdataini %>% dplyr::select(-"Rep")
infdataini <- separate(infdataini, Variable, c("Condition","Periods","Replicat"), sep = "_")
infdataini[,"Replicat"][infdataini[,"Replicat"] == "01"] <- "1"
infdataini[,"Replicat"][infdataini[,"Replicat"] == "02"] <- "2"
infdataini$Ref
## Prepare sample_df
infdataini <- separate(infdataini, Condition, c("empty","Cycle","Zone","Fraction"),sep = "")
infdataini <- infdataini %>% dplyr::select(-"empty")
infdataini<-as.data.frame(infdataini)
col <- c("Echantillon","Ref..collab","Cycle","Zone","Fraction","Periods","Replicat")
samples_df <- infdataini[,col]
samples_df$Condition <- paste(samples_df$Cycle,samples_df$Zone,samples_df$Fraction, sep = "")
for (i in row.names(samples_df)) {
if (samples_df[i,"Cycle"] == "J") { samples_df[i,"Cycle"] <- "Day"}
if (samples_df[i,"Cycle"] == "N") { samples_df[i,"Cycle"] <- "Night"}
if (samples_df[i,"Zone"] == "O") { samples_df[i,"Zone"] <- "Mixolimnion"}
if (samples_df[i,"Zone"] == "A") { samples_df[i,"Zone"] <- "Monimolimnion"}
if (samples_df[i,"Fraction"] == "G") { samples_df[i,"Fraction"] <- "Large"}
if (samples_df[i,"Fraction"] == "P") { samples_df[i,"Fraction"] <- "Small"}}
rm(infdataini)

# MetaT table -------------------------------------------------------------
  tableVinput_abundance <- tableVinput
  dt <- setDF(tableVinput_abundance,rownames = tableVinput_abundance$unigene) %>% dplyr::select(-unigene, -PfamDom, -ko, -ko_BestHitNotSignif, -GOterms)
  dt_PA <- dt
  dt_PA[dt_PA > 0 ] <- 1
  MetaT_Summary <- as.data.frame(colSums(dt_PA))
  rm(dt_PA)
  colnames(MetaT_Summary) <- "unigene_raw_count"
  MetaT_Summary$Ref..collab <- row.names(MetaT_Summary)
  MetaT_Summaryx <- merge(MetaT_Summary,samples_df,by="Ref..collab")
  rm(MetaT_Summary)
# Prepare input table -------------------------------------------------------
  
  tableVinput_abundance <- setDF(tableVinput_abundance,rownames = tableVinput_abundance$unigene) %>% dplyr::select(-unigene, -PfamDom, -ko, -ko_BestHitNotSignif, -GOterms)
  #tableVinput_abundance <- tableVinput_abundance[which(rowSums(tableVinput_abundance)>100),]
#
# Stat 2 and first filter Somme = 0 --------------------------------------------------------------------
  tableVinput_abundance$Somme <- rowSums(tableVinput_abundance)
  nrow(tableVinput_abundance)
  nrow(tableVinput_abundance %>% filter(Somme != 0))
  tableVinput_abundance0 <- tableVinput_abundance %>% dplyr::filter(Somme == 0)
  tableVinput_abundance <- tableVinput_abundance %>% dplyr::filter(Somme != 0) %>% dplyr::select(-Somme)
#

# Test SERE on real duplicat ----------------------------------------------
  # prepare table for SERE on duplicat
  tableVinput_abundance_SERE <- merge(tableVinput_abundance,tableVinput %>% dplyr::select(unigene,PfamDom,ko,ko_BestHitNotSignif,GOterms),by.x="row.names",by.y="unigene")
  interest <- tableVinput_abundance_SERE %>% filter(grepl(ko_BestHitNotSignif,pattern=",")) %>% dplyr::select(Row.names)
  tableVinput_abundance_SERE[,"ko_BestHitNotSignif"][tableVinput_abundance_SERE[,"Row.names"] %in% interest$Row.names] <- NA
  tableVinput_abundance_SERE[,"ko_BestHitNotSignif"][is.na(tableVinput_abundance_SERE[,"ko"]) == FALSE] <- NA
  tableVinput_abundance_SERE <- tableVinput_abundance_SERE %>% tidyr::unite(., col = "Kegg.Onthology",  ko, ko_BestHitNotSignif, na.rm=TRUE, sep = ",",remove=FALSE)
  tableVinput_abundance_SERE[tableVinput_abundance_SERE==""] <- NA
  #
  #tableVinput_abundance_SERE_x <- tableVinput_abundance_SERE %>% select(-PfamDom, -ko, -ko_BestHitNotSignif, -GOterms, -Row.names)
  #tableVinput_abundance_SERE_x <- separate_rows(tableVinput_abundance_SERE_x,"Kegg.Onthology",sep = ",") # Duplicate row by KO
  #tableVinput_abundance_SERE_x <- tableVinput_abundance_SERE_x %>% group_by(Kegg.Onthology) %>% summarise_all(sum)
  #
  tableVinput_abundance_SERE_x <- tableVinput_abundance_SERE # %>% filter(is.na(Kegg.Onthology) != TRUE)
  # Calculate SERE
  compare_dup <- MetaT_Summaryx
  compare_dup$Replicat_YON <- "No"
  compare_dup$EchantillonS <- compare_dup$Echantillon
  compare_dup$SERE_index <- 0
  for (i in rownames(compare_dup)) {
    condition_courante <- compare_dup[i,"Condition"]
    periode_courante <- compare_dup[i,"Periods"]
      k <- compare_dup %>% dplyr::filter(Condition == condition_courante) %>% dplyr::filter(Periods == periode_courante)
      if (nrow(k) > 1) { compare_dup[i,"Replicat_YON"] <- "Yes"
      A <- k$Ref..collab[1]
      Ax <- k$Echantillon[1]
      B <- k$Ref..collab[2]
      Bx <- k$Echantillon[2]
      compare_courant <- as.matrix(tableVinput_abundance_SERE_x %>% select(A,B))
      SERE_courant <- SERE(compare_courant)
      #SERE_courant <- (cor.test(compare_courant[,A],compare_courant[,B], alternative="greater"))$estimate
      compare_dup[i,"EchantillonS"] <- paste(Ax,Bx, sep="-")
      compare_dup[i,"SERE_index"] <- SERE_courant
      }
    }
  compare_dup <- compare_dup %>% filter(Replicat_YON == "Yes") %>% select(Condition, Periods, EchantillonS, SERE_index) %>% distinct()
  compare_dup_SERE_summary <- summary(compare_dup$SERE_index)
#
# Remove duplicat ---------------------------------------------------------
    tuniq <- unique(MetaT_Summaryx$Condition)
    puniq <- unique(MetaT_Summaryx$Periods)
    MetaT_Summary <- tibble()
    for (i in tuniq) {
      for (j in puniq) {
        k <- MetaT_Summaryx %>% dplyr::filter(Condition == i) %>% dplyr::filter(Periods == j)
        if (nrow(k) > 1) { interest <- k %>% dplyr::filter(Replicat == 1) }
        else { interest <- k }
        MetaT_Summary <- rbind(MetaT_Summary, interest)
      }
    }
    oldcol <- colnames(tableVinput_abundance)
    w <- 0
    for (i in oldcol) { f <- as.character(MetaT_Summary %>% dplyr::filter(Ref..collab == i) %>% dplyr::select(Echantillon))
      if (f != "character(0)") { 
        if (w == 0) { newcol <- f 
          selectcol <- i }
        else { newcol <- c(all_of(newcol),f)
          selectcol <- c(all_of(selectcol),i)}}
    w <- w+1 }
    tableVinput_abundance <- tableVinput_abundance %>% dplyr::select(all_of(selectcol))
    colnames(tableVinput_abundance) <- newcol
#
# Stat 3 and second filter Somme = 0 --------------------------------------------------------------------
    tableVinput_abundance$Somme <- rowSums(tableVinput_abundance)
    nrow(tableVinput_abundance)
    nrow(tableVinput_abundance %>% filter(Somme != 0))
    tableVinput_abundance0 <- tableVinput_abundance %>% dplyr::filter(Somme == 0)
    tableVinput_abundance <- tableVinput_abundance %>% dplyr::filter(Somme != 0) %>% dplyr::select(-Somme)
#
# Create sorting condition ------------------------------------------------
    #Cycle
    CycleDay <- MetaT_Summary %>% dplyr::filter(Cycle == "Day")
    CycleDay <- CycleDay$Echantillon
    CycleNight <- MetaT_Summary %>% dplyr::filter(Cycle == "Night")
    CycleNight <- CycleNight$Echantillon
    #Zone
    ZoneMixolimnion <- MetaT_Summary %>% dplyr::filter(`Zone` == "Mixolimnion")
    ZoneMixolimnion <- ZoneMixolimnion$Echantillon
    ZoneMonimolimnion <- MetaT_Summary %>% dplyr::filter(`Zone` == "Monimolimnion")
    ZoneMonimolimnion <- ZoneMonimolimnion$Echantillon
    #Fraction
    FractionSmall <- MetaT_Summary %>% dplyr::filter(`Fraction` == "Small")
    FractionSmall <- FractionSmall$Echantillon
    FractionLarge <- MetaT_Summary %>% dplyr::filter(`Fraction` == "Large")
    FractionLarge <- FractionLarge$Echantillon
    #FractionXZone
    SmallXMonimolimnion <- MetaT_Summary %>% dplyr::filter(`Fraction` == "Small") %>% dplyr::filter(`Zone` == "Monimolimnion")
    SmallXMonimolimnion <- SmallXMonimolimnion$Echantillon
    SmallXMixolimnion <- MetaT_Summary %>% dplyr::filter(`Fraction` == "Small") %>% dplyr::filter(`Zone` == "Mixolimnion")
    SmallXMixolimnion <- SmallXMixolimnion$Echantillon
    LargeXMonimolimnion <- MetaT_Summary %>% dplyr::filter(`Fraction` == "Large") %>% dplyr::filter(`Zone` == "Monimolimnion")
    LargeXMonimolimnion <- LargeXMonimolimnion$Echantillon
    LargeXMixolimnion <- MetaT_Summary %>% dplyr::filter(`Fraction` == "Large") %>% dplyr::filter(`Zone` == "Mixolimnion")
    LargeXMixolimnion <- LargeXMixolimnion$Echantillon
    #FractionXZoneXCycle
    SmallXMonimolimnionXDay <- MetaT_Summary %>% dplyr::filter(`Fraction` == "Small") %>% dplyr::filter(`Zone` == "Monimolimnion") %>% dplyr::filter(`Cycle` == "Day")
    SmallXMonimolimnionXDay <- SmallXMonimolimnionXDay$Echantillon
    SmallXMixolimnionXDay <- MetaT_Summary %>% dplyr::filter(`Fraction` == "Small") %>% dplyr::filter(`Zone` == "Mixolimnion") %>% dplyr::filter(`Cycle` == "Day")
    SmallXMixolimnionXDay <- SmallXMixolimnionXDay$Echantillon
    LargeXMonimolimnionXDay <- MetaT_Summary %>% dplyr::filter(`Fraction` == "Large") %>% dplyr::filter(`Zone` == "Monimolimnion") %>% dplyr::filter(`Cycle` == "Day")
    LargeXMonimolimnionXDay <- LargeXMonimolimnionXDay$Echantillon
    LargeXMixolimnionXDay <- MetaT_Summary %>% dplyr::filter(`Fraction` == "Large") %>% dplyr::filter(`Zone` == "Mixolimnion") %>% dplyr::filter(`Cycle` == "Day")
    LargeXMixolimnionXDay <- LargeXMixolimnionXDay$Echantillon
    
    SmallXMonimolimnionXNight <- MetaT_Summary %>% dplyr::filter(`Fraction` == "Small") %>% dplyr::filter(`Zone` == "Monimolimnion") %>% dplyr::filter(`Cycle` == "Night")
    SmallXMonimolimnionXNight <- SmallXMonimolimnionXNight$Echantillon
    SmallXMixolimnionXNight <- MetaT_Summary %>% dplyr::filter(`Fraction` == "Small") %>% dplyr::filter(`Zone` == "Mixolimnion") %>% dplyr::filter(`Cycle` == "Night")
    SmallXMixolimnionXNight <- SmallXMixolimnionXNight$Echantillon
    LargeXMonimolimnionXNight <- MetaT_Summary %>% dplyr::filter(`Fraction` == "Large") %>% dplyr::filter(`Zone` == "Monimolimnion") %>% dplyr::filter(`Cycle` == "Night")
    LargeXMonimolimnionXNight <- LargeXMonimolimnionXNight$Echantillon
    LargeXMixolimnionXNight <- MetaT_Summary %>% dplyr::filter(`Fraction` == "Large") %>% dplyr::filter(`Zone` == "Mixolimnion") %>% dplyr::filter(`Cycle` == "Night")
    LargeXMixolimnionXNight <- LargeXMixolimnionXNight$Echantillon
    
    #Periods
    Periods04 <- MetaT_Summary %>% dplyr::filter(Periods == "04")
    Periods04 <- Periods04$Echantillon
    Periods06 <- MetaT_Summary %>% dplyr::filter(Periods == "06")
    Periods06 <- Periods06$Echantillon
    Periods09 <- MetaT_Summary %>% dplyr::filter(Periods == "09")
    Periods09 <- Periods09$Echantillon
    Periods11 <- MetaT_Summary %>% dplyr::filter(Periods == "11")
    Periods11 <- Periods11$Echantillon
    #
    
    
    
# Test SERE on CDT --------------------------------------------------------------------
  # prepare table for SERE on duplicat
    tableVinput_abundance_SERE <- merge(tableVinput_abundance,tableVinput %>% dplyr::select(unigene,PfamDom,ko,ko_BestHitNotSignif,GOterms),by.x="row.names",by.y="unigene")
    interest <- tableVinput_abundance_SERE %>% filter(grepl(ko_BestHitNotSignif,pattern=",")) %>% dplyr::select(Row.names)
    tableVinput_abundance_SERE[,"ko_BestHitNotSignif"][tableVinput_abundance_SERE[,"Row.names"] %in% interest$Row.names] <- NA
    tableVinput_abundance_SERE[,"ko_BestHitNotSignif"][is.na(tableVinput_abundance_SERE[,"ko"]) == FALSE] <- NA
    tableVinput_abundance_SERE <- tableVinput_abundance_SERE %>% tidyr::unite(., col = "Kegg.Onthology",  ko, ko_BestHitNotSignif, na.rm=TRUE, sep = ",",remove=FALSE)
    tableVinput_abundance_SERE[tableVinput_abundance_SERE==""] <- NA
    #
    #tableVinput_abundance_SERE_x <- tableVinput_abundance_SERE %>% select(-PfamDom, -ko, -ko_BestHitNotSignif, -GOterms, -Row.names)
    #tableVinput_abundance_SERE_x <- separate_rows(tableVinput_abundance_SERE_x,"Kegg.Onthology",sep = ",") # Duplicate row by KO
    #tableVinput_abundance_SERE_x <- tableVinput_abundance_SERE_x %>% group_by(Kegg.Onthology) %>% summarise_all(sum)
    #
    tableVinput_abundance_SERE_x <- tableVinput_abundance_SERE # %>% filter(is.na(Kegg.Onthology) != TRUE)
    # Calculate SERE
# Calculate SERE
  #Cycle
    compare_Cycle <- MetaT_Summary
    compare_Cycle$EchantillonS <- compare_Cycle$Ref..collab
    compare_Cycle$SERE_index <- 0
    for (i in rownames(compare_Cycle)) {
      zone_courante <- compare_Cycle[i,"Zone"]
      fraction_courante <- compare_Cycle[i,"Fraction"]
      periode_courante <- compare_Cycle[i,"Periods"]
      k <- compare_Cycle %>% dplyr::filter(Zone == zone_courante) %>% dplyr::filter(Fraction == fraction_courante) %>% dplyr::filter(Periods == periode_courante)
      A <- k$Echantillon[1]
      B <- k$Echantillon[2]
      compare_courant <- as.matrix(tableVinput_abundance_SERE_x %>% select(A,B))
      SERE_courant <- SERE(compare_courant)
      compare_Cycle[i,"EchantillonS"] <- paste(A,B, sep="+")
      compare_Cycle[i,"SERE_index"] <- SERE_courant
    }
    compare_Cycle <- compare_Cycle %>% select(Zone, Fraction, Periods, EchantillonS, SERE_index) %>% distinct()
    compare_Cycle_SERE_summary <- summary(compare_Cycle$SERE_index)
  #Zone
    compare_Zone <- MetaT_Summary
    compare_Zone$EchantillonS <- compare_Zone$Ref..collab
    compare_Zone$SERE_index <- 0
    for (i in rownames(compare_Zone)) {
      cycle_courante <- compare_Zone[i,"Cycle"]
      fraction_courante <- compare_Zone[i,"Fraction"]
      periode_courante <- compare_Zone[i,"Periods"]
      k <- compare_Zone %>% dplyr::filter(Cycle == cycle_courante) %>% dplyr::filter(Fraction == fraction_courante) %>% dplyr::filter(Periods == periode_courante)
      A <- k$Echantillon[1]
      B <- k$Echantillon[2]
      compare_courant <- as.matrix(tableVinput_abundance_SERE_x %>% select(A,B))
      SERE_courant <- SERE(compare_courant)
      compare_Zone[i,"EchantillonS"] <- paste(A,B, sep="+")
      compare_Zone[i,"SERE_index"] <- SERE_courant
    }
    compare_Zone <- compare_Zone %>% select(Cycle, Fraction, Periods, EchantillonS, SERE_index) %>% distinct()
    compare_Zone_SERE_summary <- summary(compare_Zone$SERE_index)
  #Fraction
    compare_Fraction <- MetaT_Summary
    compare_Fraction$EchantillonS <- compare_Fraction$Ref..collab
    compare_Fraction$SERE_index <- 0
    for (i in rownames(compare_Fraction)) {
      cycle_courante <- compare_Fraction[i,"Cycle"]
      zone_courante <- compare_Fraction[i,"Zone"]
      periode_courante <- compare_Fraction[i,"Periods"]
      k <- compare_Fraction %>% dplyr::filter(Cycle == cycle_courante) %>% dplyr::filter(Zone == zone_courante) %>% dplyr::filter(Periods == periode_courante)
      A <- k$Echantillon[1]
      B <- k$Echantillon[2]
      compare_courant <- as.matrix(tableVinput_abundance_SERE_x %>% select(A,B))
      SERE_courant <- SERE(compare_courant)
      compare_Fraction[i,"EchantillonS"] <- paste(A,B, sep="+")
      compare_Fraction[i,"SERE_index"] <- SERE_courant
    }
    compare_Fraction <- compare_Fraction %>% select(Cycle, Zone, Periods, EchantillonS, SERE_index) %>% distinct()
    compare_Fraction_SERE_summary <- summary(compare_Fraction$SERE_index)
  # All combination
    compare_all <- as.data.frame(t(combn(MetaT_Summary$Echantillon,2)))
    compare_all$SERE_index <- 0
    compare_all$Comparaison_type <- "Mixed"
    compare_all$Comparaison <- "Mixed"
    for (i in rownames(compare_all)) {
      A <- compare_all[i,"V1"]
      B <- compare_all[i,"V2"]
      compare_courant <- as.matrix(tableVinput_abundance_SERE_x %>% select(A,B))
      SERE_courant <- SERE(compare_courant)
      #SERE_courant <- (cor.test(compare_courant[,A],compare_courant[,B], alternative="greater"))$estimate
      compare_all[i,"SERE_index"] <- SERE_courant
      ZoneA <- MetaT_Summary[,"Zone"][MetaT_Summary[,"Echantillon"] == A]
      ZoneB <- MetaT_Summary[,"Zone"][MetaT_Summary[,"Echantillon"] == B]
      FractionA <- MetaT_Summary[,"Fraction"][MetaT_Summary[,"Echantillon"] == A]
      FractionB <- MetaT_Summary[,"Fraction"][MetaT_Summary[,"Echantillon"] == B]
      CycleA <- MetaT_Summary[,"Cycle"][MetaT_Summary[,"Echantillon"] == A]
      CycleB <- MetaT_Summary[,"Cycle"][MetaT_Summary[,"Echantillon"] == B]
      PeriodsA <- MetaT_Summary[,"Periods"][MetaT_Summary[,"Echantillon"] == A]
      PeriodsB <- MetaT_Summary[,"Periods"][MetaT_Summary[,"Echantillon"] == B]
      if (ZoneA == ZoneB & FractionA == FractionB & CycleA == CycleB & PeriodsA != PeriodsB) {
        compare_all[i,"Comparaison_type"] <- "Periods"
        compare_all[i,"Comparaison"] <- paste(PeriodsA,"vs",PeriodsB,"with",CycleA,"&",FractionA,"&",ZoneA)}
      if (ZoneA == ZoneB & FractionA == FractionB & CycleA != CycleB & PeriodsA == PeriodsB) {
        compare_all[i,"Comparaison_type"] <- "Cycle"
        compare_all[i,"Comparaison"] <- paste(CycleA,"vs",CycleB,"with",PeriodsA,"&",FractionA,"&",ZoneA)}
      if (ZoneA == ZoneB & FractionA != FractionB & CycleA == CycleB & PeriodsA == PeriodsB) {
        compare_all[i,"Comparaison_type"] <- "Fraction"
        compare_all[i,"Comparaison"] <- paste(FractionA,"vs",FractionB,"with",PeriodsA,"&",CycleA,"&",ZoneA)}
      if (ZoneA != ZoneB & FractionA == FractionB & CycleA == CycleB & PeriodsA == PeriodsB) {
        compare_all[i,"Comparaison_type"] <- "Zone"
        compare_all[i,"Comparaison"] <- paste(ZoneA,"vs",ZoneB,"with",PeriodsA,"&",CycleA,"&",ZoneA)}
    }
    
  # Plot SERE ---------------------------------------------------------------
    #Prepare mix of table (duplicat and condition)
    compare_mix <- compare_all %>% select(-Comparaison)
    compare_dup_temp <- separate(compare_dup %>% select(-Periods, -Condition),EchantillonS,c("V1","V2"), sep = "-")
    compare_dup_temp$Comparaison_type <- "Biological replicat"
    compare_mix <- rbind(compare_mix,compare_dup_temp)
    #Plot
    my_comp <- as.list(as.data.frame(combn(unique(compare_mix$Comparaison_type),2)))
    palet_periodcomp <- c("#7375b5FF","#8ca252FF","#f6ae34FF","#ad494aFF","#abababFF","#abfeabFF")
    aidiv <- ggplot(compare_mix, aes(x = Comparaison_type, y = SERE_index)) + geom_violin(aes(fill = Comparaison_type),trim = FALSE,width=0.5,alpha=0.9) + geom_boxplot(width=0.075) + geom_jitter(position=position_jitter(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      scale_fill_manual(values = palet_periodcomp) + theme_unique_art() +
      scale_y_continuous(breaks =c(-1,0,1,2,3,4,5,6)) +
      labs(x="Condition comparaison",y="SERE index") + guides(fill = "none")
    print(aidiv)
    ggsave("Stat/SERE_on_all_unigene.svg", device = "svg", width = 7, height = 8)
#
# ANOSIM & ADONIS ------------------------------------------------------------------
    Unigene_table_c <- tableVinput_abundance
    dt <- t(Unigene_table_c)
    dt <- merge(dt, MetaT_Summary %>% select(Echantillon, Cycle, Periods, Fraction, Zone), by.x = "row.names", by.y = "Echantillon")
    dtx <- dt %>% select(-Cycle,-Periods, -Fraction, -Zone, -Row.names)
    dist <- vegdist(dtx, method ="bray")
#
    # Adonis
    capture.output(adonis2(formula = dist~dt$Zone+dt$Fraction+dt$Cycle+dt$Periods,distance = "bray", parallel = 10),file="Betadisp/Adonis.txt")
    # Anosim
    capture.output(anosim(dist, dt$Zone, permutations = 9999, parallel = 10),file="Betadisp/anosim_transcrit_Zone.txt")
    capture.output(anosim(dist, dt$Fraction, permutations = 9999, parallel = 10),file="Betadisp/anosim_transcrit_Fraction.txt")
    capture.output(anosim(dist, dt$Periods, permutations = 9999, parallel = 10),file="Betadisp/anosim_transcrit_Periods.txt")
    capture.output(anosim(dist, dt$Cycle, permutations = 9999, parallel = 10),file="Betadisp/anosim_transcrit_Cycle.txt")
#
# DESEQ2 on all transcript ------------------------------------------------------------------
    if (alltranscript_analysis == "yes") {
    # metadata 
    meta <- MetaT_Summary %>% select(-Ref..collab)
    row.names(meta) <- meta$Echantillon
    all(colnames(data) %in% rownames(meta))
    all(colnames(data) == rownames(meta))
    meta$Zone <- as.factor(meta$Zone)
    meta$Fraction <- as.factor(meta$Fraction)
    meta$Periods <- as.factor(meta$Periods)
    meta$Cycle <- as.factor(meta$Cycle)
    # Start normalization
    modele <- ~ Fraction+Zone+Periods
    dds <- DESeqDataSetFromMatrix(countData = tableVinput_abundance, colData = meta, 
                                  design = formula(modele))
    dds <- DESeq(dds)
#
  # Summary DESEQ2 on transcript --------------------------------------------
    ann_colors = list(
      Zone = c(Mixolimnion = "lightgrey", Monimolimnion = "firebrick"),
      Fraction = c(Small = "#1B9E77", Large = "#D95F02"),
      Periods = c(`04` = "#6d6c88", `06` = "#7b8e4b", `09` = "#e3a235", `11` = "#a146a7"),
      Cycle = c(Day = "#a14647", Night = "#6d6ea8"))
    # Summary figure
    vsd <- vst(dds, blind=FALSE)
    assay_vsd <- assay(vsd)
    sampleDists <- vegdist(t(assay_vsd), method="euclidean")
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(vsd$Zone, vsd$Fraction, vsd$Cycle, vsd$Periods, sep="-")
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    map <- pheatmap(sampleDistMatrix, show_rownames=FALSE,
                    clustering_distance_rows=sampleDists,
                    clustering_distance_cols=sampleDists,
                    annotation_row = meta %>% select(Zone,Fraction, Cycle, Periods),
                    annotation_colors = ann_colors,
                    col=colors)
    svglite("DESeq2/Summary/Heatmapsummary_all_transcripts.svg",width = 8.00,height = 6.00)
    print(map)
    dev.off()
    #PCA
    #prepare legend
    #color
    legend_color <- plotPCA(vsd, intgroup=c("Zone","Periods")) +
      geom_point(size=6,aes(color = Zone:Periods)) + 
      scale_color_manual("Zone:Periods",values = c("#720500","#AB0800","#CC0A00","#FF584F",
                                                            "#003E49","#0091AB","#00D6FC","#8FEDFE")) + 
                                                              theme(legend.title = element_text(face = "bold"))
    legend_color <- get_legend(legend_color)
    #shape
    legend_shape <- plotPCA(vsd, intgroup=c("Fraction")) + guides(color ="none") +
      geom_point(size=6,aes(shape = Fraction),color="black") + 
      theme(legend.title = element_text(face = "bold"))
    legend_shape <- get_legend(legend_shape)
    # Grid
    legend <- plot_grid(legend_color,legend_shape,ncol=1,nrow=2,align = "hv",rel_heights = c(1,1))
    # plot real PCA
    pca <- plotPCA(vsd, intgroup=c("Zone","Fraction","Periods","Cycle")) +
      theme_unique_art() + theme(legend.position = "none") +
      geom_polygon(aes(group=Zone:Fraction:Periods),color ="black") +
      geom_polygon(aes(group=Periods:Zone:Cycle),fill="#FFFFFF00",color ="black") +
      geom_point(size=4.5,aes(shape=Fraction,color = Zone:Fraction:Periods:Cycle)) +
      scale_color_manual(values = c("#720500","#720500","#AB0800","#AB0800","#CC0A00","#CC0A00","#FF584F","#FF584F","#720500","#720500","#AB0800","#AB0800","#CC0A00","#CC0A00","#FF584F","#FF584F",
                                             "#003E49","#003E49","#0091AB","#0091AB","#00D6FC","#00D6FC","#8FEDFE","#8FEDFE","#003E49","#003E49","#0091AB","#0091AB","#00D6FC","#00D6FC","#8FEDFE","#8FEDFE"))
    svglite("DESeq2/Summary/PCAsummary_all_transcripts.svg",width = 12.00,height = 6.00)
    plot_grid(pca,NULL,legend,nrow = 1,ncol = 3, rel_widths = c(1.5,0.1,0.4))
    dev.off()
#
    }
# Test SERE with CDT regroupment -----------------------------------
  # prepare tables and calculate SERE ----------------------------------------------------------
  #Cycle
    ## Day
    Day_seq <- as.data.frame(rowSums(tableVinput_abundance %>% dplyr::select(all_of(CycleDay))))
    colnames(Day_seq) <- "TotalDay"
    ## Night
    Night_seq <- as.data.frame(rowSums(tableVinput_abundance %>% dplyr::select(all_of(CycleNight))))
    colnames(Night_seq) <- "TotalNight"
    ## Merge
    Cycle_seq <- merge(x = Day_seq,y = Night_seq, by = "row.names")
    compare <- as.matrix(Cycle_seq %>% select(-Row.names))
    SERE(compare)
  # Zone
    ## Mixolimnion
    Mixolimnion_seq <- as.data.frame(rowSums(tableVinput_abundance %>% dplyr::select(all_of(ZoneMixolimnion))))
    colnames(Mixolimnion_seq) <- "TotalMixolimnion"
    ## Monimolimnion
    Monimolimnion_seq <- as.data.frame(rowSums(tableVinput_abundance %>% dplyr::select(all_of(ZoneMonimolimnion))))
    colnames(Monimolimnion_seq) <- "TotalMonimolimnion"
    ## Merge
    Zone_seq <- merge(x = Mixolimnion_seq,y = Monimolimnion_seq, by = "row.names")
    compare <- as.matrix(Zone_seq %>% select(-Row.names))
    SERE(compare)
  # Fraction
    ## Small
    Small_seq <- as.data.frame(rowSums(tableVinput_abundance %>% dplyr::select(all_of(FractionSmall))))
    colnames(Small_seq) <- "TotalSmall"
    ## Large
    Large_seq <- as.data.frame(rowSums(tableVinput_abundance %>% dplyr::select(all_of(FractionLarge))))
    colnames(Large_seq) <- "TotalLarge"
    ## Merge
    Fraction_seq <- merge(x = Small_seq,y = Large_seq, by = "row.names")
    compare <- as.matrix(Fraction_seq %>% select(-Row.names))
    SERE(compare)
#
# Merge Tables ------------------------------------------------------------
  tableVinput_abundance_annots <- merge(tableVinput_abundance,tableVinput %>% dplyr::select(unigene,PfamDom,ko,ko_BestHitNotSignif,GOterms),by.x="row.names",by.y="unigene")
  interest <- tableVinput_abundance_annots %>% filter(grepl(ko_BestHitNotSignif,pattern=",")) %>% dplyr::select(Row.names)
  tableVinput_abundance_annots[,"ko_BestHitNotSignif"][tableVinput_abundance_annots[,"Row.names"] %in% interest$Row.names] <- NA
  tableVinput_abundance_annots[,"ko_BestHitNotSignif"][is.na(tableVinput_abundance_annots[,"ko"]) == FALSE] <- NA
  tableVinput_abundance_annots <- tableVinput_abundance_annots %>% tidyr::unite(., col = "Kegg.Onthology",  ko, ko_BestHitNotSignif, na.rm=TRUE, sep = ",",remove=FALSE)
  tableVinput_abundance_annots[tableVinput_abundance_annots==""] <- NA
  #
  tableVinput_abundance_annots_tax <- merge(tableVinput_abundance_annots,Tax_table, by.x = "Row.names", by.y = "Unigene")
  rm(tableVinput_abundance_annots)
  #
    
    
# Taxonomic Annotation ---------------------------------------------------------
  # Process Taxonomic data ------------------------------------------------------------
    # Prepare Taxonomy
    uniq_tax <- as.data.frame(unique(tableVinput_abundance_annots_tax$Values)) ; colnames(uniq_tax) <- "uniq"
    uniq_tax$Domain <- NA
    uniq_tax$Supergroup <- NA
    uniq_tax$Kingdom <- NA
    uniq_tax$Division <- NA
    uniq_tax$Class <- NA
    uniq_tax$Order <- NA
    uniq_tax$Family <- NA
    uniq_tax$Genus <- NA
    uniq_tax$Species <- NA
    for (i in row.names(uniq_tax)) {
      uniq_tax[i,"Domain"] <- strsplit(strsplit(as.character(uniq_tax[i,"uniq"]), "d_")[[1]][2],";")[[1]][1]
      uniq_tax[i,"Supergroup"] <- strsplit(strsplit(strsplit(as.character(uniq_tax[i,"uniq"]), "Eukaryota;")[[1]][2],";")[[1]][1],"_")[[1]][2]
      uniq_tax[i,"Kingdom"] <- strsplit(strsplit(as.character(uniq_tax[i,"uniq"]), "k_")[[1]][2],";")[[1]][1]
      uniq_tax[i,"Division"] <- strsplit(strsplit(as.character(uniq_tax[i,"uniq"]), "p_")[[1]][2],";")[[1]][1]
      uniq_tax[i,"Class"] <- strsplit(strsplit(as.character(uniq_tax[i,"uniq"]), "c_")[[1]][2],";")[[1]][1]
      uniq_tax[i,"Order"] <- strsplit(strsplit(as.character(uniq_tax[i,"uniq"]), "o_")[[1]][2],";")[[1]][1]
      uniq_tax[i,"Family"] <- strsplit(strsplit(as.character(uniq_tax[i,"uniq"]), "f_")[[1]][2],";")[[1]][1]
      uniq_tax[i,"Genus"] <- strsplit(strsplit(as.character(uniq_tax[i,"uniq"]), "g_")[[1]][2],";")[[1]][1]
      uniq_tax[i,"Species"] <- strsplit(strsplit(as.character(uniq_tax[i,"uniq"]), "s_")[[1]][2],";")[[1]][1]
    }
    # Viridiplantae
    uniq_tax[,"Supergroup"][uniq_tax[,"Supergroup"]=="Viridiplantae"] <- "Archaeplastida"
    uniq_tax[,"Kingdom"][uniq_tax[,"Kingdom"]=="Viridiplantae"] <- "Archaeplastida"
    
    # Haptista
    uniq_tax[,"Supergroup"][uniq_tax[,"Supergroup"]=="Haptista"] <- "Hacrobia"
    # Heterolobosea
    uniq_tax[,"Class"][uniq_tax[,"Supergroup"]=="Heterolobosea"] <- "Heterolobosea"
    uniq_tax[,"Division"][uniq_tax[,"Supergroup"]=="Heterolobosea"] <- "Discoba"
    uniq_tax[,"Supergroup"][uniq_tax[,"Supergroup"]=="Heterolobosea"] <- "Excavata"
    # Katablepharidophyta
    uniq_tax[,"Division"][uniq_tax[,"Supergroup"]=="Katablepharidophyta"] <- "Katablepharidophyta"
    uniq_tax[,"Supergroup"][uniq_tax[,"Supergroup"]=="Katablepharidophyta"] <- "Hacrobia"
    # Cryptophyta
    uniq_tax[,"Division"][uniq_tax[,"Supergroup"]=="Cryptophyta"] <- "Cryptophyta"
    uniq_tax[,"Supergroup"][uniq_tax[,"Supergroup"]=="Cryptophyta"] <- "Hacrobia"
    # Malawimonadidae
    uniq_tax[,"Division"][uniq_tax[,"Supergroup"]=="Malawimonadidae"] <- "Malawimonadidae"
    uniq_tax[,"Supergroup"][uniq_tax[,"Supergroup"]=="Malawimonadidae"] <- "Excavata"
    # Rhodophyta
    uniq_tax[,"Division"][uniq_tax[,"Supergroup"]=="Rhodophyta"] <- "Rhodophyta"
    uniq_tax[,"Supergroup"][uniq_tax[,"Supergroup"]=="Rhodophyta"] <- "Archaeplastida"
    # Glaucocystophyceae
    uniq_tax[,"Class"][uniq_tax[,"Supergroup"]=="Glaucocystophyceae"] <- "Glaucocystophyceae"
    uniq_tax[,"Division"][uniq_tax[,"Supergroup"]=="Glaucocystophyceae"] <- "Glaucophyta"
    uniq_tax[,"Supergroup"][uniq_tax[,"Supergroup"]=="Glaucocystophyceae"] <- "Archaeplastida"
    # Jakobida
    uniq_tax[,"Division"][uniq_tax[,"Supergroup"]=="Jakobida"] <- "Discoba"
    uniq_tax[,"Supergroup"][uniq_tax[,"Supergroup"]=="Jakobida"] <- "Excavata"
    # Euglenozoa
    uniq_tax[,"Division"][uniq_tax[,"Supergroup"]=="Euglenozoa"] <- "Discoba"
    uniq_tax[,"Supergroup"][uniq_tax[,"Supergroup"]=="Euglenozoa"] <- "Excavata"
    # Metamonada
    uniq_tax[,"Division"][uniq_tax[,"Supergroup"]=="Metamonada"] <- "Metamonada"
    uniq_tax[,"Supergroup"][uniq_tax[,"Supergroup"]=="Metamonada"] <- "Excavata"
    # Oomycetes
    uniq_tax[,"Class"][uniq_tax[,"Class"]=="Oomycetes"] <- "Oomycota"
    # Discosea
    uniq_tax[,"Division"][(uniq_tax[,"Division"]=="Discosea")&(uniq_tax[,"Class"]=="Flabellinia")] <- "Discosea-Flabellinia"
    uniq_tax[,"Division"][(uniq_tax[,"Division"]=="Discosea")&(uniq_tax[,"Order"]=="Longamoebia")] <- "Discosea-Longamoebia"
    # Curration
    exept <- c("environmental samples", "environvironmental samples", "unclassified", "unclassified","_X","incertae sedis")
    uniq_tax[,"Last"] <- NA
    uniq_tax[,"level"] <- NA
    # Création de la colonne Last
    uniq_tax[,"Species"][is.na(uniq_tax[,"Genus"])==T] <- NA
    for (i in row.names(uniq_tax)) {
      for (j in c("Domain","Supergroup", "Kingdom","Division","Class","Order","Family","Genus","Species")) {
        if (is.na(uniq_tax[i,j])==F) { 
          if (j == "Species") {uniq_tax[i,"Last"] <- uniq_tax[i,"Genus"]}
          else {uniq_tax[i,"Last"] <- uniq_tax[i,j]}
          uniq_tax[i,"level"] <- j}
      }
    }
    
    # Re-annotation par recherche du dernier levels dans la taxonomie PR2
    PR2_tax <- as.data.frame(PR2_tax)
    PR2_tax$Supergroup <- ""
    PR2_tax[,"Supergroup"] <- sapply(strsplit(as.character(PR2_tax$Taxonomy), ";"), `[`, 2)
    #
    uniq_tax[,"PR2_taxonomy"] <- NA
    #uniq_tax <- uniq_tax %>% dplyr::filter(is.na(level) == F)
    for (i in row.names(uniq_tax)) {
      if (is.na(uniq_tax[i,"level"]) == F) {
        w <- paste0(uniq_tax[i,"Last"],";")
        if (uniq_tax[i,"level"] == "Species") {
          z <- str_remove(str_remove(paste0(gsub(" ","_",uniq_tax[i,"Species"],),";"),"\\["),"\\]")
          grepping <- strsplit(grep(pattern = z,PR2_tax$Taxonomy,value = T)[1],w)[[1]][1]}
        if (uniq_tax[i,"level"] != "Species") {
          grepping <- strsplit(grep(pattern = w,PR2_tax$Taxonomy,value = T)[1],w)[[1]][1]}
        if (is.na(grepping) == F) {
          uniq_tax[i,"PR2_taxonomy"] <- paste0(grepping,w)}
        else { uniq_tax[i,"PR2_taxonomy"] <- NA }
      }}
    # Création de la special table des annotation non re-annotées
    special_table <- uniq_tax %>% dplyr::filter(is.na(PR2_taxonomy) == T)
    special_table <- special_table %>% dplyr::filter(is.na(Last) == F)
    # df1 model table for levels in Taxonomy
    df1=data.frame(number = c(2,3,4,5,6,7,8,9,10),
                   string = c("Domain", "Supergroup", "Kingdom","Division", "Class" , "Order" , "Family", "Genus","Species"))
    # Special table (2ème et xème tours si pas re-annoté)
    for (i in row.names(special_table)) {
      k <- 0
      while (is.na(special_table[i,"PR2_taxonomy"]) == T) {
        j <- df1[,"number"][df1[,"string"] == special_table[i,"level"]]
        j <- as.numeric(j[1])-k
        w <- paste0(special_table[i,j],";")
        k <- k+1
        grepping <- strsplit(grep(pattern = w,PR2_tax$Taxonomy,value = T)[1],w)[[1]][1]
        if (is.na(grepping) == F) {
          special_table[i,"PR2_taxonomy"] <- paste0(grepping,w)}
        else { special_table[i,"PR2_taxonomy"] <- NA }
      }}
    # Replace in uniq_tax
    for (i in row.names(special_table)) {
      uniq_tax[,"PR2_taxonomy"][uniq_tax[,"uniq"]==special_table[i,"uniq"]] <- special_table[i,"PR2_taxonomy"]
    }
    # Verif
    uniq_tax_verif <- separate(uniq_tax, PR2_taxonomy, c("DomainPR2","SupergroupPR2","DivisionPR2","ClassPR2","OrderPR2","FamilyPR2","GenusPR2"), sep =";")
    Domain_t_verif <- uniq_tax_verif %>% dplyr::select(Supergroup,SupergroupPR2)
    verif_table <- uniq_tax_verif %>% dplyr::filter(Supergroup != SupergroupPR2)
    # ALAMANO
    uniq_tax[,"PR2_taxonomy"][uniq_tax[,"uniq"]=="species::Cepedea sp. Rr5::-_cellular organisms;d_Eukaryota;-_Stramenopiles;c_Bigyra;-_Opalozoa;o_Opalinata;f_Opalinidae;g_Cepedea;-_unclassified Cepedea;s_Cepedea sp. Rr5"] <- "Eukaryota;Stramenopiles;Opalozoa;Opalinata;"
    uniq_tax[,"PR2_taxonomy"][uniq_tax[,"uniq"]=="genus::Uronema::-_cellular organisms;d_Eukaryota;k_Viridiplantae;p_Chlorophyta;-_core chlorophytes;c_Chlorophyceae;-_OCC clade;o_Chaetophorales;f_Uronemataceae;g_Uronema"] <- "Eukaryota;Archaeplastida;Chlorophyta;Chlorophyceae;Chaetophorales;"
    # Remove last ";"
    uniq_tax[,"PR2_taxonomy"] <- sapply(substring(as.character(uniq_tax$PR2_taxonomy), 1, nchar(uniq_tax$PR2_taxonomy)-1), `[`, 1)
    # Add Species
    for (i in row.names(uniq_tax)) { if ((is.na(uniq_tax[i,"Species"])==F)&&(str_count(uniq_tax[i,"PR2_taxonomy"], pattern =  ";" )==6)) {uniq_tax[i,"PR2_taxonomy_complete"] <- paste0(uniq_tax[i,"PR2_taxonomy"],";",uniq_tax[i,"Species"])}
      else {uniq_tax[i,"PR2_taxonomy_complete"] <- substring(as.character(uniq_tax[i,"PR2_taxonomy"]),1,nchar(as.character(uniq_tax[i,"PR2_taxonomy"]))-0)}}
    # Merge in tableVinput_abundance_annots_tax
    tableVinput_abundance_annots_tax <- merge(tableVinput_abundance_annots_tax, uniq_tax %>% dplyr::select("uniq", "PR2_taxonomy","PR2_taxonomy_complete"), by.x = "Values", by.y = "uniq")    
#
  # General Tax Stat --------------------------------------------------------------------
    data_seq_summarize <- tableVinput_abundance_annots_tax %>% dplyr::select(Values,Row.names,PR2_taxonomy_complete)
    data_seq_summarize <- separate(data_seq_summarize, PR2_taxonomy_complete, c("Domain","Supergroup","Division","Class","Order","Family","Genus","Species"), sep =";")
    #
    data_seq_summarize[data_seq_summarize == "NA::NA::NA"] <- NA
    data_seq_summarize[data_seq_summarize == "NA::NA::-_cellular organisms;"] <- NA
    data_seq_summarize[data_seq_summarize == "NA::NA::-_root;"] <- NA
    data_seq_summarize[data_seq_summarize == "no rank::root::-_root"] <- NA
    data_seq_summarize[data_seq_summarize == "no rank::cellular organisms::-_cellular organisms"] <- NA
    data_seq_summarize[data_seq_summarize == "no rank::unclassified::NA"] <- NA
    # Taxonomy
    tax_table_TF <- as.data.frame(is.na(data_seq_summarize))
    for (level in c("Domain","Supergroup","Division","Class","Order","Family","Genus","Species")) {
      assign(paste("count", level, sep = "_"), nrow(as.data.frame(tax_table_TF[,level][tax_table_TF[,level] == FALSE])))}
    tax_stat <- as.data.frame(c("Domain","Supergroup","Division","Class","Order","Family","Genus","Species")) ; colnames(tax_stat) <- "Level"
    tax_stat$Affiliated <- c(count_Domain,count_Supergroup,count_Division,count_Class,count_Order,count_Family,count_Genus,count_Species)
    tax_stat$`Not affiliated` <- nrow(data_seq_summarize) - tax_stat$Affiliated
    tax_stat <- reshape2::melt(tax_stat, id = "Level")
    ## Plot
    tax_stat$Level <- factor(tax_stat$Level , levels=c("Domain","Supergroup","Division","Class","Order","Family","Genus","Species"))
    tax_stat$variable <- factor(tax_stat$variable , levels=c("Not affiliated","Affiliated"))
    tax_plot <- ggplot(tax_stat, mapping = aes(x= Level, y = value, fill = variable, color = variable ,linetype = variable), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
      scale_fill_manual(values = alpha(c("#1212ff","#1212ff"), c(0.25,0.5))) +
      geom_label(aes(y = value + 0.08*max(value),label = paste(value,"\n","unigenes: ", round(value*100/nrow(data_seq_summarize),1)," %",sep =""), color = variable, alpha = variable),size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
      scale_color_manual(values = c("#0000FF00","black")) +
      scale_alpha_manual(values = c(0,0.5)) +
      scale_linetype_manual(values = c("blank","dashed")) +
      guides(fill=guide_legend("Affiliation"), color = "none", linetype = "none") +
      labs(x="Taxonomic level",y="Unigenes") +
      theme(legend.title = element_text(face="bold"),
            axis.title = element_text(color = "black", face = "bold"))
    print(tax_plot)
    ggsave("Stat/Stat-Taxo.svg", device = "svg", width = 12, height = 4)
    # Save IMG ----------------------------------------------------------------
    if (inputmode == FALSE) {save.image("R_Taxonomy_no_normalizedOK_image.RData")}
#
# Prepare ko table identification -----------------------------------------
    ko2Met_Table <- as.data.frame(ko2Met_Table)
    ko2Met_Table_Interest <- ko2Met_Table %>% filter(lvl_A_val %in% c("Metabolism","Cellular_Processes"))
    #chitinaseko <- c("K01183","K20547","K13381")

    # Metabolisme #2
    ko2Met_Table_Final <- ko2Met_Table_Interest %>% filter(lvl_A_val == "Metabolism") %>% mutate('Pathway/Gene' = lvl_C_val)
    ko2Met_Table_Final <- ko2Met_Table_Final %>% filter(lvl_B_val != "Not_included_in_regular_maps")
    ko2Met_Table_Final <- ko2Met_Table_Final %>% filter(lvl_B_val != "Chemical_structure_transformation_maps")
    ko2Met_Table_Final <- ko2Met_Table_Final %>% filter(lvl_C_val != "Carbon_fixation_pathways_in_prokaryotes")
    ko2Met_Table_Final[,'Pathway/Gene'][ko2Met_Table_Final[,'lvl_B_val']=="Carbohydrate_metabolism"] <- "Other Carbohydrate metabolism"
    ko2Met_Table_Final[,'Pathway/Gene'][ko2Met_Table_Final[,'lvl_C_val']=="Glycolysis_/_Gluconeogenesis"] <- "Glycolysis and Gluconeogenesis"
    ko2Met_Table_Final[,'Pathway/Gene'][ko2Met_Table_Final[,'lvl_C_val']=="Citrate_cycle_(TCA_cycle)"] <- "TCA cycle"
    ko2Met_Table_Final[,'Pathway/Gene'][ko2Met_Table_Final[,'lvl_C_val']=="Photosynthesis_-_antenna_proteins"] <- "Photosynthesis antenna proteins"
    
    ko2Met_Table_Final[,'Pathway/Gene'][ko2Met_Table_Final[,'lvl_B_val']=="Lipid_metabolism"] <- "Lipid metabolism"

    ko2Met_Table_Final[,'Pathway/Gene'][ko2Met_Table_Final[,'lvl_B_val']=="Nucleotide_metabolism"] <- "Nucleotide metabolism"
    ko2Met_Table_Final[,'Pathway/Gene'][ko2Met_Table_Final[,'lvl_B_val']=="Amino_acid_metabolism"] <- "Amino acid metabolism"
    ko2Met_Table_Final[,'Pathway/Gene'][ko2Met_Table_Final[,'lvl_B_val']=="Metabolism_of_other_amino_acids"] <- "Amino acid metabolism"
    ko2Met_Table_Final[,'Pathway/Gene'][ko2Met_Table_Final[,'lvl_B_val']=="Glycan_biosynthesis_and_metabolism"] <- "Glycan metabolism"
    ko2Met_Table_Final[,'Pathway/Gene'][ko2Met_Table_Final[,'lvl_B_val']=="Metabolism_of_cofactors_and_vitamins"] <- "Cofactors and vitamins metabolism"
    ko2Met_Table_Final[,'Pathway/Gene'][ko2Met_Table_Final[,'lvl_B_val']=="Metabolism_of_terpenoids_and_polyketides"] <- "Terpenoids and polyketides metabolism"
    ko2Met_Table_Final[,'Pathway/Gene'][ko2Met_Table_Final[,'lvl_B_val']=="Biosynthesis_of_other_secondary_metabolites"] <- "Biosynthesis of other secondary metabolites"
    ko2Met_Table_Final[,'Pathway/Gene'][ko2Met_Table_Final[,'lvl_B_val']=="Xenobiotics_biodegradation_and_metabolism"] <- "Xenobiotics biodegradation and metabolism"
    ko2Met_Table_Final[,'lvl_B_val'][ko2Met_Table_Final[,'Pathway/Gene']=="Calvin Cycle"] <- "Energy_metabolism"
    
    ## cellular process
    #ko2Met_Table_Final <- ko2Met_Table_Interest %>% filter(lvl_A_val == "Cellular_Processes") %>% mutate('Pathway/Gene' = lvl_C_val)
    ko2Met_Table_Final <- rbind(ko2Met_Table_Final,ko2Met_Table_Interest %>% filter(lvl_A_val == "Cellular_Processes") %>% mutate('Pathway/Gene' = lvl_C_val))
    ko2Met_Table_Final <- ko2Met_Table_Final %>% filter(lvl_B_val != "Cellular_community_-_prokaryotes")
    ko2Met_Table_Final <- ko2Met_Table_Final %>% filter(lvl_C_val != "Autophagy_-_animal")
    ko2Met_Table_Final <- ko2Met_Table_Final %>% filter(lvl_C_val != "Autophagy_-_yeast")
    ko2Met_Table_Final <- ko2Met_Table_Final %>% filter(lvl_C_val != "Mitophagy_-_animal")
    ko2Met_Table_Final <- ko2Met_Table_Final %>% filter(lvl_C_val != "Mitophagy_-_yeast")
    ko2Met_Table_Final <- ko2Met_Table_Final %>% filter(lvl_C_val != "Cell_cycle_-_yeast")
    ko2Met_Table_Final <- ko2Met_Table_Final %>% filter(lvl_C_val != "Cell_cycle_-_Caulobacter")
    ko2Met_Table_Final <- ko2Met_Table_Final %>% filter(lvl_C_val != "Meiosis_-_yeast")
    ko2Met_Table_Final <- ko2Met_Table_Final %>% filter(lvl_C_val != "Apoptosis_-_fly")
    ko2Met_Table_Final <- ko2Met_Table_Final %>% filter(lvl_C_val != "Bacterial_chemotaxis")
    ko2Met_Table_Final <- ko2Met_Table_Final %>% filter(lvl_C_val != "Oocyte_meiosis")
    ko2Met_Table_Final <- ko2Met_Table_Final %>% filter(lvl_C_val != "Signaling_pathways_regulating_pluripotency_of_stem_cells")
    ko2Met_Table_Final <- ko2Met_Table_Final %>% filter(lvl_C_val != "p53_signaling_pathway")
    ko2Met_Table_Final <- ko2Met_Table_Final %>% filter(lvl_C_val != "Carbon_fixation_in_photosynthetic_organisms")
    
    
    
    ko2Met_Table_Final[,"lvl_C_val"] <- sapply(strsplit(as.character(ko2Met_Table_Final$lvl_C_val),"_-_"), `[`, 1)
    ko2Met_Table_Final[,"lvl_C_val"] <- sapply(str_replace_all(as.character(ko2Met_Table_Final$lvl_C_val),"_"," "), `[`, 1)
    
    ko2Met_Table_Final[,'Pathway/Gene'][ko2Met_Table_Final[,'lvl_C_val']=="Gap junction"] <- "Cell-Cell junction"
    ko2Met_Table_Final[,'Pathway/Gene'][ko2Met_Table_Final[,'lvl_C_val']=="Tight junction"] <- "Cell-Cell junction"
    ko2Met_Table_Final[,'Pathway/Gene'][ko2Met_Table_Final[,'lvl_C_val']=="Adherens junction"] <- "Cell-Cell junction"
    ko2Met_Table_Final[,'Pathway/Gene'][ko2Met_Table_Final[,'lvl_C_val']=="Ferroptosis"] <- "Programmed cell death"
    ko2Met_Table_Final[,'Pathway/Gene'][ko2Met_Table_Final[,'lvl_C_val']=="Apoptosis"] <- "Programmed cell death"
    ko2Met_Table_Final[,'Pathway/Gene'][ko2Met_Table_Final[,'lvl_C_val']=="Necroptosis"] <- "Programmed cell death"
    # Secondary metabolites
    ko2Met_Table_Final[,'Pathway/Gene'][ko2Met_Table_Final[,'lvl_B_val']=="Biosynthesis_of_other_secondary_metabolites"] <- "Secondary metabolites biosynthesis"
    #
    ko2Met_Table_Final_i <- ko2Met_Table_Final %>% dplyr::select('Pathway/Gene',KO_id) %>% distinct()
    ko2Met_Table_Final_j <- ko2Met_Table_Final %>% dplyr::select('Pathway/Gene','lvl_C_val',KO_id) %>% distinct()
    
# Color for ko ------------------------------------------------------------
    color_ko <- ko2Met_Table_Final %>% dplyr::select(lvl_B_val,`Pathway/Gene`) %>% distinct(lvl_B_val,`Pathway/Gene`)
    color_ko[,"palet"] <- ""
    color_ko[,"Label"] <- sapply(strsplit(as.character(color_ko$`Pathway/Gene`),"_-_"), `[`, 1)
    color_ko[,"Label"] <- sapply(str_replace_all(as.character(color_ko$Label),"_"," "), `[`, 1)
    color_ko
    for (i in 1:length(unique(ko2Met_Table_Final$lvl_B_val))) {
      B_courant <- unique(ko2Met_Table_Final$lvl_B_val)[i]
      B_countLabel <- nrow(unique(ko2Met_Table_Final %>% filter(lvl_B_val == B_courant) %>% dplyr::select('Pathway/Gene')))
      #palet_liste <- c("red_material","pink_material","purple_material","deep_purple_material","indigo_material","blue_material","light_blue_material","cyan_material","teal_material","green_material","light_green_material","lime_material","yellow_material","amber_material","orange_material","deep_orange_material","brown_material","grey_material","blue_grey_material")
      palet_liste <- c("red_material","pink_material","deep_purple_material","blue_material","teal_material","green_material","purple_material","lime_material","cyan_material","yellow_material","indigo_material","amber_material","light_blue_material","light_green_material","blue_grey_material","deep_orange_material")
      palet_name <- paste0("ggsci::",palet_liste[i])
      print(palet_name)
      palet_courante <- c(paletteer_d(palet_name,direction = -1,n=10),"black")[seq(1,floor(B_countLabel+(B_countLabel/2))+2,2)]
      print(B_courant)
      print(B_countLabel)
      print(palet_courante)
      for (j in 1:B_countLabel) {
        Label_courant <- unique(color_ko %>% filter(lvl_B_val==B_courant) %>% dplyr::select(Label))[[1]][j]
        print(Label_courant)
        color_ko[,"palet"][color_ko[,"Label"]==Label_courant] <- palet_courante[j]
      }
    }
    color_ko <- color_ko %>% arrange(lvl_B_val)
    color_ko_table <- color_ko %>% dplyr::select(palet,Label) %>% distinct()
    color_ko_table_order <- color_ko_table$Label
#
# Cross & Stat with Funct ---------------------------------------------
  ## Prepare Last column 
    data_seq_funct <- dplyr::left_join(uniq_tax,funct2Last %>% dplyr::select(Last,nbGroup),by = c("Last"="Last"))
    tableVinput_abundance_annots_tax_Funct <- dplyr::left_join(tableVinput_abundance_annots_tax,data_seq_funct %>% dplyr::select(uniq,nbGroup,Last),by = c("Values"= "uniq"))
    tableVinput_abundance_annots_tax_Funct <- separate(tableVinput_abundance_annots_tax_Funct, PR2_taxonomy_complete, c("Domain","Supergroup","Division","Class","Order","Family","Genus","Species"), sep =";")
    interest <- tableVinput_abundance_annots_tax_Funct %>% filter(grepl(ko_BestHitNotSignif,pattern=",")) %>% dplyr::select(Row.names)
    tableVinput_abundance_annots_tax_Funct[,"ko_BestHitNotSignif"][tableVinput_abundance_annots_tax_Funct[,"Row.names"] %in% interest$Row.names] <- NA
    tableVinput_abundance_annots_tax_Funct[,"ko_BestHitNotSignif"][is.na(tableVinput_abundance_annots_tax_Funct[,"ko"]) == FALSE] <- NA
    tableVinput_abundance_annots_tax_Funct <- tableVinput_abundance_annots_tax_Funct %>% tidyr::unite(., col = "Kegg.Onthology",  ko, ko_BestHitNotSignif, na.rm=TRUE, sep = ",",remove=FALSE)
    #tableVinput_abundance_annots_tax_Funct <- separate_rows(tableVinput_abundance_annots_tax_Funct,"Kegg.Onthology",sep = ",") # Duplicate row by KO 
    tableVinput_abundance_annots_tax_Funct[tableVinput_abundance_annots_tax_Funct==""] <- NA
    tableVinput_abundance_annots_tax_Funct[,"nbGroup"][tableVinput_abundance_annots_tax_Funct[,"nbGroup"]=="Unassigned"] <- NA
  #
    nrow(tableVinput_abundance_annots_tax_Funct %>% dplyr::filter(nchar(Kegg.Onthology)>6))
  ## General Grp Stat
    data_seq_summarize <- tableVinput_abundance_annots_tax_Funct %>% dplyr::select(Values,Row.names,ko,PfamDom,ko_BestHitNotSignif,Kegg.Onthology,GOterms,nbGroup)
    #("Lichen","Other multicellular Fungi", "Embryophyceae","Unaffiliated")
    data_seq_summarize[,"nbGroup"][data_seq_summarize[,"nbGroup"] %in% c("Lichen","Other multicellular Fungi", "Embryophyceae","Unaffiliated")] <- NA
    
    #
    Fun_table_TF <- as.data.frame(is.na(data_seq_summarize))
    for (level in c("ko","PfamDom","ko_BestHitNotSignif","Kegg.Onthology","GOterms","nbGroup")) {
      assign(paste("count", level, sep = "_"), nrow(as.data.frame(Fun_table_TF[,level][Fun_table_TF[,level] == FALSE])))}
    count_KOXGRP <- nrow(as.data.frame(Fun_table_TF[,"nbGroup"][Fun_table_TF[,"nbGroup"] == FALSE & Fun_table_TF[,"Kegg.Onthology"] == FALSE]))
    fun_stat <- as.data.frame(c("ko","PfamDom","ko_BestHitNotSignif","Kegg.Onthology","GOterms","nbGroup","KOXGRP")) ; colnames(fun_stat) <- "Level"
    fun_stat$Annotated <- c(count_ko,count_PfamDom,count_ko_BestHitNotSignif,count_Kegg.Onthology,count_GOterms,count_nbGroup,count_KOXGRP)
    fun_stat$`Not annotated` <- nrow(data_seq_summarize) - fun_stat$Annotated
    fun_stat <- reshape2::melt(fun_stat, id = "Level")
  ## Plot
    fun_stat$Level <- factor(fun_stat$Level , levels=c("PfamDom","ko","ko_BestHitNotSignif","Kegg.Onthology","GOterms","nbGroup","KOXGRP"))
    fun_stat$variable <- factor(fun_stat$variable , levels=c("Not annotated","Annotated"))
    fun_plot <- ggplot(fun_stat, mapping = aes(x= Level, y = value, fill = variable, color = variable ,linetype = variable), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
      scale_fill_manual(values = alpha(c("#9d161b","#9d161b"), c(0.4,0.8))) +
      geom_label(aes(y = value + 0.08*max(value),label = paste(value,"\n","unigenes: ", round(value*100/nrow(data_seq_summarize),1)," %",sep =""), color = variable, alpha = variable),size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
      scale_color_manual(values = c("#0000FF00","black")) +
      scale_alpha_manual(values = c(0,0.5)) +
      scale_linetype_manual(values = c("blank","dashed")) +
      guides(fill=guide_legend("Annotation"), color = "none", linetype = "none") +
      labs(x="Annotation type",y="Unigenes") +
      theme(legend.title = element_text(face="bold"),
            axis.title = element_text(color = "black", face = "bold"))
    print(fun_plot)
    ggsave("Stat/Stat-Functional.svg", device = "svg", width = 9, height = 4)
    ## General Grp Stat Abundance
    #
    Fun_table_TF <- tableVinput_abundance_annots_tax_Funct
    Fun_table_TF$Somme <- rowSums(Fun_table_TF %>% dplyr::select(all_of(MetaT_Summary$Echantillon)))
    Fun_table_TF <- Fun_table_TF %>% dplyr::select(Values,Row.names,ko,PfamDom,ko_BestHitNotSignif,Kegg.Onthology,GOterms,nbGroup,Somme)
    
    for (level in c("ko","PfamDom","ko_BestHitNotSignif","Kegg.Onthology","GOterms","nbGroup")) {
      assign(paste("count", level, sep = "_"), sum(as.data.frame(Fun_table_TF[,"Somme"][is.na(Fun_table_TF[,level]) == FALSE])))}
    count_KOXGRP <- colSums(as.data.frame(Fun_table_TF[,"Somme"][is.na(Fun_table_TF[,"nbGroup"]) == FALSE & is.na(Fun_table_TF[,"Kegg.Onthology"]) == FALSE]))
    fun_stat <- as.data.frame(c("ko","PfamDom","ko_BestHitNotSignif","Kegg.Onthology","GOterms","nbGroup","KOXGRP")) ; colnames(fun_stat) <- "Level"
    fun_stat$Annotated <- c(count_ko,count_PfamDom,count_ko_BestHitNotSignif,count_Kegg.Onthology,count_GOterms,count_nbGroup,count_KOXGRP)
    fun_stat$`Not annotated` <- sum(Fun_table_TF$Somme) - fun_stat$Annotated
    fun_stat <- reshape2::melt(fun_stat, id = "Level")
    ## Plot
    fun_stat$Level <- factor(fun_stat$Level , levels=c("PfamDom","ko","ko_BestHitNotSignif","Kegg.Onthology","GOterms","nbGroup","KOXGRP"))
    fun_stat$variable <- factor(fun_stat$variable , levels=c("Not annotated","Annotated"))
    fun_plot <- ggplot(fun_stat, mapping = aes(x= Level, y = value, fill = variable, color = variable ,linetype = variable), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
      scale_fill_manual(values = alpha(c("#9d161b","#9d161b"), c(0.4,0.8))) +
      geom_label(aes(y = value + 0.08*max(value),label = paste(value,"\n","transcripts: ", round(value*100/sum(Fun_table_TF$Somme),1)," %",sep =""), color = variable, alpha = variable),size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
      scale_color_manual(values = c("#0000FF00","black")) +
      scale_alpha_manual(values = c(0,0.5)) +
      scale_linetype_manual(values = c("blank","dashed")) +
      guides(fill=guide_legend("Annotation"), color = "none", linetype = "none") +
      labs(x="Annotation type",y="Transcripts") +
      theme(legend.title = element_text(face="bold"),
            axis.title = element_text(color = "black", face = "bold"))
    print(fun_plot)
    ggsave("Stat/Stat-Functional-Abundance.svg", device = "svg", width = 9, height = 4)
#
      
# Color for Tax ------------------------------------------------------------
    color_tax <- tableVinput_abundance_annots_tax_Funct %>% dplyr::select(Supergroup,Division) %>% distinct(Supergroup,Division)
    color_tax[,"palet"] <- ""
    color_tax[,"Label"] <- color_tax[,"Division"]
    for (i in row.names(color_tax)){
      if (is.na(color_tax[i,"Label"]) == TRUE && is.na(color_tax[i,"Supergroup"])==FALSE) { color_tax[i,"Label"] <- paste("Unaffiliated",color_tax[i,"Supergroup"])}
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
# Color list for DESeq2 ---------------------------------------------------
  # Color for heatmap
    # Specify colors
    ann_colors = list(
      Zone = c(Mixolimnion = "lightgrey", Monimolimnion = "firebrick"),
      Fraction = c(Small = "#1B9E77", Large = "#D95F02"),
      Periods = c(`04` = "#6d6c88", `06` = "#7b8e4b", `09` = "#e3a235", `11` = "#a146a7"),
      Cycle = c(Day = "#a14647", Night = "#6d6ea8"),
      nbGroup = c(HET = "#d84040FF", PARA = "#4080c0FF", `HCOV and SWAT` = "#98d048FF", FLAT = "#70a830FF",SAP = "#ad0635FF", END = "#d7c11cFF", MIXO = "#60d0a0FF"),
      Division = c(`Chlorophyta` = "#01579BFF",
                   `Streptophyta` = "#B3E5FCFF",
                   `Apicomplexa` = "#B71C1CFF",
                   `Choanoflagellida` = "#E65100FF",
                   `Ciliophora` = "#F44336FF",
                   `Cryptophyta`="#FBC02DFF",
                   `Dinoflagellata` = "#E57373FF",
                   `Cercozoa` = "#880E4FFF",
                   `Fungi` = "#F57C00FF",
                   `Haptophyta` = "#FFEB3BFF",
                   `Lobosa` = "#512DA8FF",
                   `Metamonada` = "#689F38FF",
                   `Ochrophyta` = "#006064FF",
                   `Pseudofungi` = "#00BCD4FF",
                   `Unaffiliated Amoebozoa` = "#673AB7FF"),
      Label = c(`Amino acid metabolism` = "#1B5E20FF",
                `Secondary metabolites biosynthesis` = "#F57F17FF",
                `Glycolysis and Gluconeogenesis` = "#B71C1CFF",
                `TCA cycle`="#D32F2FFF",
                `Other Carbohydrate metabolism` = "#F44336FF",
                `Cell cycle` = "#01579BFF",
                `Programmed cell death` = "#0288D1FF",
                `Cellular senescence` = "#03A9F4FF",
                `Focal adhesion` = "#33691EFF",
                `Cell-Cell junction` = "#689F38FF",
                `Regulation of actin cytoskeleton` = "#455A64FF",
                `Oxidative phosphorylation` = "#880E4FFF",
                `Photosynthesis` = "#C2185BFF",
                `Photosynthesis antenna proteins` = "#E91E63FF",
                `Methane metabolism` = "#F06292FF",
                `Nitrogen metabolism` = "#F8BBD0FF",
                `Sulfur metabolism` = "black",
                `Glycan metabolism` = "#4A144CFF",
                `Lipid metabolism` = "#311B92FF",
                `Cofactors and vitamins metabolism` = "#827717FF",
                `Terpenoids and polyketides metabolism` = "#006064FF",
                `Nucleotide metabolism` = "#0D47A1FF",
                `Endocytosis` = "#FF6F00FF",
                `Phagosome` = "#FFA000FF",
                `Lysosome` = "#FFC107FF",
                `Peroxisome` = "#FFD54FFF",
                `Autophagy` = "#FFECB3FF",
                `Xenobiotics biodegradation and metabolism` = "#1A237EFF",
                `Mixed` = "white")
    )
    
# KO analysis -------------------------------------------------------------
  # DEseq2 Summary ------------------------------------------------------------------
    #
    set.seed(123)
    DEseq2_KO_tableVinput <- tableVinput_abundance_annots_tax_Funct %>% dplyr::select(all_of(MetaT_Summary$Echantillon),Kegg.Onthology,nbGroup,PR2_taxonomy) %>% filter(is.na(nbGroup)==FALSE) %>% filter(is.na(Kegg.Onthology)==FALSE) %>% filter(!nbGroup %in% c("Lichen","Other multicellular Fungi", "Embryophyceae","Unaffiliated"))
    DEseq2_KO_tableVinput <- separate_rows(DEseq2_KO_tableVinput,"Kegg.Onthology",sep = ",") # Duplicate row by KO 
    #DEseq2_KO_tableVinput <- DEseq2_KO_tableVinput %>% filter(Kegg.Onthology %in% ko2Met_Table_Final_i$KO_id)
    DEseq2_KO_tableVinput_i <- DEseq2_KO_tableVinput %>% group_by(Kegg.Onthology,nbGroup,PR2_taxonomy) %>% summarise_all(sum)
    # Start ra step
    ra <- letters[seq( from = 1, to = 26 )]
    ra_table <- expand.grid(ra,ra,ra,ra)
    ra_table$Access <- paste0(ra_table$Var1,ra_table$Var2,ra_table$Var3,ra_table$Var4)
    DEseq2_KO_tableVinput_i <- merge(DEseq2_KO_tableVinput_i,ra_table %>% select(Access), by = "row.names")
    DEseq2_KO_table <- DEseq2_KO_tableVinput_i %>% select(-Kegg.Onthology,-nbGroup,-Row.names,-PR2_taxonomy)
    row.names(DEseq2_KO_table) <- DEseq2_KO_table$Access ; DEseq2_KO_table <- DEseq2_KO_table %>% select(-Access)
    # metadata 
    meta <- MetaT_Summary %>% select(-Ref..collab)
    row.names(meta) <- meta$Echantillon
    all(colnames(data) %in% rownames(meta))
    all(colnames(data) == rownames(meta))
    meta$Zone <- as.factor(meta$Zone)
    meta$Fraction <- as.factor(meta$Fraction)
    meta$Periods <- as.factor(meta$Periods)
    meta$Cycle <- as.factor(meta$Cycle)
    # Launch DEseq2
    modele <- ~ Fraction+Periods+Zone
    dds <- DESeqDataSetFromMatrix(countData = DEseq2_KO_table, colData = meta, 
                                  design = formula(modele))
    #
    dds$Fraction <- relevel(dds$Fraction, ref = "Large") ; dds$Fraction
    dds$Zone <- relevel(dds$Zone, ref = "Monimolimnion") ; dds$Zone
    dds$Periods <- relevel(dds$Periods, ref = "04") ; dds$Periods
    #
    dds <- DESeq(dds)
    # Summary figure
    vsd <- vst(dds, blind=FALSE)
    assay_vsd <- assay(vsd)
    sampleDists <- vegdist(t(assay_vsd), method="euclidean")
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(vsd$Zone, vsd$Fraction, vsd$Cycle, vsd$Periods, sep="-")
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    map <- pheatmap(sampleDistMatrix, show_rownames=FALSE,
                    clustering_distance_rows=sampleDists,
                    clustering_distance_cols=sampleDists,
                    annotation_row = meta %>% select(Zone,Fraction, Cycle, Periods),
                    annotation_colors = ann_colors,
                    col=colors)
    svglite("DESeq2/Summary/Heatmapsummary_all_KO.svg",width = 8.00,height = 6.00)
    print(map)
    dev.off()
    #PCA
    #prepare legend
    #color
    legend_color <- plotPCA(vsd, intgroup=c("Zone","Periods")) +
      geom_point(size=6,aes(color = Zone:Periods)) + 
      scale_color_manual("Zone:Periods",values = c("#720500","#AB0800","#CC0A00","#FF584F",
                                    "#003E49","#0091AB","#00D6FC","#8FEDFE")) + 
      theme(legend.title = element_text(face = "bold"))
    legend_color <- get_legend(legend_color)
    #shape
    legend_shape <- plotPCA(vsd, intgroup=c("Fraction")) + guides(color ="none") +
      geom_point(size=6,aes(shape = Fraction),color="black") + 
      theme(legend.title = element_text(face = "bold"))
    legend_shape <- get_legend(legend_shape)
    # Grid
    legend <- plot_grid(legend_color,legend_shape,ncol=1,nrow=2,align = "hv",rel_heights = c(1,1))
    # plot real PCA
    pca <- plotPCA(vsd, intgroup=c("Zone","Fraction","Periods","Cycle")) +
      theme_unique_art() + theme(legend.position = "none") +
      geom_polygon(aes(group=Zone:Fraction:Periods),color ="black") +
      geom_polygon(aes(group=Periods:Zone:Cycle),fill="#FFFFFF00",color ="black") +
      geom_point(size=4.5,aes(shape=Fraction,color = Zone:Fraction:Periods:Cycle)) +
      scale_color_manual(values = c("#720500","#720500","#AB0800","#AB0800","#CC0A00","#CC0A00","#FF584F","#FF584F","#720500","#720500","#AB0800","#AB0800","#CC0A00","#CC0A00","#FF584F","#FF584F",
                                    "#003E49","#003E49","#0091AB","#0091AB","#00D6FC","#00D6FC","#8FEDFE","#8FEDFE","#003E49","#003E49","#0091AB","#0091AB","#00D6FC","#00D6FC","#8FEDFE","#8FEDFE"))
    svglite("DESeq2/Summary/PCAsummary_all_KO.svg",width = 8.00,height = 7.00)
    plot_grid(pca,NULL,legend,nrow = 1,ncol = 3, rel_widths = c(1.5,0.1,0.4))
    dev.off()

    # Save image ------------------------------------------------------------
    save.image("R_DESEQ2_onKO_image.RData")
    save_KO_tableVinput <- tableVinput_abundance_annots_tax_Funct %>% dplyr::select(Row.names,all_of(MetaT_Summary$Echantillon),Values,PR2_taxonomy,nbGroup,Kegg.Onthology,PfamDom,GOterms) 
    save_KO_tableVinput <- separate_rows(save_KO_tableVinput,"Kegg.Onthology",sep = ",") # Duplicate row by KO 
    write.table(save_KO_tableVinput,file="tableVinput_Annot.csv",row.names = F,sep="\t")
#
    # Zone & Fraction analysis ------------------------------------------------
      # Quaternary plot on selected KO ---------------------------------------------------------
        # Prepare input table
        DEseq2_KO_tableVinput <- tableVinput_abundance_annots_tax_Funct %>% dplyr::select(all_of(MetaT_Summary$Echantillon),Kegg.Onthology,nbGroup,PR2_taxonomy) %>% filter(is.na(nbGroup)==FALSE) %>% filter(is.na(Kegg.Onthology)==FALSE) %>% filter(!nbGroup %in% c("Lichen","Other multicellular Fungi", "Embryophyceae","Unaffiliated"))
        DEseq2_KO_tableVinput <- separate_rows(DEseq2_KO_tableVinput,"Kegg.Onthology",sep = ",") # Duplicate row by KO 
        DEseq2_KO_tableVinput <- DEseq2_KO_tableVinput %>% filter(Kegg.Onthology %in% ko2Met_Table_Final_i$KO_id)
        # Taxonomy
        DEseq2_KO_tableVinput[,"Division"] <- sapply(strsplit(as.character(DEseq2_KO_tableVinput$PR2_taxonomy),";"), `[`, 3)
        DEseq2_KO_tableVinput[,"Supergroup"] <- sapply(strsplit(as.character(DEseq2_KO_tableVinput$PR2_taxonomy),";"), `[`, 2)
        DEseq2_KO_tableVinput[,"Class"] <- sapply(strsplit(as.character(DEseq2_KO_tableVinput$PR2_taxonomy),";"), `[`, 4)
        
        for (i in row.names(DEseq2_KO_tableVinput)){
          if (is.na(DEseq2_KO_tableVinput[i,"Division"]) == TRUE && is.na(DEseq2_KO_tableVinput[i,"Supergroup"])==FALSE) { DEseq2_KO_tableVinput[i,"Division"] <- paste("Unaffiliated",DEseq2_KO_tableVinput[i,"Supergroup"])}
          if (grepl(x = DEseq2_KO_tableVinput[i,"Division"], pattern = "_X") == TRUE) { DEseq2_KO_tableVinput[i,"Division"] <- paste("Unaffiliated",DEseq2_KO_tableVinput[i,"Supergroup"])}
          if (grepl(x = DEseq2_KO_tableVinput[i,"Class"], pattern = "_X") == TRUE) { DEseq2_KO_tableVinput[i,"Class"] <- paste("Unaffiliated",DEseq2_KO_tableVinput[i,"Division"])}
          if (is.na(DEseq2_KO_tableVinput[i,"Class"]) == TRUE && grepl(x = DEseq2_KO_tableVinput[i,"Division"], pattern = "Unaffiliated") == FALSE) { DEseq2_KO_tableVinput[i,"Class"] <- paste("Unaffiliated",DEseq2_KO_tableVinput[i,"Division"])}
          if (is.na(DEseq2_KO_tableVinput[i,"Class"]) == TRUE && grepl(x = DEseq2_KO_tableVinput[i,"Division"], pattern = "Unaffiliated") == TRUE) { DEseq2_KO_tableVinput[i,"Class"] <- paste("Unaffiliated",DEseq2_KO_tableVinput[i,"Supergroup"])}
        }
        #
        DEseq2_KO_tableVinput_x <- as.data.frame(DEseq2_KO_tableVinput %>% group_by(Kegg.Onthology,nbGroup,PR2_taxonomy,Class,Supergroup,Division) %>% summarise_all(sum))
        DEseq2_KO_tableVinput_x$Mixolimnion <- rowSums(DEseq2_KO_tableVinput_x %>% select(all_of(ZoneMixolimnion)))
        DEseq2_KO_tableVinput_x$Monimolimnion <- rowSums(DEseq2_KO_tableVinput_x %>% select(all_of(ZoneMonimolimnion)))
        DEseq2_KO_tableVinput_x$Small <- rowSums(DEseq2_KO_tableVinput_x %>% select(all_of(FractionSmall)))
        DEseq2_KO_tableVinput_x$Large <- rowSums(DEseq2_KO_tableVinput_x %>% select(all_of(FractionLarge)))
        DEseq2_KO_tableVinput_x <- DEseq2_KO_tableVinput_x %>% select(-all_of(MetaT_Summary$Echantillon)) %>% arrange(nbGroup)
        write.table(DEseq2_KO_tableVinput_x,file="DESeq2/Heatmap_raw/DEseq2_KO_tableVinput_x.csv",row.names = F,sep="\t")
        
        #
        DEseq2_KO_tableVinput_i <- DEseq2_KO_tableVinput %>% select(-PR2_taxonomy) %>% group_by(Kegg.Onthology,nbGroup,Class,Supergroup,Division) %>% summarise_all(sum)
        # Start ra step
        ra <- letters[seq( from = 1, to = 26 )]
        ra_table <- expand.grid(ra,ra,ra,ra)
        ra_table$Access <- paste0(ra_table$Var1,ra_table$Var2,ra_table$Var3,ra_table$Var4)
        DEseq2_KO_tableVinput_i <- merge(DEseq2_KO_tableVinput_i,ra_table %>% select(Access), by = "row.names")
        DEseq2_KO_table <- DEseq2_KO_tableVinput_i %>% select(-Kegg.Onthology,-nbGroup,-Row.names,-Class,-Division,-Supergroup)
        row.names(DEseq2_KO_table) <- DEseq2_KO_table$Access ; DEseq2_KO_table <- DEseq2_KO_table %>% select(-Access)
        # metadata 
        meta <- MetaT_Summary %>% select(-Ref..collab)
        row.names(meta) <- meta$Echantillon
        all(colnames(data) %in% rownames(meta))
        all(colnames(data) == rownames(meta))
        meta$Zone <- as.factor(meta$Zone)
        meta$Fraction <- as.factor(meta$Fraction)
        meta$Periods <- as.factor(meta$Periods)
        meta$Cycle <- as.factor(meta$Cycle)
        # Launch DEseq2
        modele <- ~ Fraction+Periods+Zone
        dds <- DESeqDataSetFromMatrix(countData = DEseq2_KO_table, colData = meta, 
                                      design = formula(modele))
        #
        dds$Fraction <- relevel(dds$Fraction, ref = "Large") ; dds$Fraction
        dds$Zone <- relevel(dds$Zone, ref = "Monimolimnion") ; dds$Zone
        dds$Periods <- relevel(dds$Periods, ref = "04") ; dds$Periods
        #
        dds <- DESeq(dds)
      # DESEQ2 result
        resultsNames(dds)
        res_zone <- results(dds, contrast= c("Zone","Monimolimnion","Mixolimnion"))
        res_fraction <- results(dds, contrast= c("Fraction","Small","Large"))
        y <- as.data.frame(res_zone[order(res_zone$padj),])
        x <- as.data.frame(res_fraction[order(res_fraction$padj),])
        #
        DEseq2_KO_tableVinput_input <- merge(DEseq2_KO_tableVinput_i %>% select(-Row.names), x %>% select(log2FoldChange,padj,baseMean), by.x = "Access", by.y = "row.names")
        DEseq2_KO_tableVinput_input <- merge(DEseq2_KO_tableVinput_input, y %>% select(log2FoldChange,padj), by.x = "Access", by.y = "row.names")
        # tax
        #DEseq2_KO_tableVinput_input[,"Division"] <- sapply(strsplit(as.character(DEseq2_KO_tableVinput_input$PR2_taxonomy),";"), `[`, 3)
        #DEseq2_KO_tableVinput_input[,"Supergroup"] <- sapply(strsplit(as.character(DEseq2_KO_tableVinput_input$PR2_taxonomy),";"), `[`, 2)
        #for (i in row.names(DEseq2_KO_tableVinput_input)){
        #  if (is.na(DEseq2_KO_tableVinput_input[i,"Division"]) == TRUE && is.na(DEseq2_KO_tableVinput_input[i,"Supergroup"])==FALSE) { DEseq2_KO_tableVinput_input[i,"Division"] <- paste("Unaffiliated",DEseq2_KO_tableVinput_input[i,"Supergroup"])}
        #  if (grepl(x = DEseq2_KO_tableVinput_input[i,"Division"], pattern = "_X") == TRUE) { DEseq2_KO_tableVinput_input[i,"Division"] <- paste("Unaffiliated",DEseq2_KO_tableVinput_input[i,"Supergroup"])}
        #}
        #
        DEseq2_KO_tableVinput_c <- DEseq2_KO_tableVinput_input %>% filter(padj.x < 1e-2 | padj.y < 1e-2) %>% filter(log2FoldChange.x > 2 | log2FoldChange.x < -2 | log2FoldChange.y > 2 | log2FoldChange.y < -2 )
        DEseq2_KO_tableVinput_c <- DEseq2_KO_tableVinput_c %>% filter(baseMean > 10)
        # CDT
        DEseq2_KO_tableVinput_c$Large <- rowSums(DEseq2_KO_tableVinput_c %>% select(all_of(FractionLarge)))
        DEseq2_KO_tableVinput_c$Small <- rowSums(DEseq2_KO_tableVinput_c %>% select(all_of(FractionSmall)))
        DEseq2_KO_tableVinput_c$Mixolimnion <- rowSums(DEseq2_KO_tableVinput_c %>% select(all_of(ZoneMixolimnion)))
        DEseq2_KO_tableVinput_c$Monimolimnion <- rowSums(DEseq2_KO_tableVinput_c %>% select(all_of(ZoneMonimolimnion)))
        #
        #DEseq2_KO_tableVinput_c[,"nbGroup"][DEseq2_KO_tableVinput_c[,"nbGroup"] == "PARA" & DEseq2_KO_tableVinput_c[,"Division"] == "Ochrophyta"] <- "MIXO"
        #
        DEseq2_KO_tableVinput_c[,"TXT"] <- NA
        for (i in row.names(DEseq2_KO_tableVinput_c)) {
          nbGroupi <- DEseq2_KO_tableVinput_c[i,"nbGroup"]
          valuM <- sort((DEseq2_KO_tableVinput_c %>% filter(nbGroup == nbGroupi) %>% select(baseMean))$baseMean, decreasing = TRUE)[10]
          if (DEseq2_KO_tableVinput_c[i,"baseMean"]>=valuM & DEseq2_KO_tableVinput_c[i,"padj.x"]< 1e-2) { DEseq2_KO_tableVinput_c[i,"TXT"] <- DEseq2_KO_tableVinput_c[i,"Kegg.Onthology"]}
          if (DEseq2_KO_tableVinput_c[i,"baseMean"]>=valuM & DEseq2_KO_tableVinput_c[i,"padj.y"]< 1e-2) { DEseq2_KO_tableVinput_c[i,"TXT"] <- DEseq2_KO_tableVinput_c[i,"Kegg.Onthology"]}
        }
        #
        color_tax_i <- color_tax %>% filter(Label %in% unique(DEseq2_KO_tableVinput_c$Division)) %>% filter(is.na(Label)==FALSE)
        order_tax_i <- color_tax_i$Label
        #
        Quaternary_plot <- ggplot(DEseq2_KO_tableVinput_c, mapping = aes(y= log2FoldChange.y, x = log2FoldChange.x)) + theme_unique_art() + 
          theme(legend.text=element_text(size=10),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "right") + 
          geom_point(aes(size=baseMean, color = factor(Division,level = order_tax_i)), shape=21, stroke = 1) +
          scale_size(range = c(2,10),breaks = c(10, 50, 100, 500, 1000)) +
          labs(color="Division", fill = "none", size = "DESeq2 baseMean") + 
          guides(y.sec = "axis", size = guide_legend(title.position = "top",title.hjust = 0.5,ncol=1),
                 color = guide_legend(override.aes = list(shape = 16, size=5),title.position = "top"),
                 fill = "none") + 
          facet_wrap(nbGroup~.) +
          annotate("rect", xmin = -10, xmax = 10, ymin = -10, ymax = 10, color = "black", fill = "#00000000") +#,linetype = "dashed") +
          annotate("rect", xmin = -6, xmax = 6, ymin = -6, ymax = 6, color = "darkgrey", fill = "#00000000",linetype = "dashed") +
          annotate("rect", xmin = -5, xmax = 5, ymin = -5, ymax = 5, color = "grey", fill = "#00000000",linetype = "dashed") +
          annotate("rect", xmin = -4, xmax = 4, ymin = -4, ymax = 4, color = "lightgrey", fill = "#00000000",linetype = "dashed") +
          annotate("text", x = -11, y = 0, label = "Large", angle=90, fontface = "bold") +
          annotate("text", x = 11, y = 0, label = "Small", angle = -90, fontface = "bold") +
          annotate("text", x = 0, y = -11, label = "Mixolimnion", fontface = "bold") +
          annotate("text", x = 0, y = 11, label = "Monimolimnion", fontface = "bold") +
          scale_x_continuous(limits = c(-11,11)) +
          scale_y_reverse(limits = c(11,-11)) +
          geom_label_repel(data = DEseq2_KO_tableVinput_c %>% filter(is.na(TXT) == FALSE),
                           size = 3.2,aes(label = TXT,fill=factor(Division,level = order_tax_i)), 
                           box.padding = unit(0.45, "lines"),nudge_y = 1, nudge_x = 1, alpha = 0.9,
                           max.overlaps = 1000) +
          scale_color_manual(values = color_tax_i$palet) + theme(legend.position = "bottom") + theme(legend.position = c(0.7, 0.15),legend.direction = "horizontal") +
          scale_fill_manual(values = (color_tax_i %>% filter(Division %in% (DEseq2_KO_tableVinput_c %>% filter(is.na(TXT) ==FALSE) %>% select(Division) %>% distinct())$Division) %>% select("palet"))$palet)
        print(Quaternary_plot)
        ggsave("DESeq2/Quaternaryplot/Quaterneryplot_KO_1e-02_log2_2_withlabel_only_10first.svg", device = "svg", width = 11, height = 10)
        #
    #
      # Heatmap on selected KO -----------------------------------------------------------
        # DESEQ2 result
        resultsNames(dds)
        res_zone <- results(dds, contrast= c("Zone","Monimolimnion","Mixolimnion"))
        res_fraction <- results(dds, contrast= c("Fraction","Small","Large"))
        # 
        y <- as.data.frame(res_zone[order(res_zone$padj),])
        x <- as.data.frame(res_fraction[order(res_fraction$padj),])
        #
        vsd <- vst(dds, blind=FALSE)
        vsd <- assay(vsd)
        vsd <- merge(DEseq2_KO_tableVinput_i %>% select(Class,Division,Supergroup,Kegg.Onthology,nbGroup, Access),vsd, by.x = "Access", by.y = "row.names")
        vsd <- merge(vsd, x %>% select(log2FoldChange,padj,baseMean), by.x = "Access", by.y = "row.names")
        vsd <- merge(vsd, y %>% select(log2FoldChange,padj), by.x = "Access", by.y = "row.names")
        
        #vsd[,"Division"] <- sapply(strsplit(as.character(vsd$PR2_taxonomy),";"), `[`, 3)
        #vsd[,"Supergroup"] <- sapply(strsplit(as.character(vsd$PR2_taxonomy),";"), `[`, 2)
        #vsd[,"Class"] <- sapply(strsplit(as.character(vsd$PR2_taxonomy),";"), `[`, 4)
        
        #for (i in row.names(vsd)){
        #  if (is.na(vsd[i,"Division"]) == TRUE && is.na(vsd[i,"Supergroup"])==FALSE) { vsd[i,"Division"] <- paste("Unaffiliated",vsd[i,"Supergroup"])}
        #  if (grepl(x = vsd[i,"Division"], pattern = "_X") == TRUE) { vsd[i,"Division"] <- paste("Unaffiliated",vsd[i,"Supergroup"])}
        #  if (grepl(x = vsd[i,"Class"], pattern = "_X") == TRUE) { vsd[i,"Class"] <- paste("Unaffiliated",vsd[i,"Division"])}
        #  if (is.na(vsd[i,"Class"]) == TRUE && grepl(x = vsd[i,"Division"], pattern = "Unaffiliated") == FALSE) { vsd[i,"Class"] <- paste("Unaffiliated",vsd[i,"Division"])}
        #  if (is.na(vsd[i,"Class"]) == TRUE && grepl(x = vsd[i,"Division"], pattern = "Unaffiliated") == TRUE) { vsd[i,"Class"] <- paste("Unaffiliated",vsd[i,"Supergroup"])}
        #}
        
        DEseq2_KO_tableVinput_heat <- vsd
        #
        row.names(DEseq2_KO_tableVinput_heat) <- paste0(DEseq2_KO_tableVinput_heat$Access,"_",DEseq2_KO_tableVinput_heat$Kegg.Onthology)
        #
        DEseq2_KO_tableVinput_heat$Label <- "Mixed"
        for (i in row.names(DEseq2_KO_tableVinput_heat)) {
          KO <- DEseq2_KO_tableVinput_heat[i,"Kegg.Onthology"]
          Label <- ko2Met_Table_Final_i %>% filter(KO_id == KO)
          if (nrow(Label) > 1) { Labeli <- "Mixed"}
          else {Labeli <- Label$`Pathway/Gene` }
          DEseq2_KO_tableVinput_heat[i,"Label"] <- Labeli
        }
        DEseq2_KO_tableVinput_heat[,"Label"] <- sapply(strsplit(as.character(DEseq2_KO_tableVinput_heat$Label),"_-_"), `[`, 1)
        DEseq2_KO_tableVinput_heat[,"Label"] <- sapply(str_replace_all(as.character(DEseq2_KO_tableVinput_heat$Label),"_"," "), `[`, 1)
        DEseq2_KO_tableVinput_heat <- DEseq2_KO_tableVinput_heat %>% arrange(nbGroup,Label)
        DEseq2_KO_tableVinput_heat_i <- DEseq2_KO_tableVinput_heat
    #
        # Plot log2 > 2 up 10 and more than 5 sharing ko --------------------------------------------------------------------
          # Zone --------------------------------------------------------------------
        df <- DEseq2_KO_tableVinput_heat_i %>% arrange(Label, Kegg.Onthology) %>% filter(padj.y < 5e-2) %>% filter(log2FoldChange.y > 2 | log2FoldChange.y < -2) %>% filter(baseMean > 10) #%>% filter(Label != "Mixed")
        display_DESEQ2 <- reshape2::dcast(df,formula = Class+Division+nbGroup~Kegg.Onthology,fun.aggregate = mean,value.var = "log2FoldChange.y")
        display_DESEQ2[is.na(display_DESEQ2) == TRUE] <- 0
        display_DESEQ2 <- display_DESEQ2 %>% group_by(nbGroup,Division,Class) %>% summarise_all(mean)
        display_DESEQ2 <- as.data.frame(display_DESEQ2)
        row.names(display_DESEQ2) <- paste0(display_DESEQ2$Class,"_",display_DESEQ2$nbGroup)
        vect <- c()
        for (i in colnames(display_DESEQ2[,c(-1,-2,-3)])) { 
          lni <- length(display_DESEQ2[,i][display_DESEQ2[,i] != 0])
          if (lni < 6) { vect <- c(vect,i)}}
        display_DESEQ2 <- display_DESEQ2 %>% select(-all_of(vect))
        pattrow <- display_DESEQ2 %>% select(-"nbGroup",-"Class",-"Division")
        pattrow[pattrow != 0] <- 1
        pattrow <- row.names(pattrow)[rowSums(pattrow)==0]
        display_DESEQ2 <- display_DESEQ2 %>% filter(!row.names(display_DESEQ2) %in% pattrow)
        metaKO <- df %>% select(Kegg.Onthology,Label) %>% filter(!Kegg.Onthology %in% vect) %>% distinct() %>% arrange(Label, Kegg.Onthology)
        # Set color
        colors <- colorRampPalette(c("dodgerblue3", "white", "firebrick3"))(50)
        paletteLength <- 50
        myBreaks <- c(seq(min(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), 0, length.out=ceiling(paletteLength/2) + 1), 
                      seq(max(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology)))/paletteLength, max(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), length.out=floor(paletteLength/2)))
        
        matrixstar <- t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology)))
        matrixstar[as.numeric(matrixstar) <= -6 | as.numeric(matrixstar) >= 6 ] <- "∗∗"
        matrixstar[as.numeric(matrixstar) <= -5 | as.numeric(matrixstar) >= 5 ] <- "∗"
        matrixstar[as.numeric(matrixstar) <= -4 | as.numeric(matrixstar) >= 4 ] <- "⋅"
        matrixstar[is.na(as.numeric(matrixstar)) == FALSE ] <- ""
        
        map <- pheatmap(t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))),
                        cluster_cols=FALSE, cluster_rows=FALSE, show_colnames=TRUE,
                        display_numbers = matrixstar,
                        annotation_row = metaKO %>% select(Label),
                        annotation_col = display_DESEQ2 %>% select(Division),
                        annotation_colors = ann_colors,
                        row_split = metaKO$Label,
                        column_split = display_DESEQ2$nbGroup,
                        color = colors,
                        breaks = myBreaks,
                        border_color = "lightgrey",
                        angle_col = "90")
        heighti <- round(nrow(t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))))*10/48,1)
        if (heighti < 8 & heighti > 4 ) { heighti<- heighti+heighti/1.5 }
        if (heighti < 4) { heighti<- heighti+heighti*1.5 }
        svglite("DESeq2/Heatmap_raw/Heatmap_Zone_KO_5e-02_log2_2_up10_morethan5sharingko.svg",width = 12.00,height = heighti)
        print(map)
        dev.off()
        # Table
        matrix_pheatmap_i <- reshape2::melt(t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))))
        matrix_pheatmap_i <- matrix_pheatmap_i %>% filter(value > 5 | value < -5)
        matrix_pheatmap_i <- separate(matrix_pheatmap_i,Var2,c("Class","Group"),sep="_")
        matrix_pheatmap_i$Zone <- ""
        matrix_pheatmap_i[,"Zone"][matrix_pheatmap_i[,"value"]> 5] <- "Monimolimnion"
        matrix_pheatmap_i[,"Zone"][matrix_pheatmap_i[,"value"]< -5] <- "Mixolimnion"
        matrix_pheatmap_i <- merge(matrix_pheatmap_i,display_DESEQ2 %>% select(Division, Class), by="Class")
        matrix_pheatmap_i <- merge(matrix_pheatmap_i,metaKO, by.x="Var1",by.y="Kegg.Onthology")
        matrix_pheatmap_i <- merge(matrix_pheatmap_i, ko2Met_Table_Interest %>% select(KO_id,KO_val),by.x="Var1",by.y="KO_id")
        matrix_pheatmap_i <- matrix_pheatmap_i %>% distinct() %>% arrange(Label,Group,Zone)
        write.table(matrix_pheatmap_i,file="DESeq2/Heatmap_raw/Heatmap_Zone_5e-02_log2_5_up10_morethan5sharingko.csv",row.names = F,sep="\t")
      #
          # Fraction ----------------------------------------------------------------
        df <- DEseq2_KO_tableVinput_heat_i %>% arrange(Label, Kegg.Onthology) %>% filter(padj.x < 5e-2) %>% filter(log2FoldChange.x > 2 | log2FoldChange.x < -2) %>% filter(baseMean > 10) #%>% filter(Label != "Mixed")
        #df[,"nbGroup"][df[,"nbGroup"] == "PARA" & df[,"Division"] == "Ochrophyta"] <- "MIXO"
        display_DESEQ2 <- reshape2::dcast(df,formula = Class+Division+nbGroup~Kegg.Onthology,fun.aggregate = mean,value.var = "log2FoldChange.x")
        display_DESEQ2[is.na(display_DESEQ2) == TRUE] <- 0
        display_DESEQ2 <- display_DESEQ2 %>% group_by(nbGroup,Division,Class) %>% summarise_all(mean)
        display_DESEQ2 <- as.data.frame(display_DESEQ2)
        row.names(display_DESEQ2) <- paste0(display_DESEQ2$Class,"_",display_DESEQ2$nbGroup)
        vect <- c()
        for (i in colnames(display_DESEQ2[,c(-1,-2,-3)])) { 
          lni <- length(display_DESEQ2[,i][display_DESEQ2[,i] != 0])
          if (lni < 6) { vect <- c(vect,i)}}
        display_DESEQ2 <- display_DESEQ2 %>% select(-all_of(vect))
        metaKO <- df %>% select(Kegg.Onthology,Label) %>% filter(!Kegg.Onthology %in% vect) %>% distinct() %>% arrange(Label, Kegg.Onthology)
        # Set color
        colors <- colorRampPalette(c("dodgerblue3", "white", "firebrick3"))(50)
        paletteLength <- 50
        myBreaks <- c(seq(min(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), 0, length.out=ceiling(paletteLength/2) + 1), 
                      seq(max(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology)))/paletteLength, max(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), length.out=floor(paletteLength/2)))
        
        matrixstar <- t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology)))
        matrixstar[as.numeric(matrixstar) <= -6 | as.numeric(matrixstar) >= 6 ] <- "∗∗"
        matrixstar[as.numeric(matrixstar) <= -5 | as.numeric(matrixstar) >= 5 ] <- "∗"
        matrixstar[as.numeric(matrixstar) <= -4 | as.numeric(matrixstar) >= 4 ] <- "⋅"
        matrixstar[is.na(as.numeric(matrixstar)) == FALSE ] <- ""
        
        map <- pheatmap(t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), 
                        cluster_cols=FALSE, cluster_rows=FALSE, show_colnames=TRUE,
                        annotation_row = metaKO %>% select(Label), 
                        display_numbers = matrixstar,
                        annotation_col = display_DESEQ2 %>% select(Division), 
                        annotation_colors = ann_colors,
                        row_split = metaKO$Label,
                        column_split = display_DESEQ2$nbGroup,
                        color = colors,
                        breaks = myBreaks,
                        border_color = "lightgrey",
                        angle_col = "90")
        heighti <- round(nrow(t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))))*10/48,1)
        if (heighti < 8 & heighti > 4 ) { heighti<- heighti+heighti/1.5 }
        if (heighti < 4) { heighti<- heighti+heighti*1.5 }
        svglite("DESeq2/Heatmap_raw/Heatmap_Fraction_KO_5e-02_log2_4_up10_morethan5sharingko.svg",width = 12.00,height = heighti)
        print(map)
        dev.off()
        # Table
        matrix_pheatmap_i <- reshape2::melt(t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))))
        matrix_pheatmap_i <- matrix_pheatmap_i %>% filter(value >= 5 | value <= -5)
        matrix_pheatmap_i <- separate(matrix_pheatmap_i,Var2,c("Class","Group"),sep="_")
        matrix_pheatmap_i$Fraction <- ""
        matrix_pheatmap_i[,"Fraction"][matrix_pheatmap_i[,"value"]> 5] <- "Small"
        matrix_pheatmap_i[,"Fraction"][matrix_pheatmap_i[,"value"]< -5] <- "Large"
        matrix_pheatmap_i <- merge(matrix_pheatmap_i,display_DESEQ2 %>% select(Division, Class), by="Class")
        matrix_pheatmap_i <- merge(matrix_pheatmap_i,metaKO, by.x="Var1",by.y="Kegg.Onthology")
        matrix_pheatmap_i <- merge(matrix_pheatmap_i, ko2Met_Table_Interest %>% select(KO_id,KO_val),by.x="Var1",by.y="KO_id")
        matrix_pheatmap_i <- matrix_pheatmap_i %>% distinct() %>% arrange(Label,Group,Fraction)
        write.table(matrix_pheatmap_i,file="DESeq2/Heatmap_raw/Heatmap_Fraction_5e-02_log2_5_up10_morethan5sharingko.csv",row.names = F,sep="\t")
        #
        # Plot log2 > 2 all --------------------------------------------------------------------
          # Zone --------------------------------------------------------------------
        df <- DEseq2_KO_tableVinput_heat_i %>% arrange(Label, Kegg.Onthology) %>% filter(padj.y < 5e-2) %>% filter(log2FoldChange.y > 2 | log2FoldChange.y < -2) %>% filter(baseMean > 10) #%>% filter(Label != "Mixed")
        #df[,"nbGroup"][df[,"nbGroup"] == "PARA" & df[,"Division"] == "Ochrophyta"] <- "MIXO"
        display_DESEQ2 <- reshape2::dcast(df,formula = Class+Division+nbGroup~Kegg.Onthology,fun.aggregate = mean,value.var = "log2FoldChange.y")
        display_DESEQ2[is.na(display_DESEQ2) == TRUE] <- 0
        display_DESEQ2 <- display_DESEQ2 %>% group_by(nbGroup,Division,Class) %>% summarise_all(mean)
        display_DESEQ2 <- as.data.frame(display_DESEQ2)
        row.names(display_DESEQ2) <- paste0(display_DESEQ2$Class,"_",display_DESEQ2$nbGroup)
        metaKO <- df %>% select(Kegg.Onthology,Label) %>% distinct() %>% arrange(Label, Kegg.Onthology)
        # Set color
        colors <- colorRampPalette(c("dodgerblue3", "white", "firebrick3"))(50)
        paletteLength <- 50
        myBreaks <- c(seq(min(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), 0, length.out=ceiling(paletteLength/2) + 1), 
                      seq(max(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology)))/paletteLength, max(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), length.out=floor(paletteLength/2)))
        
        matrixstar <- t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology)))
        matrixstar[as.numeric(matrixstar) <= -6 | as.numeric(matrixstar) >= 6 ] <- "∗∗"
        matrixstar[as.numeric(matrixstar) <= -5 | as.numeric(matrixstar) >= 5 ] <- "∗"
        matrixstar[as.numeric(matrixstar) <= -4 | as.numeric(matrixstar) >= 4 ] <- "⋅"
        matrixstar[is.na(as.numeric(matrixstar)) == FALSE ] <- ""
        
        map <- pheatmap(t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), 
                        cluster_cols=FALSE, cluster_rows=FALSE, show_colnames=TRUE,
                        annotation_row = metaKO %>% select(Label), 
                        display_numbers = matrixstar,
                        annotation_col = display_DESEQ2 %>% select(Division), 
                        annotation_colors = ann_colors,
                        row_split = metaKO$Label,
                        column_split = display_DESEQ2$nbGroup,
                        color = colors,
                        breaks = myBreaks,
                        border_color = "lightgrey",
                        angle_col = "90")
        heighti <- round(nrow(t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))))*10/48,1)
        if (heighti < 8 & heighti > 4 ) { heighti<- heighti+heighti/1.5 }
        if (heighti < 4) { heighti<- heighti+heighti*1.5 }
        svglite("DESeq2/Heatmap_raw/Heatmap_Zone_KO_5e-02_log2_2.svg",width = 12.00,height = heighti)
        print(map)
        dev.off()
        # Table
        matrix_pheatmap_i <- reshape2::melt(t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))))
        matrix_pheatmap_i <- matrix_pheatmap_i %>% filter(value > 5 | value < -5)
        matrix_pheatmap_i <- separate(matrix_pheatmap_i,Var2,c("Class","Group"),sep="_")
        matrix_pheatmap_i$Zone <- ""
        matrix_pheatmap_i[,"Zone"][matrix_pheatmap_i[,"value"]> 5] <- "Monimolimnion"
        matrix_pheatmap_i[,"Zone"][matrix_pheatmap_i[,"value"]< -5] <- "Mixolimnion"
        matrix_pheatmap_i <- merge(matrix_pheatmap_i,display_DESEQ2 %>% select(Division, Class), by="Class")
        matrix_pheatmap_i <- merge(matrix_pheatmap_i,metaKO, by.x="Var1",by.y="Kegg.Onthology")
        matrix_pheatmap_i <- merge(matrix_pheatmap_i, ko2Met_Table_Interest %>% select(KO_id,KO_val),by.x="Var1",by.y="KO_id")
        matrix_pheatmap_i <- matrix_pheatmap_i %>% distinct() %>% arrange(Label,Group,Zone)
        write.table(matrix_pheatmap_i,file="DESeq2/Heatmap_raw/Heatmap_Zone_5e-02_log2_5_up10_all.csv",row.names = F,sep="\t")
        #
#
          # Fraction ----------------------------------------------------------------
        df <- DEseq2_KO_tableVinput_heat_i %>% arrange(Label, Kegg.Onthology) %>% filter(padj.x < 5e-2) %>% filter(log2FoldChange.x > 2 | log2FoldChange.x < -2) %>% filter(baseMean > 10) #%>% filter(Label != "Mixed")
        #df[,"nbGroup"][df[,"nbGroup"] == "PARA" & df[,"Division"] == "Ochrophyta"] <- "MIXO"
        display_DESEQ2 <- reshape2::dcast(df,formula = Class+Division+nbGroup~Kegg.Onthology,fun.aggregate = mean,value.var = "log2FoldChange.x")
        display_DESEQ2[is.na(display_DESEQ2) == TRUE] <- 0
        display_DESEQ2 <- display_DESEQ2 %>% group_by(nbGroup,Division,Class) %>% summarise_all(mean)
        display_DESEQ2 <- as.data.frame(display_DESEQ2)
        row.names(display_DESEQ2) <- paste0(display_DESEQ2$Class,"_",display_DESEQ2$nbGroup)
        metaKO <- df %>% select(Kegg.Onthology,Label) %>% distinct() %>% arrange(Label, Kegg.Onthology)
        # Set color
        colors <- colorRampPalette(c("dodgerblue3", "white", "firebrick3"))(50)
        paletteLength <- 50
        myBreaks <- c(seq(min(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), 0, length.out=ceiling(paletteLength/2) + 1), 
                      seq(max(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology)))/paletteLength, max(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), length.out=floor(paletteLength/2)))
        
        matrixstar <- t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology)))
        matrixstar[as.numeric(matrixstar) <= -6 | as.numeric(matrixstar) >= 6 ] <- "∗∗"
        matrixstar[as.numeric(matrixstar) <= -5 | as.numeric(matrixstar) >= 5 ] <- "∗"
        matrixstar[as.numeric(matrixstar) <= -4 | as.numeric(matrixstar) >= 4 ] <- "⋅"
        matrixstar[is.na(as.numeric(matrixstar)) == FALSE ] <- ""
        
        map <- pheatmap(t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), 
                        cluster_cols=FALSE, cluster_rows=FALSE, show_colnames=TRUE,
                        annotation_row = metaKO %>% select(Label), 
                        annotation_col = display_DESEQ2 %>% select(Division), 
                        annotation_colors = ann_colors,
                        display_numbers = matrixstar,
                        row_split = metaKO$Label,
                        column_split = display_DESEQ2$nbGroup,
                        color = colors,
                        breaks = myBreaks,
                        border_color = "lightgrey",
                        angle_col = "90")
        heighti <- round(nrow(t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))))*10/48,1)
        if (heighti < 8 & heighti > 4 ) { heighti<- heighti+heighti/1.5 }
        if (heighti < 4) { heighti<- heighti+heighti*1.5 }
        svglite("DESeq2/Heatmap_raw/Heatmap_Fraction_KO_5e-02_log2_2.svg",width = 12.00,height = heighti)
        print(map)
        dev.off()
        # Table
        matrix_pheatmap_i <- reshape2::melt(t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))))
        matrix_pheatmap_i <- matrix_pheatmap_i %>% filter(value >= 5 | value <= -5)
        matrix_pheatmap_i <- separate(matrix_pheatmap_i,Var2,c("Class","Group"),sep="_")
        matrix_pheatmap_i$Fraction <- ""
        matrix_pheatmap_i[,"Fraction"][matrix_pheatmap_i[,"value"]> 5] <- "Small"
        matrix_pheatmap_i[,"Fraction"][matrix_pheatmap_i[,"value"]< -5] <- "Large"
        matrix_pheatmap_i <- merge(matrix_pheatmap_i,display_DESEQ2 %>% select(Division, Class), by="Class")
        matrix_pheatmap_i <- merge(matrix_pheatmap_i,metaKO, by.x="Var1",by.y="Kegg.Onthology")
        matrix_pheatmap_i <- merge(matrix_pheatmap_i, ko2Met_Table_Interest %>% select(KO_id,KO_val),by.x="Var1",by.y="KO_id")
        matrix_pheatmap_i <- matrix_pheatmap_i %>% distinct() %>% arrange(Label,Group,Fraction)
        write.table(matrix_pheatmap_i,file="DESeq2/Heatmap_raw/Heatmap_Fraction_5e-02_log2_5_up10_all.csv",row.names = F,sep="\t")
        #
#
        # Summarize log2 > 2 all --------------------------------------------------------------------
          # Zone --------------------------------------------------------------------
        df_y <- DEseq2_KO_tableVinput_heat_i %>% arrange(Label, Kegg.Onthology) %>% filter(padj.y < 5e-2) %>% filter(log2FoldChange.y > 2 | log2FoldChange.y < -2) %>% filter(baseMean > 10) #%>% filter(Label != "Mixed")
        display_DESEQ2_y <- reshape2::dcast(df_y,formula = Class+Division+nbGroup~Kegg.Onthology,fun.aggregate = mean,value.var = "log2FoldChange.y")
        display_DESEQ2_y[is.na(display_DESEQ2_y) == TRUE] <- 0
        display_DESEQ2_y <- display_DESEQ2_y %>% group_by(nbGroup,Division,Class) %>% summarise_all(mean)
        display_DESEQ2_y <- as.data.frame(display_DESEQ2_y)
        row.names(display_DESEQ2_y) <- paste0(display_DESEQ2_y$Class,"_",display_DESEQ2_y$nbGroup)
        metaKO <- df_y %>% select(Kegg.Onthology,Label) %>% distinct() %>% arrange(Label, Kegg.Onthology)
        summarize_DESEQ2_y <- t(as.matrix(display_DESEQ2_y %>% select(metaKO$Kegg.Onthology)))
        summarize_DESEQ2_y[summarize_DESEQ2_y > 1] <- 1
        summarize_DESEQ2_y[summarize_DESEQ2_y < -1] <- -1
        summarize_DESEQ2_y <- merge(summarize_DESEQ2_y, metaKO,by.x="row.names",by.y="Kegg.Onthology") %>% arrange(Label)
        row.names(summarize_DESEQ2_y) <- summarize_DESEQ2_y$Row.names
        summarize_DESEQ2_y<- summarize_DESEQ2_y %>% select(-Row.names)
        summarize_DESEQ2_y<- reshape2::melt(summarize_DESEQ2_y,id ="Label")  
        summarize_DESEQ2_y<- summarize_DESEQ2_y  %>% group_by(Label,variable) %>%  
          mutate(Monimolimnion = sum(value>0),
                 Mixolimnion = sum(value<0)) %>% select(-value) %>% distinct()
        #
        patternclgr <- as.data.frame(summarize_DESEQ2_y) %>% select(-Label) %>% group_by(variable) %>% summarise_all(sum)
        patternclgr <- patternclgr %>% filter(Mixolimnion+Monimolimnion<=5) %>% select(variable)
        summarize_DESEQ2_y <- summarize_DESEQ2_y %>% filter(!variable %in% c(patternclgr)$variable)
        #
        summarize_DESEQ2_y <- separate(summarize_DESEQ2_y,variable,c("Class","Group"),"_")
        summarize_DESEQ2_y <- reshape2::melt(summarize_DESEQ2_y,id=c("Class","Group","Label"))
        summarize_DESEQ2_y$mode <- paste(summarize_DESEQ2_y$Group,summarize_DESEQ2_y$Class,sep="_")
        
        # GGplot
        b_Zone <- ggplot(summarize_DESEQ2_y, mapping = aes(x= variable, y = Label, 
                                                  size = value, 
                                                  color = variable), 
                    Rowv = NA, col = colMain, scale = "column") +
          geom_point(stat="identity", aes(group=Label),alpha=0.5) +
          scale_size(range = c(1,25), name = "Abundance") +
          geom_text(mapping = aes(label = value), 
                    size = 2.5,
                    color = "black",
                    fontface = "plain") +
          theme_bw() +
          scale_color_manual(values=c("#ad494aFF","#7375b5FF")) + 
          theme(axis.title = element_text(face="bold", size=12), 
                axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
                axis.title.x = element_text(vjust = -1), 
                title = element_text(face="bold", size=14),
                legend.title = element_text(face="bold"),
                legend.position = "right",
                legend.text = element_text(size=10),
                strip.background = element_rect(color="black", fill="white", size=1, linetype="blank"),
                strip.text = element_text(face="bold")) + 
          ggh4x::facet_nested_wrap(~ Group+Class,nest_line = element_line(linetype = 1),nrow=3) + 
          labs(color = "Zone") + guides(size="none") + guides(color = guide_legend(override.aes = list(size=5))) +
          theme(axis.title = element_blank(),legend.position = "bottom",axis.text.x=element_blank(),axis.ticks.x = element_blank(),title=element_blank())

        svglite("DESeq2/Heatmap_raw/Bubble_Zone_KO_5e-02_log2_2.svg",width = 19.00,height = 11.00)
        print(b_Zone)
        dev.off()
        
          # Fraction --------------------------------------------------------------------
        df_x <- DEseq2_KO_tableVinput_heat_i %>% arrange(Label, Kegg.Onthology) %>% filter(padj.x < 5e-2) %>% filter(log2FoldChange.x > 2 | log2FoldChange.x < -2) %>% filter(baseMean > 10) #%>% filter(Label != "Mixed")
        display_DESEQ2_x <- reshape2::dcast(df_x,formula = Class+Division+nbGroup~Kegg.Onthology,fun.aggregate = mean,value.var = "log2FoldChange.x")
        display_DESEQ2_x[is.na(display_DESEQ2_x) == TRUE] <- 0
        display_DESEQ2_x <- display_DESEQ2_x %>% group_by(nbGroup,Division,Class) %>% summarise_all(mean)
        display_DESEQ2_x <- as.data.frame(display_DESEQ2_x)
        row.names(display_DESEQ2_x) <- paste0(display_DESEQ2_x$Class,"_",display_DESEQ2_x$nbGroup)
        metaKO <- df_x %>% select(Kegg.Onthology,Label) %>% distinct() %>% arrange(Label, Kegg.Onthology)
        summarize_DESEQ2_x <- t(as.matrix(display_DESEQ2_x %>% select(metaKO$Kegg.Onthology)))
        summarize_DESEQ2_x[summarize_DESEQ2_x > 1] <- 1
        summarize_DESEQ2_x[summarize_DESEQ2_x < -1] <- -1
        summarize_DESEQ2_x <- merge(summarize_DESEQ2_x, metaKO,by.x="row.names",by.y="Kegg.Onthology") %>% arrange(Label)
        row.names(summarize_DESEQ2_x) <- summarize_DESEQ2_x$Row.names
        summarize_DESEQ2_x<- summarize_DESEQ2_x %>% select(-Row.names)
        summarize_DESEQ2_x<- reshape2::melt(summarize_DESEQ2_x,id ="Label")  
        summarize_DESEQ2_x<- summarize_DESEQ2_x  %>% group_by(Label,variable) %>%  
          mutate(Small = sum(value>0),
                 Large = sum(value<0)) %>% select(-value) %>% distinct()
        #
        patternclgr <- as.data.frame(summarize_DESEQ2_x) %>% select(-Label) %>% group_by(variable) %>% summarise_all(sum)
        patternclgr <- patternclgr %>% filter(Small+Large<=5) %>% select(variable)
        summarize_DESEQ2_x <- summarize_DESEQ2_x %>% filter(!variable %in% c(patternclgr)$variable)
        #
        summarize_DESEQ2_x <- separate(summarize_DESEQ2_x,variable,c("Class","Group"),"_")
        summarize_DESEQ2_x <- reshape2::melt(summarize_DESEQ2_x,id=c("Class","Group","Label"))
        summarize_DESEQ2_x$mode <- paste(summarize_DESEQ2_x$Group,summarize_DESEQ2_x$Class,sep="_")
        
        # GGplot
        b_Zone <- ggplot(summarize_DESEQ2_x, mapping = aes(x= variable, y = Label, 
                                                         size = value, 
                                                         color = variable), 
                         Rowv = NA, col = colMain, scale = "column") +
          geom_point(stat="identity", aes(group=Label),alpha=0.5) +
          scale_size(range = c(1,25), name = "Abundance") +
          geom_text(mapping = aes(label = value), 
                    size = 2.5,
                    color = "black",
                    fontface = "plain") +
          theme_bw() +
          scale_color_manual(values=c("#ad494aFF","#7375b5FF")) + 
          theme(axis.title = element_text(face="bold", size=12), 
                axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
                axis.title.x = element_text(vjust = -1), 
                title = element_text(face="bold", size=14),
                legend.title = element_text(face="bold"),
                legend.position = "right",
                legend.text = element_text(size=10),
                strip.background = element_rect(color="black", fill="white", size=1, linetype="blank"),
                strip.text = element_text(face="bold")) + 
          ggh4x::facet_nested_wrap(~ Group+Class,nest_line = element_line(linetype = 1),nrow=3) + 
          labs(color = "Fraction") + guides(size="none") + guides(color = guide_legend(override.aes = list(size=5))) +
          theme(axis.title = element_blank(),legend.position = "bottom",axis.text.x=element_blank(),axis.ticks.x = element_blank(),title=element_blank())

        svglite("DESeq2/Heatmap_raw/Bubble_Fraction_KO_5e-02_log2_2.svg",width = 19.00,height = 11.00)
        print(b_Zone)
        dev.off()
        
        
          # Mix -------------------------------------------------------------------
            # Fraction
        df_x <- DEseq2_KO_tableVinput_heat_i %>% arrange(Label, Kegg.Onthology) %>% filter(padj.x < 5e-2) %>% filter(log2FoldChange.x > 2 | log2FoldChange.x < -2) %>% filter(baseMean > 10) #%>% filter(Label != "Mixed")
        display_DESEQ2_x <- reshape2::dcast(df_x,formula = Class+Division+nbGroup~Kegg.Onthology,fun.aggregate = mean,value.var = "log2FoldChange.x")
        display_DESEQ2_x[is.na(display_DESEQ2_x) == TRUE] <- 0
        display_DESEQ2_x <- display_DESEQ2_x %>% group_by(nbGroup,Division,Class) %>% summarise_all(mean)
        display_DESEQ2_x <- as.data.frame(display_DESEQ2_x)
        row.names(display_DESEQ2_x) <- paste0(display_DESEQ2_x$Class,"_",display_DESEQ2_x$nbGroup)
        metaKO <- df_x %>% select(Kegg.Onthology,Label) %>% distinct() %>% arrange(Label, Kegg.Onthology)
        summarize_DESEQ2_x <- t(as.matrix(display_DESEQ2_x %>% select(metaKO$Kegg.Onthology)))
        summarize_DESEQ2_x[summarize_DESEQ2_x > 1] <- 1
        summarize_DESEQ2_x[summarize_DESEQ2_x < -1] <- -1
        summarize_DESEQ2_x <- merge(summarize_DESEQ2_x, metaKO,by.x="row.names",by.y="Kegg.Onthology") %>% arrange(Label)
        row.names(summarize_DESEQ2_x) <- summarize_DESEQ2_x$Row.names
        summarize_DESEQ2_x<- summarize_DESEQ2_x %>% select(-Row.names)
        summarize_DESEQ2_x<- reshape2::melt(summarize_DESEQ2_x,id ="Label")  
        summarize_DESEQ2_x<- summarize_DESEQ2_x  %>% group_by(Label,variable) %>%  
          mutate(Small = sum(value>0),
                 Large = sum(value<0)) %>% select(-value) %>% distinct()
        #
        summarize_DESEQ2_x <- separate(summarize_DESEQ2_x,variable,c("Class","Group"),"_")
        summarize_DESEQ2_x <- reshape2::melt(summarize_DESEQ2_x,id=c("Class","Group","Label"))
        summarize_DESEQ2_x$mode <- "Fraction"
        #Zone
        df_y <- DEseq2_KO_tableVinput_heat_i %>% arrange(Label, Kegg.Onthology) %>% filter(padj.y < 5e-2) %>% filter(log2FoldChange.y > 2 | log2FoldChange.y < -2) %>% filter(baseMean > 10) #%>% filter(Label != "Mixed")
        display_DESEQ2_y <- reshape2::dcast(df_y,formula = Class+Division+nbGroup~Kegg.Onthology,fun.aggregate = mean,value.var = "log2FoldChange.y")
        display_DESEQ2_y[is.na(display_DESEQ2_y) == TRUE] <- 0
        display_DESEQ2_y <- display_DESEQ2_y %>% group_by(nbGroup,Division,Class) %>% summarise_all(mean)
        display_DESEQ2_y <- as.data.frame(display_DESEQ2_y)
        row.names(display_DESEQ2_y) <- paste0(display_DESEQ2_y$Class,"_",display_DESEQ2_y$nbGroup)
        metaKO <- df_y %>% select(Kegg.Onthology,Label) %>% distinct() %>% arrange(Label, Kegg.Onthology)
        summarize_DESEQ2_y <- t(as.matrix(display_DESEQ2_y %>% select(metaKO$Kegg.Onthology)))
        summarize_DESEQ2_y[summarize_DESEQ2_y > 1] <- 1
        summarize_DESEQ2_y[summarize_DESEQ2_y < -1] <- -1
        summarize_DESEQ2_y <- merge(summarize_DESEQ2_y, metaKO,by.x="row.names",by.y="Kegg.Onthology") %>% arrange(Label)
        row.names(summarize_DESEQ2_y) <- summarize_DESEQ2_y$Row.names
        summarize_DESEQ2_y<- summarize_DESEQ2_y %>% select(-Row.names)
        summarize_DESEQ2_y<- reshape2::melt(summarize_DESEQ2_y,id ="Label")  
        summarize_DESEQ2_y<- summarize_DESEQ2_y  %>% group_by(Label,variable) %>%  
          mutate(Monimolimnion = sum(value>0),
                 Mixolimnion = sum(value<0)) %>% select(-value) %>% distinct()
        #
        summarize_DESEQ2_y <- separate(summarize_DESEQ2_y,variable,c("Class","Group"),"_")
        summarize_DESEQ2_y <- reshape2::melt(summarize_DESEQ2_y,id=c("Class","Group","Label"))
        summarize_DESEQ2_y$mode <- "Zone"
        
        # Join
        Mix_summarize_DESEQ2 <- rbind(summarize_DESEQ2_y,summarize_DESEQ2_x)
        Mix_summarize_DESEQ2[is.na(Mix_summarize_DESEQ2)==TRUE] <- 0
        #
        patternclgr <- as.data.frame(Mix_summarize_DESEQ2) %>% select(-Label,-mode,-variable) %>% group_by(Group,Class) %>% summarise_all(sum)
        patternclgr <- patternclgr %>% filter(value<=15) %>% select(Class,Group)
        Mix_summarize_DESEQ2 <- Mix_summarize_DESEQ2 %>% filter(!(Class %in% c(patternclgr)$Class & Group %in% c(patternclgr)$Group))
        #
        Mix_summarize_DESEQ2[Mix_summarize_DESEQ2==0] <- NA
        orderLBL <- rev(unique(sort(Mix_summarize_DESEQ2$Label)))
        orderGrp <- c("END","SAP","FLAT","HCOV and SWAT","HET","MIXO","PARA") 
        # GGplot
        b_Zone <- ggplot(Mix_summarize_DESEQ2, mapping = aes(x= variable, y = factor(Label,levels=orderLBL), 
                                                         size = value, 
                                                         color = variable), 
                         Rowv = NA, col = colMain, scale = "column") +
          geom_point(stat="identity", aes(group=Label),alpha=0.5) +
          scale_size(range = c(1,25), name = "Abundance") +
          geom_text(mapping = aes(label = value), 
                    size = 2.5,
                    color = "black",
                    fontface = "plain") +
          theme_bw() +
          scale_color_manual(values=c("#ad494aFF","#7375b5FF","orange","lightgreen")) + 
          theme(axis.title = element_text(face="bold", size=12), 
                axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
                axis.title.x = element_text(vjust = -1), 
                title = element_text(face="bold", size=14),
                legend.title = element_text(face="bold"),
                legend.position = "right",
                legend.text = element_text(size=10),
                strip.background = element_rect(color="black", fill="white", size=1, linetype="blank"),
                strip.text = element_text(face="bold")) + 
          geom_vline(xintercept = 2.5,linewidth = 0.3,linetype=2) + 
          ggh4x::facet_nested_wrap( ~ factor(Group,levels=orderGrp)+Class,nest_line = element_line(linetype = 1),nrow=2) + 
          labs(color = "Condition") + guides(size="none") + guides(color = guide_legend(override.aes = list(size=5))) +
          theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "bottom",axis.text.x=element_blank(),axis.ticks.x = element_blank(),title=element_blank())
        
        svglite("DESeq2/Heatmap_raw/Bubble_Mix_KO_5e-02_log2_2.svg",width = 16.50,height = 8.00)
        print(b_Zone)
        dev.off()
        
        # GGplot
        b_Zone <- ggplot(Mix_summarize_DESEQ2, mapping = aes(x= factor(Label,levels=rev(orderLBL)), y = variable, 
                                                             size = value, 
                                                             color = variable), 
                         Rowv = NA, col = colMain, scale = "column") +
          geom_point(stat="identity", aes(group=Label),alpha=0.5) +
          scale_size(range = c(1,25), name = "Abundance") +
          geom_text(mapping = aes(label = value), 
                    size = 2.5,
                    color = "black",
                    fontface = "plain") +
          theme_bw() +
          scale_color_manual(values=c("#ad494aFF","#7375b5FF","orange","lightgreen")) + 
          theme(axis.title = element_text(face="bold", size=12), 
                axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0), 
                axis.title.x = element_text(vjust = -1), 
                title = element_text(face="bold", size=14),
                legend.title = element_text(face="bold"),
                legend.position = "right",
                legend.text = element_text(size=10),
                strip.background = element_rect(color="black", fill="white", size=1, linetype="blank"),
                strip.text = element_text(face="bold")) + 
          geom_hline(yintercept = 2.5,linewidth = 0.5,linetype="dotted") + 
          ggh4x::facet_nested_wrap(factor(Group,levels=orderGrp)+Class ~ .,nest_line = element_line(linetype = 1),dir = "v", strip.position = "left",ncol = 2) + 
          labs(color = "Condition") + guides(size="none") + guides(color = guide_legend(override.aes = list(size=5))) +
          theme(axis.title.y = element_blank(),axis.title.x = element_blank(),legend.position = "right",axis.text.y = element_blank(), axis.ticks.y = element_blank(),title=element_blank())

        svglite("DESeq2/Heatmap_raw/Bubble_Mix_KO_5e-02_log2_2_vertical.svg",width = 11.00,height = 17.0)
        print(b_Zone)
        dev.off()
        
      # Periods without interaction analysis ------------------------------------
      # DESEQ2 result
        resultsNames(dds)
        res_periods04_06 <- results(dds, contrast =c("Periods","04","06"))
        res_periods04_09 <- results(dds, contrast =c("Periods","04","09"))
        res_periods04_11 <- results(dds, contrast =c("Periods","04","11"))
        res_periods06_09 <- results(dds, contrast =c("Periods","06","09"))
        res_periods06_11 <- results(dds, contrast =c("Periods","06","11"))
        res_periods09_11 <- results(dds, contrast =c("Periods","09","11"))
        #
        a <- as.data.frame(res_periods04_06[order(res_periods04_06$padj),])
        colnames(a) <- c("baseMean.a","log2FoldChange.a","lfcSE.a","stat.a","pvalue.a","padj.a")
        b <- as.data.frame(res_periods04_09[order(res_periods04_09$padj),])
        colnames(b) <- c("baseMean.b","log2FoldChange.b","lfcSE.b","stat.b","pvalue.b","padj.b")
        c <- as.data.frame(res_periods04_11[order(res_periods04_11$padj),])
        colnames(c) <- c("baseMean.c","log2FoldChange.c","lfcSE.c","stat.c","pvalue.c","padj.c")
        d <- as.data.frame(res_periods06_09[order(res_periods06_09$padj),])
        colnames(d) <- c("baseMean.d","log2FoldChange.d","lfcSE.d","stat.d","pvalue.d","padj.d")
        e <- as.data.frame(res_periods06_11[order(res_periods06_11$padj),])
        colnames(e) <- c("baseMean.e","log2FoldChange.e","lfcSE.e","stat.e","pvalue.e","padj.e")
        f <- as.data.frame(res_periods09_11[order(res_periods09_11$padj),])
        colnames(f) <- c("baseMean.f","log2FoldChange.f","lfcSE.f","stat.f","pvalue.f","padj.f")
        #
        vsd <- vst(dds, blind=FALSE)
        vsd <- assay(vsd)
        vsd <- merge(DEseq2_KO_tableVinput_i %>% select(Class,Division,Supergroup,Kegg.Onthology,nbGroup, Access),vsd, by.x = "Access", by.y = "row.names")
        vsd <- merge(vsd, a %>% select(log2FoldChange.a,padj.a,baseMean.a), by.x = "Access", by.y = "row.names")
        vsd <- merge(vsd, b %>% select(log2FoldChange.b,padj.b,baseMean.b), by.x = "Access", by.y = "row.names")
        vsd <- merge(vsd, c %>% select(log2FoldChange.c,padj.c,baseMean.c), by.x = "Access", by.y = "row.names")
        vsd <- merge(vsd, d %>% select(log2FoldChange.d,padj.d,baseMean.d), by.x = "Access", by.y = "row.names")
        vsd <- merge(vsd, e %>% select(log2FoldChange.e,padj.e,baseMean.e), by.x = "Access", by.y = "row.names")
        vsd <- merge(vsd, f %>% select(log2FoldChange.f,padj.f,baseMean.f), by.x = "Access", by.y = "row.names")
        #
        #vsd[,"Division"] <- sapply(strsplit(as.character(vsd$PR2_taxonomy),";"), `[`, 3)
        #vsd[,"Supergroup"] <- sapply(strsplit(as.character(vsd$PR2_taxonomy),";"), `[`, 2)
        #vsd[,"Class"] <- sapply(strsplit(as.character(vsd$PR2_taxonomy),";"), `[`, 4)
        #
        #for (i in row.names(vsd)){
        #  if (is.na(vsd[i,"Division"]) == TRUE && is.na(vsd[i,"Supergroup"])==FALSE) { vsd[i,"Division"] <- paste("Unaffiliated",vsd[i,"Supergroup"])}
        #  if (grepl(x = vsd[i,"Division"], pattern = "_X") == TRUE) { vsd[i,"Division"] <- paste("Unaffiliated",vsd[i,"Supergroup"])}
        #  if (grepl(x = vsd[i,"Class"], pattern = "_X") == TRUE) { vsd[i,"Class"] <- paste("Unaffiliated",vsd[i,"Division"])}
        #  if (is.na(vsd[i,"Class"]) == TRUE && grepl(x = vsd[i,"Division"], pattern = "Unaffiliated") == FALSE) { vsd[i,"Class"] <- paste("Unaffiliated",vsd[i,"Division"])}
        #  if (is.na(vsd[i,"Class"]) == TRUE && grepl(x = vsd[i,"Division"], pattern = "Unaffiliated") == TRUE) { vsd[i,"Class"] <- paste("Unaffiliated",vsd[i,"Supergroup"])}
        #}
        #
        DEseq2_KO_tableVinput_heat <- vsd #%>% filter(padj.x < 1e-8 | padj.y < 1e-8) %>% filter(log2FoldChange.x > 5 | log2FoldChange.x < -5 | log2FoldChange.y > 5 | log2FoldChange.y < -5 )
        #
        row.names(DEseq2_KO_tableVinput_heat) <- paste0(DEseq2_KO_tableVinput_heat$Access,"_",DEseq2_KO_tableVinput_heat$Kegg.Onthology)
        #
        DEseq2_KO_tableVinput_heat$Label <- "Mixed"
        for (i in row.names(DEseq2_KO_tableVinput_heat)) {
          KO <- DEseq2_KO_tableVinput_heat[i,"Kegg.Onthology"]
          Label <- ko2Met_Table_Final_i %>% filter(KO_id == KO)
          if (nrow(Label) > 1) { Labeli <- "Mixed"}
          else {Labeli <- Label$`Pathway/Gene` }
          DEseq2_KO_tableVinput_heat[i,"Label"] <- Labeli
        }
        DEseq2_KO_tableVinput_heat[,"Label"] <- sapply(strsplit(as.character(DEseq2_KO_tableVinput_heat$Label),"_-_"), `[`, 1)
        DEseq2_KO_tableVinput_heat[,"Label"] <- sapply(str_replace_all(as.character(DEseq2_KO_tableVinput_heat$Label),"_"," "), `[`, 1)
        DEseq2_KO_tableVinput_heat <- DEseq2_KO_tableVinput_heat %>% arrange(nbGroup,Label)
        DEseq2_KO_tableVinput_heat <- DEseq2_KO_tableVinput_heat %>% select(-baseMean.b,-baseMean.c,-baseMean.d,-baseMean.e,-baseMean.f)
        colnames(DEseq2_KO_tableVinput_heat)[colnames(DEseq2_KO_tableVinput_heat) == 'baseMean.a'] <- 'baseMean'
        DEseq2_KO_tableVinput_heat_i <- DEseq2_KO_tableVinput_heat #%>% filter(nbGroup %in% c("PARA")) %>% filter(Division == "Fungi")
        #
        for (i in c("a","b","c","d","e","f")) {
          if (i == "a") { comparaison <- "04vs06" }
          if (i == "b") { comparaison <- "04vs09" }
          if (i == "c") { comparaison <- "04vs11" }
          if (i == "d") { comparaison <- "06vs09" }
          if (i == "e") { comparaison <- "06vs11" }
          if (i == "f") { comparaison <- "09vs11" }
          print(comparaison)
          df <- DEseq2_KO_tableVinput_heat_i %>% arrange(Label, Kegg.Onthology) %>% filter(get(paste("padj",i,sep=".")) < 5e-2) %>% filter(get(paste("log2FoldChange",i,sep=".")) > 2 | get(paste("log2FoldChange",i,sep=".")) < -2) %>% filter(baseMean > 10) #%>% filter(Label != "Mixed")
          #df[,"nbGroup"][df[,"nbGroup"] == "PARA" & df[,"Division"] == "Ochrophyta"] <- "MIXO"
          display_DESEQ2 <- reshape2::dcast(df,formula = Class+Division+nbGroup~Kegg.Onthology,fun.aggregate = mean,value.var = paste("log2FoldChange",i,sep="."))
          display_DESEQ2[is.na(display_DESEQ2) == TRUE] <- 0
          display_DESEQ2 <- display_DESEQ2 %>% group_by(nbGroup,Division,Class) %>% summarise_all(mean)
          display_DESEQ2 <- as.data.frame(display_DESEQ2)
          row.names(display_DESEQ2) <- paste0(display_DESEQ2$Class,"_",display_DESEQ2$nbGroup)
          vect <- c()
          for (j in colnames(display_DESEQ2[,c(-1,-2,-3)])) {
            lni <- length(display_DESEQ2[,j][display_DESEQ2[,j] != 0])
            if (lni < 6) { vect <- c(vect,j)}}
          display_DESEQ2 <- display_DESEQ2 %>% select(-all_of(vect))
          metaKO <- df %>% select(Kegg.Onthology,Label) %>% filter(!Kegg.Onthology %in% vect) %>% distinct() %>% arrange(Label, Kegg.Onthology)
          # Set color
          colors <- colorRampPalette(c("dodgerblue3", "white", "firebrick3"))(50)
          paletteLength <- 50
          myBreaks <- c(seq(min(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), 0, length.out=ceiling(paletteLength/2) + 1), 
                        seq(max(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology)))/paletteLength, max(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), length.out=floor(paletteLength/2)))
          
          matrixstar <- t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology)))
          matrixstar[as.numeric(matrixstar) <= -6 | as.numeric(matrixstar) >= 6 ] <- "∗∗"
          matrixstar[as.numeric(matrixstar) <= -5 | as.numeric(matrixstar) >= 5 ] <- "∗"
          matrixstar[as.numeric(matrixstar) <= -4 | as.numeric(matrixstar) >= 4 ] <- "⋅"
          matrixstar[is.na(as.numeric(matrixstar)) == FALSE ] <- ""
          
          map <- pheatmap(t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), 
                          cluster_cols=FALSE, cluster_rows=FALSE, show_colnames=TRUE,
                          annotation_row = metaKO %>% select(Label), 
                          display_numbers = matrixstar,
                          annotation_col = display_DESEQ2 %>% select(Division), 
                          annotation_colors = ann_colors,
                          row_split = metaKO$Label,
                          column_split = display_DESEQ2$nbGroup,
                          color = colors,
                          breaks = myBreaks,
                          border_color = "lightgrey",
                          angle_col = "90")
          heighti <- round(nrow(t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))))*10/48,1)
          print(heighti)
          if (heighti < 8 & heighti >= 4 ) { heighti<- heighti+heighti/1.5 }
          if (heighti < 4 & heighti >= 2 ) { heighti<- heighti+heighti*1.5 }
          if (heighti < 2 & heighti >= 0.5 ) { heighti<- heighti+heighti*6 }
          if (heighti < 0.5) { heighti <- heighti+heighti*8 }
          
          svglite(paste0("DESeq2/Heatmap_raw/Heatmap_Periods_",comparaison,"_5e-02_log2_2_up10_morethan5sharingko.svg"),width = 12.00,height = heighti)
          print(map)
          dev.off()
          #
          # Table
          matrix_pheatmap_i <- reshape2::melt(t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))))
          matrix_pheatmap_i <- matrix_pheatmap_i %>% filter(value > 6 | value < -6)
          matrix_pheatmap_i <- separate(matrix_pheatmap_i,Var2,c("Class","Group"),sep="_")
          matrix_pheatmap_i$Periods <- ""
          matrix_pheatmap_i[,"Periods"][matrix_pheatmap_i[,"value"]> 5] <- str_split(comparaison,"vs")[[1]][1]
          matrix_pheatmap_i[,"Periods"][matrix_pheatmap_i[,"value"]< -5] <- str_split(comparaison,"vs")[[1]][2]
          matrix_pheatmap_i <- merge(matrix_pheatmap_i,display_DESEQ2 %>% select(Division, Class), by="Class")
          matrix_pheatmap_i <- merge(matrix_pheatmap_i,metaKO, by.x="Var1",by.y="Kegg.Onthology")
          matrix_pheatmap_i <- merge(matrix_pheatmap_i, ko2Met_Table_Interest %>% select(KO_id,KO_val),by.x="Var1",by.y="KO_id")
          matrix_pheatmap_i <- matrix_pheatmap_i %>% distinct() %>% arrange(Label,Group,Periods)
          write.table(matrix_pheatmap_i,file=paste0("DESeq2/Heatmap_raw/Heatmap_Periods_",comparaison,"_5e-02_log2_6_up10_morethan5sharingko.csv"),row.names = F,sep="\t")
          #
          #Bubble
          df <- DEseq2_KO_tableVinput_heat_i %>% arrange(Label, Kegg.Onthology) %>% filter(get(paste("padj",i,sep=".")) < 5e-2) %>% filter(get(paste("log2FoldChange",i,sep=".")) > 2 | get(paste("log2FoldChange",i,sep=".")) < -2) %>% filter(baseMean > 10) #%>% filter(Label != "Mixed")
          #df[,"nbGroup"][df[,"nbGroup"] == "PARA" & df[,"Division"] == "Ochrophyta"] <- "MIXO"
          display_DESEQ2 <- reshape2::dcast(df,formula = Class+Division+nbGroup~Kegg.Onthology,fun.aggregate = mean,value.var = paste("log2FoldChange",i,sep="."))
          display_DESEQ2[is.na(display_DESEQ2) == TRUE] <- 0
          display_DESEQ2 <- display_DESEQ2 %>% group_by(nbGroup,Division,Class) %>% summarise_all(mean)
          display_DESEQ2 <- as.data.frame(display_DESEQ2)
          row.names(display_DESEQ2) <- paste0(display_DESEQ2$Class,"_",display_DESEQ2$nbGroup)
          metaKO <- df %>% select(Kegg.Onthology,Label) %>% distinct() %>% arrange(Label, Kegg.Onthology)
          summarize_DESEQ2 <- t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology)))
          summarize_DESEQ2[summarize_DESEQ2 > 1] <- 1
          summarize_DESEQ2[summarize_DESEQ2 < -1] <- -1
          summarize_DESEQ2 <- merge(summarize_DESEQ2, metaKO,by.x="row.names",by.y="Kegg.Onthology") %>% arrange(Label)
          row.names(summarize_DESEQ2) <- summarize_DESEQ2$Row.names
          summarize_DESEQ2<- summarize_DESEQ2 %>% select(-Row.names)
          summarize_DESEQ2<- reshape2::melt(summarize_DESEQ2,id ="Label")  
          summarize_DESEQ2<- summarize_DESEQ2  %>% group_by(Label,variable) %>%  
            mutate('{str_split(comparaison,"vs")[[1]][1]}' := sum(value>0),
                   '{str_split(comparaison,"vs")[[1]][2]}' := sum(value<0)) %>% select(-value) %>% distinct()
          #
          patternclgr <- as.data.frame(summarize_DESEQ2) %>% select(-Label) %>% group_by(variable) %>% summarise_all(sum)
          patternclgr <- patternclgr %>% filter(get(str_split(comparaison,"vs")[[1]][1])+get(str_split(comparaison,"vs")[[1]][2])<=5) %>% select(variable)
          summarize_DESEQ2 <- summarize_DESEQ2 %>% filter(!variable %in% c(patternclgr)$variable)
          #
          summarize_DESEQ2 <- separate(summarize_DESEQ2,variable,c("Class","Group"),"_")
          summarize_DESEQ2 <- reshape2::melt(summarize_DESEQ2,id=c("Class","Group","Label"))
          summarize_DESEQ2$mode <- paste(summarize_DESEQ2$Group,summarize_DESEQ2$Class,sep="_")
          
          # GGplot
          b_Periods <- ggplot(summarize_DESEQ2, mapping = aes(x= variable, y = Label, 
                                                           size = value, 
                                                           color = variable), 
                           Rowv = NA, col = colMain, scale = "column") +
            geom_point(stat="identity", aes(group=Label),alpha=0.5) +
            scale_size(range = c(1,25), name = "Abundance") +
            geom_text(mapping = aes(label = value), 
                      size = 2.5,
                      color = "black",
                      fontface = "plain") +
            theme_bw() +
            scale_color_manual(values=c("#ad494aFF","#7375b5FF")) + 
            theme(axis.title = element_text(face="bold", size=12), 
                  axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
                  axis.title.x = element_text(vjust = -1), 
                  title = element_text(face="bold", size=14),
                  legend.title = element_text(face="bold"),
                  legend.position = "right",
                  legend.text = element_text(size=10),
                  strip.background = element_rect(color="black", fill="white", size=1, linetype="blank"),
                  strip.text = element_text(face="bold")) + 
            ggh4x::facet_nested_wrap(~ Group+Class,nest_line = element_line(linetype = 1),nrow=3) + 
            labs(color = paste0("Periods",comparaison,sep = ": ")) + guides(size="none") + guides(color = guide_legend(override.aes = list(size=5))) +
            theme(axis.title.x = element_blank(),legend.position = "bottom",axis.text.x=element_blank(),axis.ticks.x = element_blank(),title=element_blank())
          print(b_Periods)
          
          svglite(paste0("DESeq2/Heatmap_raw/Bubble_Periods",comparaison,"_KO_5e-02_log2_2.svg"),width = 19.00,height = 11.00)
          print(b_Periods)
          dev.off()
        }
        # log2foldchnage > 2
        for (i in c("a","b","c","d","e","f")) {
          if (i == "a") { comparaison <- "04vs06" }
          if (i == "b") { comparaison <- "04vs09" }
          if (i == "c") { comparaison <- "04vs11" }
          if (i == "d") { comparaison <- "06vs09" }
          if (i == "e") { comparaison <- "06vs11" }
          if (i == "f") { comparaison <- "09vs11" }
          print(comparaison)
          df <- DEseq2_KO_tableVinput_heat_i %>% arrange(Label, Kegg.Onthology) %>% filter(get(paste("padj",i,sep=".")) < 5e-2) %>% filter(get(paste("log2FoldChange",i,sep=".")) > 2 | get(paste("log2FoldChange",i,sep=".")) < -2) %>% filter(baseMean > 10) #%>% filter(Label != "Mixed")
          #df[,"nbGroup"][df[,"nbGroup"] == "PARA" & df[,"Division"] == "Ochrophyta"] <- "MIXO"
          display_DESEQ2 <- reshape2::dcast(df,formula = Class+Division+nbGroup~Kegg.Onthology,fun.aggregate = mean,value.var = paste("log2FoldChange",i,sep="."))
          display_DESEQ2[is.na(display_DESEQ2) == TRUE] <- 0
          display_DESEQ2 <- display_DESEQ2 %>% group_by(nbGroup,Division,Class) %>% summarise_all(mean)
          display_DESEQ2 <- as.data.frame(display_DESEQ2)
          row.names(display_DESEQ2) <- paste0(display_DESEQ2$Class,"_",display_DESEQ2$nbGroup)
          metaKO <- df %>% select(Kegg.Onthology,Label) %>% distinct() %>% arrange(Label, Kegg.Onthology)
          # Set color
          colors <- colorRampPalette(c("dodgerblue3", "white", "firebrick3"))(50)
          paletteLength <- 50
          myBreaks <- c(seq(min(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), 0, length.out=ceiling(paletteLength/2) + 1), 
                        seq(max(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology)))/paletteLength, max(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), length.out=floor(paletteLength/2)))
          
          matrixstar <- t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology)))
          matrixstar[as.numeric(matrixstar) <= -6 | as.numeric(matrixstar) >= 6 ] <- "∗∗"
          matrixstar[as.numeric(matrixstar) <= -5 | as.numeric(matrixstar) >= 5 ] <- "∗"
          matrixstar[as.numeric(matrixstar) <= -4 | as.numeric(matrixstar) >= 4 ] <- "⋅"
          matrixstar[is.na(as.numeric(matrixstar)) == FALSE ] <- ""
          
          map <- pheatmap(t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), 
                          cluster_cols=FALSE, cluster_rows=FALSE, show_colnames=TRUE,
                          annotation_row = metaKO %>% select(Label), 
                          display_numbers = matrixstar,
                          annotation_col = display_DESEQ2 %>% select(Division), 
                          annotation_colors = ann_colors,
                          row_split = metaKO$Label,
                          column_split = display_DESEQ2$nbGroup,
                          color = colors,
                          breaks = myBreaks,
                          border_color = "lightgrey",
                          angle_col = "90")
          heighti <- round(nrow(t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))))*10/48,1)
          if (heighti < 8 & heighti > 4 ) { heighti<- heighti+heighti/1.5 }
          if (heighti < 4) { heighti<- heighti+heighti*1.5 }
          svglite(paste0("DESeq2/Heatmap_raw/Heatmap_Periods_",comparaison,"_KO_5e-02_log2_2.svg"),width = 12.00,height = heighti)
          print(map)
          dev.off()
          # Table
          matrix_pheatmap_i <- reshape2::melt(t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))))
          matrix_pheatmap_i <- matrix_pheatmap_i %>% filter(value > 2 | value < -2)
          matrix_pheatmap_i <- separate(matrix_pheatmap_i,Var2,c("Class","Group"),sep="_")
          matrix_pheatmap_i$Periods <- ""
          matrix_pheatmap_i[,"Periods"][matrix_pheatmap_i[,"value"]> 2] <- str_split(comparaison,"vs")[[1]][1]
          matrix_pheatmap_i[,"Periods"][matrix_pheatmap_i[,"value"]< -2] <- str_split(comparaison,"vs")[[1]][2]
          matrix_pheatmap_i <- merge(matrix_pheatmap_i,display_DESEQ2 %>% select(Division, Class), by="Class")
          matrix_pheatmap_i <- merge(matrix_pheatmap_i,metaKO, by.x="Var1",by.y="Kegg.Onthology")
          matrix_pheatmap_i <- merge(matrix_pheatmap_i, ko2Met_Table_Interest %>% select(KO_id,KO_val),by.x="Var1",by.y="KO_id")
          matrix_pheatmap_i <- matrix_pheatmap_i %>% distinct() %>% arrange(Group,Class,Label,Periods)
          write.table(matrix_pheatmap_i,file=paste0("DESeq2/Heatmap_raw/Heatmap_Periods_",comparaison,"_5e-02_log2_2_all.csv"),row.names = F,sep="\t")
          #
        }
        
      # Periods summarize log2 > 2 all ------------------------------------------
        Mix_summarize_DESEQ2 <- data.frame() 
        for (i in c("a","d","f")) {
          if (i == "a") { comparaison <- "04vs06" }
          if (i == "b") { comparaison <- "04vs09" }
          if (i == "c") { comparaison <- "04vs11" }
          if (i == "d") { comparaison <- "06vs09" }
          if (i == "e") { comparaison <- "06vs11" }
          if (i == "f") { comparaison <- "09vs11" }
          print(comparaison)
        df <- DEseq2_KO_tableVinput_heat_i %>% arrange(Label, Kegg.Onthology) %>% filter(get(paste("padj",i,sep=".")) < 5e-2) %>% filter(get(paste("log2FoldChange",i,sep=".")) > 2 | get(paste("log2FoldChange",i,sep=".")) < -2) %>% filter(baseMean > 10) #%>% filter(Label != "Mixed")
        display_DESEQ2_x <- reshape2::dcast(df,formula = Class+Division+nbGroup~Kegg.Onthology,fun.aggregate = mean,value.var = paste("log2FoldChange",i,sep="."))
        display_DESEQ2_x[is.na(display_DESEQ2_x) == TRUE] <- 0
        display_DESEQ2_x <- display_DESEQ2_x %>% group_by(nbGroup,Division,Class) %>% summarise_all(mean)
        display_DESEQ2_x <- as.data.frame(display_DESEQ2_x)
        row.names(display_DESEQ2_x) <- paste0(display_DESEQ2_x$Class,"_",display_DESEQ2_x$nbGroup)
        metaKO <- df %>% select(Kegg.Onthology,Label) %>% distinct() %>% arrange(Label, Kegg.Onthology)
        summarize_DESEQ2_x <- t(as.matrix(display_DESEQ2_x %>% select(metaKO$Kegg.Onthology)))
        summarize_DESEQ2_x[summarize_DESEQ2_x > 1] <- 1
        summarize_DESEQ2_x[summarize_DESEQ2_x < -1] <- -1
        summarize_DESEQ2_x <- merge(summarize_DESEQ2_x, metaKO,by.x="row.names",by.y="Kegg.Onthology") %>% arrange(Label)
        row.names(summarize_DESEQ2_x) <- summarize_DESEQ2_x$Row.names
        summarize_DESEQ2_x<- summarize_DESEQ2_x %>% select(-Row.names)
        summarize_DESEQ2_x<- reshape2::melt(summarize_DESEQ2_x,id ="Label")  
        summarize_DESEQ2_x<- summarize_DESEQ2_x  %>% group_by(Label,variable) %>%  
          mutate('{paste0(comparaison,": over-expression in ",str_split(comparaison,"vs")[[1]][1])}' := sum(value>0),
                 '{paste0(comparaison,": over-expression in ",str_split(comparaison,"vs")[[1]][2])}' := sum(value<0)) %>% select(-value) %>% distinct()
        #
        summarize_DESEQ2_x <- separate(summarize_DESEQ2_x,variable,c("Class","Group"),"_")
        summarize_DESEQ2_x <- reshape2::melt(summarize_DESEQ2_x,id=c("Class","Group","Label"))
        summarize_DESEQ2_x$mode <- comparaison
        Mix_summarize_DESEQ2 <- rbind(Mix_summarize_DESEQ2,summarize_DESEQ2_x)
        }
        Mix_summarize_DESEQ2[is.na(Mix_summarize_DESEQ2)==TRUE] <- 0
        #
        patternclgr <- as.data.frame(Mix_summarize_DESEQ2) %>% select(-Label,-mode,-variable) %>% group_by(Group,Class) %>% summarise_all(sum)
        patternclgr <- patternclgr %>% filter(value<=20) %>% select(Class,Group)
        Mix_summarize_DESEQ2 <- Mix_summarize_DESEQ2 %>% filter(!(Class %in% c(patternclgr)$Class & Group %in% c(patternclgr)$Group))
        #
        Mix_summarize_DESEQ2[Mix_summarize_DESEQ2==0] <- NA
        orderLBL <- rev(unique(sort(Mix_summarize_DESEQ2$Label)))
        orderGrp <- c("PARA","FLAT","HCOV and SWAT","HET","MIXO","END","SAP") 
        # GGplot
        b_Zone <- ggplot(Mix_summarize_DESEQ2, mapping = aes(x= variable, y = factor(Label,levels=orderLBL), 
                                                             size = value, 
                                                             color = variable), 
                         Rowv = NA, col = colMain, scale = "column") +
          geom_point(stat="identity", aes(group=Label),alpha=0.5) +
          scale_size(range = c(1,25), name = "Abundance") +
          geom_text(mapping = aes(label = value), 
                    size = 2.5,
                    color = "black",
                    fontface = "plain") +
          theme_bw() +
          scale_color_manual(values=c("#f4817e","#a51015","#6baed6","#25519b","lightgreen","darkgreen","#f6ae6b","darkorange","#d4b9da","#ce1256","grey","#252525")) + 
          theme(axis.title = element_text(face="bold", size=12), 
                axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
                axis.title.x = element_text(vjust = -1), 
                title = element_text(face="bold", size=14),
                legend.title = element_text(face="bold"),
                legend.position = "bottom",
                legend.text = element_text(size=10),
                strip.background = element_rect(color="black", fill="white", size=1, linetype="blank"),
                strip.text = element_text(face="bold")) + 
          geom_vline(xintercept = c(2.5,4.5),linewidth = 0.3,linetype=2) + 
          ggh4x::facet_nested_wrap( ~ factor(Group,levels=orderGrp)+Class,nest_line = element_line(linetype = 1),nrow=2) + 
          labs(color = "Condition") + guides(size="none") + guides(color = guide_legend(override.aes = list(size=5),ncol=3)) +
          theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "bottom",axis.text.x=element_blank(),axis.ticks.x = element_blank(),title=element_blank())
        
        svglite("DESeq2/Heatmap_raw/Bubble_Periods_KO_5e-02_log2_2.svg",width = 18.0,height = 9.00)
        print(b_Zone)
        dev.off()
#
    # Periods and zone interaction effect analysis with DESEQ2 --------------------------------------------------------
      # Periods effect by Zone --------------------------------------------------
  # relaunch DESeq2
    for (levelzone in c("Mixolimnion","Monimolimnion")) {
    modele <- ~ Zone+Fraction+Periods+Zone:Periods
    dds <- DESeqDataSetFromMatrix(countData = DEseq2_KO_table, colData = meta, 
                                  design = formula(modele))
    #
    dds$Fraction <- relevel(dds$Fraction, ref = "Large") ; dds$Fraction
    dds$Zone <- relevel(dds$Zone, ref = levelzone) ; dds$Zone 
    dds$Periods <- relevel(dds$Periods, ref = "04") ; dds$Periods
    #
    dds <- DESeq(dds)
    # DESEQ2 result
    resultsNames(dds)
    res_periods04_06 <- results(dds, contrast =c("Periods","04","06"))
    res_periods04_09 <- results(dds, contrast =c("Periods","04","09"))
    res_periods04_11 <- results(dds, contrast =c("Periods","04","11"))
    res_periods06_09 <- results(dds, contrast =c("Periods","06","09"))
    res_periods06_11 <- results(dds, contrast =c("Periods","06","11"))
    res_periods09_11 <- results(dds, contrast =c("Periods","09","11"))
    #
    a <- as.data.frame(res_periods04_06[order(res_periods04_06$padj),])
    colnames(a) <- c("baseMean.a","log2FoldChange.a","lfcSE.a","stat.a","pvalue.a","padj.a")
    b <- as.data.frame(res_periods04_09[order(res_periods04_09$padj),])
    colnames(b) <- c("baseMean.b","log2FoldChange.b","lfcSE.b","stat.b","pvalue.b","padj.b")
    c <- as.data.frame(res_periods04_11[order(res_periods04_11$padj),])
    colnames(c) <- c("baseMean.c","log2FoldChange.c","lfcSE.c","stat.c","pvalue.c","padj.c")
    d <- as.data.frame(res_periods06_09[order(res_periods06_09$padj),])
    colnames(d) <- c("baseMean.d","log2FoldChange.d","lfcSE.d","stat.d","pvalue.d","padj.d")
    e <- as.data.frame(res_periods06_11[order(res_periods06_11$padj),])
    colnames(e) <- c("baseMean.e","log2FoldChange.e","lfcSE.e","stat.e","pvalue.e","padj.e")
    f <- as.data.frame(res_periods09_11[order(res_periods09_11$padj),])
    colnames(f) <- c("baseMean.f","log2FoldChange.f","lfcSE.f","stat.f","pvalue.f","padj.f")
    #
    vsd <- vst(dds, blind=FALSE)
    vsd <- assay(vsd)
    vsd <- merge(DEseq2_KO_tableVinput_i %>% select(Class,Division,Supergroup,Kegg.Onthology,nbGroup, Access),vsd, by.x = "Access", by.y = "row.names")
    vsd <- merge(vsd, a %>% select(log2FoldChange.a,padj.a,baseMean.a), by.x = "Access", by.y = "row.names")
    vsd <- merge(vsd, b %>% select(log2FoldChange.b,padj.b,baseMean.b), by.x = "Access", by.y = "row.names")
    vsd <- merge(vsd, c %>% select(log2FoldChange.c,padj.c,baseMean.c), by.x = "Access", by.y = "row.names")
    vsd <- merge(vsd, d %>% select(log2FoldChange.d,padj.d,baseMean.d), by.x = "Access", by.y = "row.names")
    vsd <- merge(vsd, e %>% select(log2FoldChange.e,padj.e,baseMean.e), by.x = "Access", by.y = "row.names")
    vsd <- merge(vsd, f %>% select(log2FoldChange.f,padj.f,baseMean.f), by.x = "Access", by.y = "row.names")
    #
    #vsd[,"Division"] <- sapply(strsplit(as.character(vsd$PR2_taxonomy),";"), `[`, 3)
    #vsd[,"Supergroup"] <- sapply(strsplit(as.character(vsd$PR2_taxonomy),";"), `[`, 2)
    #vsd[,"Class"] <- sapply(strsplit(as.character(vsd$PR2_taxonomy),";"), `[`, 4)
    #
    #for (i in row.names(vsd)){
    #  if (is.na(vsd[i,"Division"]) == TRUE && is.na(vsd[i,"Supergroup"])==FALSE) { vsd[i,"Division"] <- paste("Unaffiliated",vsd[i,"Supergroup"])}
    #  if (grepl(x = vsd[i,"Division"], pattern = "_X") == TRUE) { vsd[i,"Division"] <- paste("Unaffiliated",vsd[i,"Supergroup"])}
    #  if (grepl(x = vsd[i,"Class"], pattern = "_X") == TRUE) { vsd[i,"Class"] <- paste("Unaffiliated",vsd[i,"Division"])}
    #  if (is.na(vsd[i,"Class"]) == TRUE && grepl(x = vsd[i,"Division"], pattern = "Unaffiliated") == FALSE) { vsd[i,"Class"] <- paste("Unaffiliated",vsd[i,"Division"])}
    #  if (is.na(vsd[i,"Class"]) == TRUE && grepl(x = vsd[i,"Division"], pattern = "Unaffiliated") == TRUE) { vsd[i,"Class"] <- paste("Unaffiliated",vsd[i,"Supergroup"])}
    #}
    #
    DEseq2_KO_tableVinput_heat <- vsd #%>% filter(padj.x < 1e-8 | padj.y < 1e-8) %>% filter(log2FoldChange.x > 5 | log2FoldChange.x < -5 | log2FoldChange.y > 5 | log2FoldChange.y < -5 )
    #
    row.names(DEseq2_KO_tableVinput_heat) <- paste0(DEseq2_KO_tableVinput_heat$Access,"_",DEseq2_KO_tableVinput_heat$Kegg.Onthology)
    #
    DEseq2_KO_tableVinput_heat$Label <- "Mixed"
    for (i in row.names(DEseq2_KO_tableVinput_heat)) {
      KO <- DEseq2_KO_tableVinput_heat[i,"Kegg.Onthology"]
      Label <- ko2Met_Table_Final_i %>% filter(KO_id == KO)
      if (nrow(Label) > 1) { Labeli <- "Mixed"}
      else {Labeli <- Label$`Pathway/Gene` }
      DEseq2_KO_tableVinput_heat[i,"Label"] <- Labeli
    }
    DEseq2_KO_tableVinput_heat[,"Label"] <- sapply(strsplit(as.character(DEseq2_KO_tableVinput_heat$Label),"_-_"), `[`, 1)
    DEseq2_KO_tableVinput_heat[,"Label"] <- sapply(str_replace_all(as.character(DEseq2_KO_tableVinput_heat$Label),"_"," "), `[`, 1)
    DEseq2_KO_tableVinput_heat <- DEseq2_KO_tableVinput_heat %>% arrange(nbGroup,Label)
    DEseq2_KO_tableVinput_heat <- DEseq2_KO_tableVinput_heat %>% select(-baseMean.b,-baseMean.c,-baseMean.d,-baseMean.e,-baseMean.f)
    colnames(DEseq2_KO_tableVinput_heat)[colnames(DEseq2_KO_tableVinput_heat) == 'baseMean.a'] <- 'baseMean'
    DEseq2_KO_tableVinput_heat_i <- DEseq2_KO_tableVinput_heat #%>% filter(nbGroup %in% c("PARA")) %>% filter(Division == "Fungi")
    #
      # Heatmap Plot multiple comparison ------------------------------------------------
    # log2foldchange > 2 up 100
    for (i in c("a","b","c","d","e","f")) {
      if (i == "a") { comparaison <- "04vs06" }
      if (i == "b") { comparaison <- "04vs09" }
      if (i == "c") { comparaison <- "04vs11" }
      if (i == "d") { comparaison <- "06vs09" }
      if (i == "e") { comparaison <- "06vs11" }
      if (i == "f") { comparaison <- "09vs11" }
    print(comparaison)
    df <- DEseq2_KO_tableVinput_heat_i %>% arrange(Label, Kegg.Onthology) %>% filter(get(paste("padj",i,sep=".")) < 1e-2) %>% filter(get(paste("log2FoldChange",i,sep=".")) > 2 | get(paste("log2FoldChange",i,sep=".")) < -2) %>% filter(baseMean > 10) #%>% filter(Label != "Mixed")
    #df[,"nbGroup"][df[,"nbGroup"] == "PARA" & df[,"Division"] == "Ochrophyta"] <- "MIXO"
    display_DESEQ2 <- reshape2::dcast(df,formula = Class+Division+nbGroup~Kegg.Onthology,fun.aggregate = mean,value.var = paste("log2FoldChange",i,sep="."))
    display_DESEQ2[is.na(display_DESEQ2) == TRUE] <- 0
    display_DESEQ2 <- display_DESEQ2 %>% group_by(nbGroup,Division,Class) %>% summarise_all(mean)
    display_DESEQ2 <- as.data.frame(display_DESEQ2)
    row.names(display_DESEQ2) <- paste0(display_DESEQ2$Class,"_",display_DESEQ2$nbGroup)
    vect <- c()
    for (i in colnames(display_DESEQ2[,c(-1,-2,-3)])) { 
      lni <- length(display_DESEQ2[,i][display_DESEQ2[,i] != 0])
      if (lni < 2) { vect <- c(vect,i)}}
    display_DESEQ2 <- display_DESEQ2 %>% select(-all_of(vect))
    metaKO <- df %>% select(Kegg.Onthology,Label) %>% filter(!Kegg.Onthology %in% vect) %>% distinct() %>% arrange(Label, Kegg.Onthology)
    # Set color
    colors <- colorRampPalette(c("dodgerblue3", "white", "firebrick3"))(50)
    paletteLength <- 50
    myBreaks <- c(seq(min(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology)))/paletteLength, max(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), length.out=floor(paletteLength/2)))
    
    matrixstar <- t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology)))
    matrixstar[as.numeric(matrixstar) <= -6 | as.numeric(matrixstar) >= 6 ] <- "\u2217\u2217"
    matrixstar[as.numeric(matrixstar) <= -5 | as.numeric(matrixstar) >= 5 ] <- "\u2217"
    matrixstar[as.numeric(matrixstar) <= -4 | as.numeric(matrixstar) >= 4 ] <- "\u00B7"
    matrixstar[is.na(as.numeric(matrixstar)) == FALSE ] <- ""
    
    map <- pheatmap(t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), 
                    cluster_cols=FALSE, cluster_rows=FALSE, show_colnames=TRUE,
                    annotation_row = metaKO %>% select(Label), 
                    display_numbers = matrixstar,
                    annotation_col = display_DESEQ2 %>% select(Division), 
                    annotation_colors = ann_colors,
                    row_split = metaKO$Label,
                    column_split = display_DESEQ2$nbGroup,
                    color = colors,
                    breaks = myBreaks,
                    border_color = "lightgrey",
                    angle_col = "90")
    heighti <- round(nrow(t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))))*10/48,1)
    if (heighti < 20) { heighti<- heighti+heighti}
    if (heighti < 8) { heighti<- heighti+heighti/2}
    svglite(paste0("DESeq2/Heatmap_Interaction/Heatmap_Periods_",comparaison,"_in_",levelzone,"_KO_1e-02_log2_2_up100.svg"),width = 12.00,height = heighti)
    print(map)
    dev.off()
    }
    # log2foldchnage > 2
    for (i in c("a","b","c","d","e","f")) {
      if (i == "a") { comparaison <- "04vs06" }
      if (i == "b") { comparaison <- "04vs09" }
      if (i == "c") { comparaison <- "04vs11" }
      if (i == "d") { comparaison <- "06vs09" }
      if (i == "e") { comparaison <- "06vs11" }
      if (i == "f") { comparaison <- "09vs11" }
      print(comparaison)
      df <- DEseq2_KO_tableVinput_heat_i %>% arrange(Label, Kegg.Onthology) %>% filter(get(paste("padj",i,sep=".")) < 1e-2) %>% filter(get(paste("log2FoldChange",i,sep=".")) > 2 | get(paste("log2FoldChange",i,sep=".")) < -2) #%>% filter(baseMean > 10) #%>% filter(Label != "Mixed")
      #df[,"nbGroup"][df[,"nbGroup"] == "PARA" & df[,"Division"] == "Ochrophyta"] <- "MIXO"
      display_DESEQ2 <- reshape2::dcast(df,formula = Class+Division+nbGroup~Kegg.Onthology,fun.aggregate = mean,value.var = paste("log2FoldChange",i,sep="."))
      display_DESEQ2[is.na(display_DESEQ2) == TRUE] <- 0
      display_DESEQ2 <- display_DESEQ2 %>% group_by(nbGroup,Division,Class) %>% summarise_all(mean)
      display_DESEQ2 <- as.data.frame(display_DESEQ2)
      row.names(display_DESEQ2) <- paste0(display_DESEQ2$Class,"_",display_DESEQ2$nbGroup)
      metaKO <- df %>% select(Kegg.Onthology,Label) %>% distinct() %>% arrange(Label, Kegg.Onthology)
      # Set color
      colors <- colorRampPalette(c("dodgerblue3", "white", "firebrick3"))(50)
      paletteLength <- 50
      myBreaks <- c(seq(min(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), 0, length.out=ceiling(paletteLength/2) + 1), 
                    seq(max(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology)))/paletteLength, max(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), length.out=floor(paletteLength/2)))
      
      matrixstar <- t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology)))
      matrixstar[as.numeric(matrixstar) <= -6 | as.numeric(matrixstar) >= 6 ] <- "\u2217\u2217"
      matrixstar[as.numeric(matrixstar) <= -5 | as.numeric(matrixstar) >= 5 ] <- "\u2217"
      matrixstar[as.numeric(matrixstar) <= -4 | as.numeric(matrixstar) >= 4 ] <- "\u00B7"
      matrixstar[is.na(as.numeric(matrixstar)) == FALSE ] <- ""
      
      map <- pheatmap(t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), 
                      cluster_cols=FALSE, cluster_rows=FALSE, show_colnames=TRUE,
                      annotation_row = metaKO %>% select(Label), 
                      display_numbers = matrixstar,
                      annotation_col = display_DESEQ2 %>% select(Division), 
                      annotation_colors = ann_colors,
                      row_split = metaKO$Label,
                      column_split = display_DESEQ2$nbGroup,
                      color = colors,
                      breaks = myBreaks,
                      border_color = "lightgrey",
                      angle_col = "90")
      heighti <- round(nrow(t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))))*10/48,1)
      if (heighti < 20) { heighti<- heighti+heighti}
      if (heighti < 8) { heighti<- heighti+heighti/2}
      svglite(paste0("DESeq2/Heatmap_Interaction/Heatmap_Periods_",comparaison,"_in_",levelzone,"_KO_1e-02_log2_2.svg"),width = 12.00,height = heighti)
      print(map)
      dev.off()
    }
    }
# 
      # Zone by Periods effect --------------------------------------------------
        # relaunch DESeq2
        for (levelperiods in c("04","06","09","11")) {
          modele <- ~ Zone+Fraction+Periods+Zone:Periods
          dds <- DESeqDataSetFromMatrix(countData = DEseq2_KO_table, colData = meta, 
                                        design = formula(modele))
          #
          dds$Fraction <- relevel(dds$Fraction, ref = "Large") ; dds$Fraction
          dds$Zone <- relevel(dds$Zone, ref = "Monimolimnion") ; dds$Zone 
          dds$Periods <- relevel(dds$Periods, ref = levelperiods) ; dds$Periods
          #
          dds <- DESeq(dds)
          # DESEQ2 result
          resultsNames(dds)
          res_ZonebyPeriods <- results(dds, contrast =c("Zone","Monimolimnion","Mixolimnion"))
          #
          a <- as.data.frame(res_ZonebyPeriods[order(res_ZonebyPeriods$padj),])
          colnames(a) <- c("baseMean.a","log2FoldChange.a","lfcSE.a","stat.a","pvalue.a","padj.a")
          #
          vsd <- vst(dds, blind=FALSE)
          vsd <- assay(vsd)
          vsd <- merge(DEseq2_KO_tableVinput_i %>% select(Class,Division,Supergroup,Kegg.Onthology,nbGroup, Access),vsd, by.x = "Access", by.y = "row.names")
          vsd <- merge(vsd, a %>% select(log2FoldChange.a,padj.a,baseMean.a), by.x = "Access", by.y = "row.names")
          #
          #vsd[,"Division"] <- sapply(strsplit(as.character(vsd$PR2_taxonomy),";"), `[`, 3)
          #vsd[,"Supergroup"] <- sapply(strsplit(as.character(vsd$PR2_taxonomy),";"), `[`, 2)
          #vsd[,"Class"] <- sapply(strsplit(as.character(vsd$PR2_taxonomy),";"), `[`, 4)
          #
          #for (i in row.names(vsd)){
          #  if (is.na(vsd[i,"Division"]) == TRUE && is.na(vsd[i,"Supergroup"])==FALSE) { vsd[i,"Division"] <- paste("Unaffiliated",vsd[i,"Supergroup"])}
          #  if (grepl(x = vsd[i,"Division"], pattern = "_X") == TRUE) { vsd[i,"Division"] <- paste("Unaffiliated",vsd[i,"Supergroup"])}
          #  if (grepl(x = vsd[i,"Class"], pattern = "_X") == TRUE) { vsd[i,"Class"] <- paste("Unaffiliated",vsd[i,"Division"])}
          #  if (is.na(vsd[i,"Class"]) == TRUE && grepl(x = vsd[i,"Division"], pattern = "Unaffiliated") == FALSE) { vsd[i,"Class"] <- paste("Unaffiliated",vsd[i,"Division"])}
          #  if (is.na(vsd[i,"Class"]) == TRUE && grepl(x = vsd[i,"Division"], pattern = "Unaffiliated") == TRUE) { vsd[i,"Class"] <- paste("Unaffiliated",vsd[i,"Supergroup"])}
          #}
          #
          DEseq2_KO_tableVinput_heat <- vsd #%>% filter(padj.x < 1e-8 | padj.y < 1e-8) %>% filter(log2FoldChange.x > 5 | log2FoldChange.x < -5 | log2FoldChange.y > 5 | log2FoldChange.y < -5 )
          #
          row.names(DEseq2_KO_tableVinput_heat) <- paste0(DEseq2_KO_tableVinput_heat$Access,"_",DEseq2_KO_tableVinput_heat$Kegg.Onthology)
          #
          DEseq2_KO_tableVinput_heat$Label <- "Mixed"
          for (i in row.names(DEseq2_KO_tableVinput_heat)) {
            KO <- DEseq2_KO_tableVinput_heat[i,"Kegg.Onthology"]
            Label <- ko2Met_Table_Final_i %>% filter(KO_id == KO)
            if (nrow(Label) > 1) { Labeli <- "Mixed"}
            else {Labeli <- Label$`Pathway/Gene` }
            DEseq2_KO_tableVinput_heat[i,"Label"] <- Labeli
          }
          DEseq2_KO_tableVinput_heat[,"Label"] <- sapply(strsplit(as.character(DEseq2_KO_tableVinput_heat$Label),"_-_"), `[`, 1)
          DEseq2_KO_tableVinput_heat[,"Label"] <- sapply(str_replace_all(as.character(DEseq2_KO_tableVinput_heat$Label),"_"," "), `[`, 1)
          DEseq2_KO_tableVinput_heat <- DEseq2_KO_tableVinput_heat %>% arrange(nbGroup,Label)
          DEseq2_KO_tableVinput_heat <- DEseq2_KO_tableVinput_heat
          colnames(DEseq2_KO_tableVinput_heat)[colnames(DEseq2_KO_tableVinput_heat) == 'baseMean.a'] <- 'baseMean'
          DEseq2_KO_tableVinput_heat_i <- DEseq2_KO_tableVinput_heat #%>% filter(nbGroup %in% c("PARA")) %>% filter(Division == "Fungi")
          #
          # Heatmap Plot multiple comparison ------------------------------------------------
            # log2foldchange > 2 up100
            print(paste(levelperiods,"with log2foldchange > 2 up100"))
            df <- DEseq2_KO_tableVinput_heat_i %>% arrange(Label, Kegg.Onthology) %>% filter(padj.a < 1e-2) %>% filter(log2FoldChange.a > 2 | log2FoldChange.a < -2) %>% filter(baseMean > 10) #%>% filter(Label != "Mixed")
            #df[,"nbGroup"][df[,"nbGroup"] == "PARA" & df[,"Division"] == "Ochrophyta"] <- "MIXO"
            display_DESEQ2 <- reshape2::dcast(df,formula = Class+Division+nbGroup~Kegg.Onthology,fun.aggregate = mean,value.var = "log2FoldChange.a")
            display_DESEQ2[is.na(display_DESEQ2) == TRUE] <- 0
            display_DESEQ2 <- display_DESEQ2 %>% group_by(nbGroup,Division,Class) %>% summarise_all(mean)
            display_DESEQ2 <- as.data.frame(display_DESEQ2)
            row.names(display_DESEQ2) <- paste0(display_DESEQ2$Class,"_",display_DESEQ2$nbGroup)
            vect <- c()
            for (i in colnames(display_DESEQ2[,c(-1,-2,-3)])) { 
              lni <- length(display_DESEQ2[,i][display_DESEQ2[,i] != 0])
              if (lni < 4) { vect <- c(vect,i)}}
            display_DESEQ2 <- display_DESEQ2 %>% select(-all_of(vect))
            metaKO <- df %>% select(Kegg.Onthology,Label) %>% filter(!Kegg.Onthology %in% vect) %>% distinct() %>% arrange(Label, Kegg.Onthology)
            
            # Set color
            colors <- colorRampPalette(c("dodgerblue3", "white", "firebrick3"))(50)
            paletteLength <- 50
            myBreaks <- c(seq(min(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), 0, length.out=ceiling(paletteLength/2) + 1), 
                          seq(max(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology)))/paletteLength, max(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), length.out=floor(paletteLength/2)))
            
            map <- pheatmap(t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), 
                            cluster_cols=FALSE, cluster_rows=FALSE, show_colnames=TRUE,
                            annotation_row = metaKO %>% select(Label), 
                            annotation_col = display_DESEQ2 %>% select(Division), 
                            annotation_colors = ann_colors,
                            row_split = metaKO$Label,
                            column_split = display_DESEQ2$nbGroup,
                            color = colors,
                            breaks = myBreaks,
                            border_color = "lightgrey",
                            angle_col = "90")
            svglite(paste0("DESeq2/Heatmap_Interaction/Heatmap_ZonebyPeriods_in_",levelperiods,"_KO_1e-02_log2_2_up100.svg"),width = 12.00,height = 40.00)
            print(map)
            dev.off()
            # log2foldchange = all
            print(paste(levelperiods,"with log2foldchange > 2"))
            df <- DEseq2_KO_tableVinput_heat_i %>% arrange(Label, Kegg.Onthology) %>% filter(padj.a < 1e-2) %>% filter(log2FoldChange.a > 2 | log2FoldChange.a < -2) #%>% filter(baseMean > 10) #%>% filter(Label != "Mixed")
            #df[,"nbGroup"][df[,"nbGroup"] == "PARA" & df[,"Division"] == "Ochrophyta"] <- "MIXO"
            display_DESEQ2 <- reshape2::dcast(df,formula = Class+Division+nbGroup~Kegg.Onthology,fun.aggregate = mean,value.var = "log2FoldChange.a")
            display_DESEQ2[is.na(display_DESEQ2) == TRUE] <- 0
            display_DESEQ2 <- display_DESEQ2 %>% group_by(nbGroup,Division,Class) %>% summarise_all(mean)
            display_DESEQ2 <- as.data.frame(display_DESEQ2)
            row.names(display_DESEQ2) <- paste0(display_DESEQ2$Class,"_",display_DESEQ2$nbGroup)
            metaKO <- df %>% select(Kegg.Onthology,Label) %>% distinct() %>% arrange(Label, Kegg.Onthology)
            # Set color
            colors <- colorRampPalette(c("dodgerblue3", "white", "firebrick3"))(50)
            paletteLength <- 50
            myBreaks <- c(seq(min(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), 0, length.out=ceiling(paletteLength/2) + 1), 
                          seq(max(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology)))/paletteLength, max(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), length.out=floor(paletteLength/2)))
            
            map <- pheatmap(t(as.matrix(display_DESEQ2 %>% select(metaKO$Kegg.Onthology))), 
                            cluster_cols=FALSE, cluster_rows=FALSE, show_colnames=TRUE,
                            annotation_row = metaKO %>% select(Label), 
                            annotation_col = display_DESEQ2 %>% select(Division), 
                            annotation_colors = ann_colors,
                            row_split = metaKO$Label,
                            column_split = display_DESEQ2$nbGroup,
                            color = colors,
                            breaks = myBreaks,
                            border_color = "lightgrey",
                            angle_col = "90")
            svglite(paste0("DESeq2/Heatmap_Interaction/Heatmap_ZonebyPeriods_in_",levelperiods,"_KO_1e-02_log2_2.svg"),width = 12.00,height = 140.00)
            print(map)
            dev.off()
        }
        # 
# PCoA on KO --------------------------------------------------------------
  # ZoneXFractionXCycle PCoA -----------------------------------------------------------
        # Start
        set.seed(123)
        DEseq2_PcoA <- tableVinput_abundance_annots_tax_Funct %>% dplyr::select(all_of(MetaT_Summary$Echantillon),Kegg.Onthology,nbGroup) %>% filter(is.na(nbGroup)==FALSE) %>% filter(is.na(Kegg.Onthology)==FALSE) %>% filter(!nbGroup %in% c("Lichen","Other multicellular Fungi", "Embryophyceae","Unaffiliated"))
        DEseq2_PcoA <- separate_rows(DEseq2_PcoA,"Kegg.Onthology",sep = ",") # Duplicate row by KO 
        DEseq2_PcoA <- DEseq2_PcoA %>% group_by(Kegg.Onthology,nbGroup) %>% summarise_all(sum)
        # 
        DEseq2_PcoA_i <- reshape2::melt(DEseq2_PcoA,id.vars = c("nbGroup","Kegg.Onthology"))
        DEseq2_PcoA_i <- merge(DEseq2_PcoA_i,MetaT_Summary %>% select(Echantillon,Fraction,Zone,Cycle,Periods),by.x = "variable",by.y="Echantillon")
        DEseq2_PcoA_i$variable <- paste(DEseq2_PcoA_i$Fraction,DEseq2_PcoA_i$Zone,DEseq2_PcoA_i$Cycle,DEseq2_PcoA_i$Periods,sep="X") ; DEseq2_PcoA_i <- DEseq2_PcoA_i %>% select(-Fraction,-Zone,-Cycle,-Periods)
        DEseq2_PcoA_i$CDT <- paste(DEseq2_PcoA_i$nbGroup,DEseq2_PcoA_i$variable,sep=" in ")
        dataCDT <- reshape2::dcast(data = DEseq2_PcoA_i,formula = Kegg.Onthology~CDT,fun.aggregate = sum,value.var = "value")
        row.names(dataCDT) <- dataCDT$Kegg.Onthology ; dataCDT <- dataCDT %>% select(-Kegg.Onthology)
    # Plot --------------------------------------------------------------------
      # Abondance
        pcoa <- ape::pcoa(vegdist(t(dataCDT), method = "bray"))
        coord <- pcoa$vectors[,1:3] ; colnames(coord) <- c("PCoA 1","PCoA 2","PCoA 3")
        Dim1Seq <- paste0(round(pcoa$values[1,"Relative_eig"],3)*100,"%")
        Dim2Seq <- paste0(round(pcoa$values[2,"Relative_eig"],3)*100,"%")
        Dim3Seq <- paste0(round(pcoa$values[3,"Relative_eig"],3)*100,"%")
        coord <- as.data.frame(coord)
        coord$CDT <- row.names(coord)
        coord[,"nbGroup"] <- sapply(strsplit(as.character(coord$CDT)," in "), `[`, 1)
        coord[,"variable"] <- sapply(strsplit(as.character(coord$CDT)," in "), `[`, 2)
        coord <- separate(coord,"variable",c("Fraction","Zone","Cycle","Periods"),sep="X")
        # ggplot2
        palet_tree_grp <- c("#ad0635FF","#70a830FF","#98d048FF","#d84040FF","#60d0a0FF","#4080c0FF","#d7c11cFF","#40a078FF","#333b43FF")
        pcoa_plot_c <- ggplot(coord) + theme_unique_art() + 
          theme(legend.text=element_text(size=10),legend.position = "right") + 
          geom_point(aes(y= `PCoA 2`, x = `PCoA 1`,color = nbGroup, fill = Cycle, shape = Zone, size = Fraction), stroke =2) + 
          #geom_circle(aes(x0 = `PCoA 1`, y0 = `PCoA 3`, fill= Fraction, r=0.005), colour = "black")  + 
          scale_shape_manual(values = c(21, 24)) +
          scale_size_manual(values = c(4, 2)) +
          scale_color_manual(values = palet_tree_grp) + 
          scale_fill_manual(values = c("white","black")) + 
          labs(title="Abundance",color="Functional groups", shape = "Zone (shape)", size = "Fraction (size)", fill = "Cycle (fill)", x= paste0("Dim 1 [",Dim1Seq,"]"), y=paste0("Dim2 [",Dim2Seq,"]")) +
          facet_grid(Periods~.) +
          guides(fill = guide_legend(override.aes = list(shape = 21, size = 3)), shape = guide_legend(override.aes = list(size = 3))) +
          scale_y_continuous(position = 'left')
        print(pcoa_plot_c)
        LEGEND_c <- get_legend(pcoa_plot_c)
        pcoa_plot_c <- pcoa_plot_c + theme(legend.position = "none")
        #
      # Présence absence
        dataCDT[dataCDT>1] <- 1
        pcoa <- ape::pcoa(vegdist(t(dataCDT), method = "bray"))
        coord <- pcoa$vectors[,1:3] ; colnames(coord) <- c("PCoA 1","PCoA 2","PCoA 3")
        Dim1Seq <- paste0(round(pcoa$values[1,"Relative_eig"],3)*100,"%")
        Dim2Seq <- paste0(round(pcoa$values[2,"Relative_eig"],3)*100,"%")
        Dim3Seq <- paste0(round(pcoa$values[3,"Relative_eig"],3)*100,"%")
        coord <- as.data.frame(coord)
        coord$CDT <- row.names(coord)
        coord[,"nbGroup"] <- sapply(strsplit(as.character(coord$CDT)," in "), `[`, 1)
        coord[,"variable"] <- sapply(strsplit(as.character(coord$CDT)," in "), `[`, 2)
        coord <- separate(coord,"variable",c("Fraction","Zone","Cycle","Periods"),sep="X")
        # ggplot2
        palet_tree_grp <- c("#ad0635FF","#70a830FF","#98d048FF","#d84040FF","#60d0a0FF","#4080c0FF","#d7c11cFF","#40a078FF","#333b43FF")
        pcoa_plot_a <- ggplot(coord) + theme_unique_art() + 
          theme(legend.text=element_text(size=10),legend.position = "right") + 
          geom_point(aes(y= `PCoA 2`, x = `PCoA 1`,color = nbGroup, fill = Cycle, shape = Zone, size = Fraction), stroke =2) + 
          #geom_circle(aes(x0 = `PCoA 1`, y0 = `PCoA 3`, fill= Fraction, r=0.005), colour = "black")  + 
          scale_shape_manual(values = c(21, 24)) +
          scale_size_manual(values = c(4, 2)) +
          scale_color_manual(values = palet_tree_grp) + 
          scale_fill_manual(values = c("white","black")) + 
          labs(title="Presence/absence",color="Functional groups", shape = "Zone (shape)", size = "Fraction (size)", fill = "Cycle (fill)", x= paste0("Dim 1 [",Dim1Seq,"]"), y=paste0("Dim2 [",Dim2Seq,"]")) +
          facet_grid(Periods~.) +
          guides(fill = guide_legend(override.aes = list(shape = 21, size = 3)), shape = guide_legend(override.aes = list(size = 3))) +
          scale_y_continuous(position = 'left') + theme(legend.position = "none")
        print(pcoa_plot_a)
        #
        # Coplot
        svglite("DESeq2/PCoA/GRPandCDTS_PCoA_bray_noNormalization.svg",width = 12.00,height = 9.00)
        b_plot <- plot_grid(pcoa_plot_a,pcoa_plot_c,LEGEND_c, ncol = 3, nrow = 1, rel_widths = c(3,3,2),rel_heights = c(3))
        print(b_plot)
        dev.off()
        #
  # ZoneXFractionXCycle DESEQ2 PCoA -----------------------------------------------------------
    #Launch DESEQ2
        #
        DEseq2_KO_tableVinput <- tableVinput_abundance_annots_tax_Funct %>% dplyr::select(all_of(MetaT_Summary$Echantillon),Kegg.Onthology,nbGroup) %>% filter(is.na(nbGroup)==FALSE) %>% filter(is.na(Kegg.Onthology)==FALSE) %>% filter(!nbGroup %in% c("Lichen","Other multicellular Fungi", "Embryophyceae","Unaffiliated"))
        DEseq2_KO_tableVinput <- separate_rows(DEseq2_KO_tableVinput,"Kegg.Onthology",sep = ",") # Duplicate row by KO 
        DEseq2_KO_tableVinput_i <- DEseq2_KO_tableVinput %>% group_by(Kegg.Onthology,nbGroup) %>% summarise_all(sum)
        # Start ra step
        ra <- letters[seq( from = 1, to = 26 )]
        ra_table <- expand.grid(ra,ra,ra,ra)
        ra_table$Access <- paste0(ra_table$Var1,ra_table$Var2,ra_table$Var3,ra_table$Var4)
        DEseq2_KO_tableVinput_i <- merge(DEseq2_KO_tableVinput_i,ra_table %>% select(Access), by = "row.names")
        DEseq2_KO_table <- DEseq2_KO_tableVinput_i %>% select(-Kegg.Onthology,-nbGroup,-Row.names)
        row.names(DEseq2_KO_table) <- DEseq2_KO_table$Access ; DEseq2_KO_table <- DEseq2_KO_table %>% select(-Access)
        # Start DESEQ2
        modele <- ~ Zone+Fraction+Periods
        dds <- DESeqDataSetFromMatrix(countData = DEseq2_KO_table, colData = meta, 
                                      design = formula(modele))
        #
        dds$Fraction <- relevel(dds$Fraction, ref = "Large") ; dds$Fraction
        dds$Zone <- relevel(dds$Zone, ref = "Monimolimnion") ; dds$Zone 
        dds$Periods <- relevel(dds$Periods, ref = "04") ; dds$Periods
        #
        dds <- DESeq(dds)
        #Normalization
        #vsd <- vst(dds, blind=FALSE)
        #normalized_dds <- assay(vsd)
        dds <- estimateSizeFactors(dds)
        normalized_dds <- counts(dds, normalized=TRUE)
        #
        DEseq2_PCoA <- merge(DEseq2_KO_tableVinput_i %>% select(Kegg.Onthology,nbGroup, Access),normalized_dds, by.x = "Access", by.y = "row.names")
        #
        DEseq2_PCoA_i <- DEseq2_PCoA %>% select(-Access)
        #
        DEseq2_PCoA_i <- reshape2::melt(DEseq2_PCoA_i,id.vars = c("nbGroup","Kegg.Onthology"))
        DEseq2_PCoA_i <- merge(DEseq2_PCoA_i,MetaT_Summary %>% select(Echantillon,Fraction,Zone,Cycle,Periods),by.x = "variable",by.y="Echantillon")
        DEseq2_PCoA_i$variable <- paste(DEseq2_PCoA_i$Fraction,DEseq2_PCoA_i$Zone,DEseq2_PCoA_i$Cycle,DEseq2_PCoA_i$Periods,sep="X") ; DEseq2_PCoA_i <- DEseq2_PCoA_i %>% select(-Fraction,-Zone,-Cycle,-Periods)
        DEseq2_PCoA_i$CDT <- paste(DEseq2_PCoA_i$nbGroup,DEseq2_PCoA_i$variable,sep=" in ")
        dataCDT <- reshape2::dcast(data = DEseq2_PCoA_i,formula = Kegg.Onthology~CDT,fun.aggregate = sum,value.var = "value")
        row.names(dataCDT) <- dataCDT$Kegg.Onthology ; dataCDT <- dataCDT %>% select(-Kegg.Onthology)
    # plot --------------------------------------------------------------------
      # Abondance
        pcoa <- ape::pcoa(vegdist(t(dataCDT), method = "bray"))
        coord <- pcoa$vectors[,1:3] ; colnames(coord) <- c("PCoA 1","PCoA 2","PCoA 3")
        Dim1Seq <- paste0(round(pcoa$values[1,"Relative_eig"],3)*100,"%")
        Dim2Seq <- paste0(round(pcoa$values[2,"Relative_eig"],3)*100,"%")
        Dim3Seq <- paste0(round(pcoa$values[3,"Relative_eig"],3)*100,"%")
        coord <- as.data.frame(coord)
        coord$CDT <- row.names(coord)
        coord[,"nbGroup"] <- sapply(strsplit(as.character(coord$CDT)," in "), `[`, 1)
        coord[,"variable"] <- sapply(strsplit(as.character(coord$CDT)," in "), `[`, 2)
        coord <- separate(coord,"variable",c("Fraction","Zone","Cycle","Periods"),sep="X")
        # ggplot2
        palet_tree_grp <- c("#ad0635FF","#70a830FF","#98d048FF","#d84040FF","#60d0a0FF","#4080c0FF","#d7c11cFF","#40a078FF","#333b43FF")
        pcoa_plot_c <- ggplot(coord) + theme_unique_art() + 
          theme(legend.text=element_text(size=10),legend.position = "right") + 
          geom_point(aes(y= `PCoA 2`, x = `PCoA 1`,color = nbGroup, fill = Cycle, shape = Zone, size = Fraction), stroke =2) + 
          #geom_circle(aes(x0 = `PCoA 1`, y0 = `PCoA 3`, fill= Fraction, r=0.005), colour = "black")  + 
          scale_shape_manual(values = c(21, 24)) +
          scale_size_manual(values = c(4, 2)) +
          scale_color_manual(values = palet_tree_grp) + 
          scale_fill_manual(values = c("white","black")) + 
          labs(title="Abundance",color="Functional groups", shape = "Zone (shape)", size = "Fraction (size)", fill = "Cycle (fill)", x= paste0("Dim 1 [",Dim1Seq,"]"), y=paste0("Dim2 [",Dim2Seq,"]")) +
          facet_grid(Periods~.) +
          guides(fill = guide_legend(override.aes = list(shape = 21, size = 3)), shape = guide_legend(override.aes = list(size = 3))) +
          scale_y_continuous(position = 'left')
        print(pcoa_plot_c)
        LEGEND_c <- get_legend(pcoa_plot_c)
        pcoa_plot_c <- pcoa_plot_c + theme(legend.position = "none")
        #
      # Présence absence
        dataCDT[dataCDT>1] <- 1
        pcoa <- ape::pcoa(vegdist(t(dataCDT), method = "bray"))
        coord <- pcoa$vectors[,1:3] ; colnames(coord) <- c("PCoA 1","PCoA 2","PCoA 3")
        Dim1Seq <- paste0(round(pcoa$values[1,"Relative_eig"],3)*100,"%")
        Dim2Seq <- paste0(round(pcoa$values[2,"Relative_eig"],3)*100,"%")
        Dim3Seq <- paste0(round(pcoa$values[3,"Relative_eig"],3)*100,"%")
        coord <- as.data.frame(coord)
        coord$CDT <- row.names(coord)
        coord[,"nbGroup"] <- sapply(strsplit(as.character(coord$CDT)," in "), `[`, 1)
        coord[,"variable"] <- sapply(strsplit(as.character(coord$CDT)," in "), `[`, 2)
        coord <- separate(coord,"variable",c("Fraction","Zone","Cycle","Periods"),sep="X")
        # ggplot2
        palet_tree_grp <- c("#ad0635FF","#70a830FF","#98d048FF","#d84040FF","#60d0a0FF","#4080c0FF","#d7c11cFF","#40a078FF","#333b43FF")
        pcoa_plot_a <- ggplot(coord) + theme_unique_art() + 
          theme(legend.text=element_text(size=10),legend.position = "right") + 
          geom_point(aes(y= `PCoA 2`, x = `PCoA 1`,color = nbGroup, fill = Cycle, shape = Zone, size = Fraction), stroke =2) + 
          #geom_circle(aes(x0 = `PCoA 1`, y0 = `PCoA 3`, fill= Fraction, r=0.005), colour = "black")  + 
          scale_shape_manual(values = c(21, 24)) +
          scale_size_manual(values = c(4, 2)) +
          scale_color_manual(values = palet_tree_grp) + 
          scale_fill_manual(values = c("white","black")) + 
          labs(title="Presence/absence",color="Functional groups", shape = "Zone (shape)", size = "Fraction (size)", fill = "Cycle (fill)", x= paste0("Dim 1 [",Dim1Seq,"]"), y=paste0("Dim2 [",Dim2Seq,"]")) +
          facet_grid(Periods~.) +
          guides(fill = guide_legend(override.aes = list(shape = 21, size = 3)), shape = guide_legend(override.aes = list(size = 3))) +
          scale_y_continuous(position = 'left') + theme(legend.position = "none")
        print(pcoa_plot_a)
        #
        # Coplot
        svglite("DESeq2/PCoA/GRPandCDTS_PCoA_bray_Deseq_standart_Normalization.svg",width = 12.00,height = 9.00)
        b_plot <- plot_grid(pcoa_plot_a,pcoa_plot_c,LEGEND_c, ncol = 3, nrow = 1, rel_widths = c(3,3,2),rel_heights = c(3))
        print(b_plot)
        dev.off()
        #
# Taxonomy All fct° without Cycle ---------------------------------------------------------------
  Hist_DivisionXKO <- tableVinput_abundance_annots_tax_Funct %>% dplyr::select(all_of(MetaT_Summary$Echantillon),Kegg.Onthology,nbGroup,Supergroup,Division) %>% filter(is.na(nbGroup)==FALSE) %>% filter(is.na(Kegg.Onthology)==FALSE) %>% filter(!nbGroup %in% c("Lichen","Other multicellular Fungi", "Embryophyceae","Unaffiliated"))
  Hist_DivisionXKO <- Hist_DivisionXKO %>% mutate(Division = ifelse(is.na(Division) == TRUE, paste0('Unaffiliated ',Supergroup),Division))
  Hist_DivisionXKO[Hist_DivisionXKO=="Unaffiliated_NA"] <- NA
  color_tax_table_order <- unique(Hist_DivisionXKO %>% arrange(Supergroup,Division) %>% select(Division))
  #
  Hist_DivisionXKO <- Hist_DivisionXKO %>% dplyr::select(-"Kegg.Onthology",-"Supergroup") %>% group_by(nbGroup,Division) %>% summarise_all(sum)
  Hist_DivisionXKO <- reshape2::melt(Hist_DivisionXKO,id.vars = c("nbGroup","Division"))
  Hist_DivisionXKO <- merge(Hist_DivisionXKO, MetaT_Summary %>% select(Echantillon, Zone, Fraction, Periods),by.x = "variable",by.y = "Echantillon")
  Hist_DivisionXKO <- Hist_DivisionXKO %>% select(-variable) %>% group_by(nbGroup,Division,Zone,Fraction,Periods) %>% summarise_all(sum)
  #
  Hist_DivisionXKO <- as.data.frame(Hist_DivisionXKO)
  for (i in row.names(Hist_DivisionXKO)) {
    Fractionx <- Hist_DivisionXKO[i,"Fraction"]
    Zonex <- Hist_DivisionXKO[i,"Zone"]
    Periodsx <- Hist_DivisionXKO[i,"Periods"]
    nbGroupx <- Hist_DivisionXKO[i,"nbGroup"]
    Hist_DivisionXKO[i,"Sum"] <- sum(Hist_DivisionXKO %>% filter(Periods == Periodsx) %>% filter(Fraction == Fractionx) %>% filter(Zone == Zonex) %>% filter(nbGroup == nbGroupx) %>% dplyr::select(value))/nrow(Hist_DivisionXKO %>% filter(Periods == Periodsx) %>% filter(Fraction == Fractionx) %>% filter(Zone == Zonex) %>% filter(nbGroup == nbGroupx))
    }
  #Plot
  color_tax_table_order <- color_tax_table_order %>% filter(Division %in% Hist_DivisionXKO$Division)
  color_tax_table_order <- c(color_tax_table_order$Division)
  paletspe_i <- rep(c("#b71c1d","#f24336","#e57373",
                               "#512da8","#6437b2",
                               "#26579b","#3b88d1",
                               "#689f38",
                               "#f6c02d","#f9eb3b",#"#fdd835",
                               "#212121","#616161",
                               "#880e4f",
                               "#3c9cac","#48bcd4"),
                               4)
  Hist_DivisionXKO_fig <- ggplot() + 
    geom_bar(data = Hist_DivisionXKO,
             mapping = aes(y= value, x = Periods, fill = factor(Division, levels = color_tax_table_order), group = nbGroup),
             stat="identity",alpha=1) + 
    geom_bar(data = Hist_DivisionXKO %>% select(-Division) %>% group_by(nbGroup, Zone, Fraction,Periods) %>% summarise_all(sum), 
             mapping = aes(y= Sum, x = Periods, color = nbGroup), fill = NA,
             stat="identity") + 
    geom_text(data = Hist_DivisionXKO %>% select(-Division) %>% group_by(nbGroup, Zone, Fraction,Periods) %>% summarise_all(sum) %>% arrange(desc(nbGroup)),
              mapping = aes(y= Sum, x = Periods, label = nbGroup), 
              size = 2.5, position = position_stack(vjust = 0.5),
              fontface = "bold") +
    facet_grid(~Fraction+Zone,scales="free_x") + theme_unique_art() + 
    scale_fill_manual(values = paletspe_i, na.value = NA) + 
    scale_color_manual(values = rep("black",15)) + 
    labs(x="Periods",y="Transcripts count",fill="Divisions") + guides(fill=guide_legend(ncol=1)) +
    guides(color="none")
  print(Hist_DivisionXKO_fig)
  ggsave("GrpXDivisions_conditions_AllFunction*_abundance.svg", device = "svg", width = 16, height = 7)
  write.table(Hist_DivisionXKO,"FunctionalGroupContribXCDTS.tsv",row.names = FALSE,col.names = TRUE,quote = FALSE, sep ="\t")
  
  #
        