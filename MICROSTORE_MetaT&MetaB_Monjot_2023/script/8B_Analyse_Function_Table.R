# !/usr/bin/env Rscript
#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 31/05/2022
#
# Script Composition
# Set input and output -----------------------------------------------------------
# Detect R or Rstudio
  se <- Sys.getenv()
  if (length(se[grepl("rstudio",se,ignore.case=TRUE)]) > 0 ) { inputmode <- TRUE }
  if (length(se[grepl("rstudio",se,ignore.case=TRUE)]) == 0 ) { inputmode <- FALSE }
#
# Input argument if using Rstudio
  if (inputmode == TRUE) {
  input <- "Table_Supp_1.tsv"
  output <- "V4-unified-correct-paired-out-compo"
  region <- "V4"
  }
#
# Input argument if using R
  args = commandArgs(trailingOnly=TRUE)
  if ( inputmode == FALSE ) {
    input <- args[1]
    output <- args[2]
    region <- args[3]
  }
#
# Import package -----------------------------------------------------------
  pkg <- c("ggplot2","vegan","dplyr","tidyr","cowplot","FactoMineR","factoextra","reshape2","varhandle","ggrepel","ggpubr","ggsci","scales","hrbrthemes","GUniFrac","svglite","treemap", "VennDiagram","stringr","cluster","NbClust","tibble","treemapify","psych","gplots","ggExtra","paletteer","elementalist")
  lapply(pkg, require, character.only = TRUE)
#
# Theme unique Dark perso -------------------------------------------------------
  theme_unique_dark <- function (base_size = 12, base_family = "") {
    ret <- (theme_bw(base_size = base_size, base_family = base_family) +
              theme(text = element_text(colour = "black"),
                    title = element_text(color = "black", face = "bold"),
                    #axis.ticks = element_blank(),
                    line = element_line(color = "black"),
                    rect = element_rect(fill = "white", color = "black"),
                    axis.title = element_blank(),
                    #axis.text.y = element_blank(),
                    #axis.text.x = element_blank(),
                    axis.text.x = element_text(color = "black", size = 8, vjust = 1, hjust=1,angle = 45,face="bold"),
                    axis.line = element_line(color = "#969696", linetype = 1),
                    legend.background = element_rect(fill = NULL, color = NULL),
                    legend.position = "bottom",
                    legend.key = element_rect(fill = NA, color = NA, linetype = 0),
                    legend.text = element_text(size = 8),
                    legend.title = element_text(size = 8, face="bold"),
                    #strip.background = element_rect(fill=NA,colour=NA,size=NA,linetype = NULL),
                    strip.text = element_text(color="black",face="bold",vjust=.5,hjust=.5),
                    strip.text.y = element_text(angle=0),
                    panel.background = element_rect(fill = "white", color = NULL),
                    #panel.border = element_blank(),
                    panel.grid = element_line(color = "#252525"),
                    panel.grid.major = element_line(color = "white",linetype="blank"),
                    panel.grid.minor = element_line(color = "#252525",linetype="dotted"),
                    plot.background = element_rect(fill = "white", colour = "#252525", linetype = 0)))
    ret
  } 
  theme_unique_darkbis <- function (base_size = 12, base_family = "") {
    ret <- (theme_bw(base_size = base_size, base_family = base_family) +
              theme(text = element_text(colour = "black"),
                    title = element_text(color = "black", face = "bold", vjust = 0.5, hjust = 0.5,size = 10),
                    axis.ticks = element_blank(),
                    line = element_line(color = "black"),
                    rect = element_rect(fill = "white", color = "black"),
                    axis.title = element_text(color = "black", face = "bold"),
                    axis.text.y = element_blank(),
                    axis.text.x = element_blank(),
                    axis.line = element_blank(),
                    #axis.line = element_line(color = "#969696", linetype = 1),
                    legend.background = element_rect(fill = NULL, color = NULL),
                    legend.position = "none",
                    legend.key = element_rect(fill = NA, color = NA, linetype = 0),
                    legend.text = element_text(size = 8),
                    legend.title = element_text(size = 8, face="bold"),
                    strip.background = element_rect(fill=NA,colour=NA,size=NA,linetype = NULL),
                    strip.text = element_text(color="black",face="bold",vjust=.5,hjust=.5,size = 12),
                    panel.background = element_rect(fill = "white", color = NULL),
                    panel.border = element_blank(),
                    panel.grid = element_line(color = "#252525"),
                    panel.grid.major = element_line(color = "white"),
                    panel.grid.minor = element_line(color = "white"),
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
                    legend.title = element_text(size = 8, face="bold"),
                    legend.position = "none",
                    #panel.background = element_rect_round(radius = unit(0.5, "cm"),fill = "white", color = NULL,size = 1),
                    #strip.background = element_rect(fill=NA,colour=NA,size=NA,linetype = NULL),
                    strip.text = element_text(color="black",face="bold",vjust=.5,hjust=.5),
                    panel.background = element_rect(fill = "white", color = NULL),
                    panel.border = element_blank(),
                    panel.grid = element_line(color = "#252525"),
                    panel.grid.major = element_line(color = "grey86"),
                    panel.grid.minor = element_line(color = "grey86"),
                    plot.background = element_rect(fill = "white", colour = "white", linetype = 0)))
    ret
  } 
  theme_unique_darkcris <- function (base_size = 12, base_family = "") {
    ret <- (theme_bw(base_size = base_size, base_family = base_family) +
              theme(text = element_text(colour = "black"),
                    title = element_text(color = "black", face = "bold"),
                    #axis.ticks = element_blank(),
                    line = element_line(color = "black"),
                    rect = element_rect(fill = "white", color = "black"),
                    axis.title = element_blank(),
                    #axis.text.y = element_blank(),
                    #axis.text.x = element_blank(),
                    axis.text.x = element_text(color = "black", size = 8, vjust = 1, hjust=1,angle = 45,face="bold"),
                    axis.line = element_line(color = "#969696", linetype = 1),
                    legend.background = element_rect(fill = NULL, color = NULL),
                    #legend.position = "bottom",
                    legend.key = element_rect(fill = NA, color = NA, linetype = 0),
                    legend.text = element_text(size = 8),
                    legend.title = element_text(size = 8, face="bold"),
                    #strip.background = element_rect(fill=NA,colour=NA,size=NA,linetype = NULL),
                    strip.text = element_text(color="black",face="bold",vjust=.5,hjust=.5),
                    strip.text.y = element_text(angle=0),
                    panel.background = element_rect(fill = "white", color = NULL),
                    #panel.border = element_blank(),
                    panel.grid = element_line(color = "#252525"),
                    panel.grid.major = element_line(color = "white",linetype="blank"),
                    panel.grid.minor = element_line(color = "#252525",linetype="dotted"),
                    plot.background = element_rect(fill = "white", colour = "#252525", linetype = 0)))
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
  
#
# Set directory, create file result -----------------------------------------------------------
  if (inputmode == TRUE) {
    current <- dirname(rstudioapi::getSourceEditorContext()$path)
    setwd(current)
    output <- paste(current,"../dataDADA2/result",output, sep = "/")
  }
  if (inputmode == FALSE) {
    output <- paste("../dataDADA2/result",output, sep = "/")
  }
# Set working directory
  if (dir.exists(output) == FALSE) { dir.create(output, recursive = TRUE) }
  if (dir.exists(output) == TRUE) { setwd(output) }
#
# Import 8A.Rdata ---------------------------------------------------------
  load("Functional-Analyse/8A_Kmean_Clusterization.RData")
#
# Code-Table_Last.Group ---------------------------------------------------
    ## Code Table manually if necessary
    #data[,"nbGroup"][data[,"nbGroup"]=="Cluster 1"]<-"PARA"
    #data[,"nbGroup"][data[,"nbGroup"]=="Cluster 2"]<-"SAP"
    #data[,"nbGroup"][data[,"nbGroup"]=="Cluster 3"]<-"HCOV and SWAT"
    #data[,"nbGroup"][data[,"nbGroup"]=="Cluster 4"]<-"HET"
    #data[,"nbGroup"][data[,"nbGroup"]=="Cluster 5"]<-"PARA"
    #data[,"nbGroup"][data[,"nbGroup"]=="Cluster 6"]<-"MIXO"
    #data[,"nbGroup"][data[,"nbGroup"]=="Cluster 7"]<-"HET"
    #data[,"nbGroup"][data[,"nbGroup"]=="Cluster 8"]<-"FLAT"
    #data[,"nbGroup"][data[,"nbGroup"]=="Cluster 9"]<-"END"
    assoc_table <- read.table("../../../rawdata/cluster_name.tsv",header=TRUE,sep="\t")
    for (i in row.names(data)) { data[i,"nbGroup"] <- assoc_table %>% filter(nbGroup == data[i,"nbGroup"]) %>% select(name)}
    #
    ## Write table
    write.table(data,"Functional-Analyse/Group_FUNCTION/data_Last.Group.csv",row.names = FALSE,sep=";")
    #
# Input ASV_table ---------------------------------------------------------
    seq_mat_rare <- read.csv("ASV-Table/Seq-mat-rare.csv", row.names = 1, header= TRUE,sep="\t")
    asv_mat_rare <- read.csv("ASV-Table/ASV-mat-rare.csv", row.names = 1, header= TRUE,sep="\t")
    statRarefy<- read.csv("Stat-Analyse/StatRarefy_withoutDuplicat.csv", row.names = 1, header= TRUE,sep="\t")
    colnames(statRarefy)<-c("avRarefy-Sequence","apRarefy-Sequence","avRarefy-ASV","apRarefy-ASV")
    #
    # Prepare data inf --------------------------------------------------------
    infdataini <- read.table(file = "../../../rawdata/metadata_metaB", sep = "_", header = FALSE, col.names = c("Id","Condition","Period","Region","Replicat","empty"))
    infdataini$OSTA <- "OSTA"
    infdataini$Technologie  <- "Hiseq"
    infdataini$Amplicon <- paste(infdataini$Id,infdataini$OSTA, sep = "")
    infdataini <- infdataini %>% select(-"Id",-"OSTA")
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
      if (samples_df[i,"Period"] == "4") { samples_df[i,"Period"] <- "04"}
      if (samples_df[i,"Period"] == "6") { samples_df[i,"Period"] <- "06"}
      if (samples_df[i,"Period"] == "9") { samples_df[i,"Period"] <- "09"}
      if (samples_df[i,"Period"] == "11") { samples_df[i,"Period"] <- "11"}
    }
    #
    #
  # Select V4 or V9 ---------------------------------------------------------
    samples_df<-filter(samples_df, Region == region)
    #
# Create dataframe Condition --------------------------------------------
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
    # Sequence ----------------------------------------------------------------
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
  # Sequence ----------------------------------------------------------------
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
# Taxon-Group -------------------------------------------------------------
    matrix_data <- seq_mat_rare
    matrix_data$ASV_Id <- row.names(matrix_data)
    profile_tax_table <- as.data.frame(matrix_data %>% select(ASV_Id,Taxonomy))
    profile_tax_table <- as.data.frame(profile_tax_table[order(profile_tax_table$Taxonomy),])
    profile_tax_table <- separate(profile_tax_table, all_of(Taxonomy), c("Domain","Supergroup","Division","Class","Order","Family","Genus","Species"), sep =";")
    profile_tax_table[profile_tax_table=='NA'] <- NA
    profile_tax_table[,"Last"] <- NA
    ## Input exeption for "Last" Annotation
    exept <- c("environmental samples", "environvironmental samples", "unclassified", "unclassified","_X")
    ## If Last annotation doesn't exist with exeption, use prior annotation
    for (i in row.names(profile_tax_table)) {
      w <- length(profile_tax_table[i,][is.na(profile_tax_table[i,]) != TRUE])
      if (colnames(profile_tax_table[w])=="Species"){w=8}
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
    #
  # tax-Sequence ----------------------------------------------------------------
    data_seq_tax <- data_seq
    data_seq_tax<-merge(data_seq_tax,profile_tax_table,by="ASV_Id")
    data_seq_tax[is.na(data_seq_tax) == TRUE] <- "Not Affiliated"
    data$Last<-data$Row.names
    data_seq_tax <- left_join(x = data_seq_tax, y = data %>% select(Last,nbGroup), by = "Last")
    data_seq_tax[,"nbGroup"][data_seq_tax[,"Domain"]=="Bacteria"]<-"Bacteria"
    data_seq_tax[,"nbGroup"][data_seq_tax[,"Division"]=="Metazoa"]<-"Metazoa"
    data_seq_tax[,"nbGroup"][data_seq_tax[,"Class"]=="Embryophyceae"]<-"Embryophyceae"
    data_seq_tax[,"nbGroup"][data_seq_tax[,"Family"]=="Lecanoromycetes"]<-"Lichen"
    data_seq_tax[,"nbGroup"][data_seq_tax[,"Family"]=="Lichinomycetes"]<-"Lichen"
    data_seq_tax[,"nbGroup"][(data_seq_tax[,"Class"]=="Basidiomycota")&(is.na(data_seq_tax[,"nbGroup"])==TRUE)]<-"Other multicellular Fungi"
    data_seq_tax[,"nbGroup"][is.na(data_seq_tax[,"nbGroup"])==TRUE]<-"Unassigned"
    data_seq_tax$RangInterest <- data_seq_tax$nbGroup
    write.table(data_seq_tax, file = "Functional-Analyse/Table/ASV-Table-GrpFunct.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    #
  # tax-ASV ----------------------------------------------------------------
    data_asv_tax <- data_asv
    data_asv_tax<-merge(data_asv_tax,profile_tax_table,by="ASV_Id")
    data_asv_tax[is.na(data_asv_tax) == TRUE] <- "Not Affiliated"
    data$Last<-data$Row.names
    data_asv_tax <- left_join(x = data_asv_tax, y = data %>% select(Last,nbGroup), by = "Last")
    data_asv_tax[,"nbGroup"][data_asv_tax[,"Domain"]=="Bacteria"]<-"Bacteria"
    data_asv_tax[,"nbGroup"][data_asv_tax[,"Division"]=="Metazoa"]<-"Metazoa"
    data_asv_tax[,"nbGroup"][data_asv_tax[,"Class"]=="Embryophyceae"]<-"Embryophyceae"
    data_asv_tax[,"nbGroup"][data_asv_tax[,"Family"]=="Lecanoromycetes"]<-"Lichen"
    data_asv_tax[,"nbGroup"][data_asv_tax[,"Family"]=="Lichinomycetes"]<-"Lichen"
    data_asv_tax[,"nbGroup"][(data_asv_tax[,"Class"]=="Basidiomycota")&(is.na(data_asv_tax[,"nbGroup"])==TRUE)]<-"Other multicellular Fungi"
    data_asv_tax[,"nbGroup"][is.na(data_asv_tax[,"nbGroup"])==TRUE]<-"Unassigned"
    data_asv_tax$RangInterest <- data_asv_tax$nbGroup
    write.table(data_asv_tax, file = "Functional-Analyse/Table/ASV(asv)-Table-GrpFunct.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
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
    #
    paletspe <- rep(c("#e53a35","#f44336","#ef5350","#e57373","#ef9a9a","#ffcdd2",
                      "#7b1fa2","#8e24aa","#9c27b0",
                      "#303f9f",
                      "#0288d1","#049be5",
                      #"#01796b",
                      #"#689f38","#7cb342","#8bc34a",
                      "#689f38","#7cb342",
                      "#fbc02c","#fdd835","#ffeb3a","#ffee58","#fff176",
                      "#f57c00",
                      "#616161","#757575","#9e9e9e","#bdbdbd",
                      "#1876d2",
                      #"#455b64","#546e7a","#607d8b","#78909c","#90a4ae"),
                      #"#2e897b","#349688","#3ba69a","#4db6ac","#80cbc4"),
                      "#3897a7","#41acc1","#4cc6da","#80deea","#b2ebf2"),
                    #"#c2195a","#d2195c","#d4195d","#e61f61","#ea4079"),
                    4)
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
      ZoneXPeriods_tax <- merge(data_seq_tax %>% select(ASV_Id,RangInterest,Supergroup,Division), ZoneXPeriods, by = "ASV_Id")
      ZoneXPeriods_tax <- ZoneXPeriods_tax %>% select(-ASV_Id)
      for (i in rownames(ZoneXPeriods_tax)) { if (ZoneXPeriods_tax[i,"Division"]=="Not Affiliated" && ZoneXPeriods_tax[i,"Supergroup"]!="Not Affiliated") { ZoneXPeriods_tax[i,"Division"] <- paste("Unaffiliated",ZoneXPeriods_tax[i,"Supergroup"],sep=" ")}
        if (ZoneXPeriods_tax[i,"Division"]=="Alveolata_X" ) { ZoneXPeriods_tax[i,"Division"] <- paste("Unaffiliated",ZoneXPeriods_tax[i,"Supergroup"],sep=" ")}}
      ZoneXPeriods_tax <- ZoneXPeriods_tax %>% group_by(RangInterest,Supergroup,Division) %>% summarise_all(sum)

      ZoneXPeriods_tax_melt <-melt(ZoneXPeriods_tax,id.vars = c("RangInterest","Supergroup","Division"))
      ZoneXPeriods_tax_melt <- ZoneXPeriods_tax_melt %>% filter(!RangInterest %in% c("Lichen","Other multicellular Fungi", "Embryophyceae","Unaffiliated","Metazoa"))
      ZoneXPeriods_tax_melt <- separate(ZoneXPeriods_tax_melt,"variable",c("Zone","Periods"),sep="X")
      for (i in row.names(ZoneXPeriods_tax_melt)) {
        Periodsx <- ZoneXPeriods_tax_melt[i,"Periods"]
        Zonex <- ZoneXPeriods_tax_melt[i,"Zone"]
        RangInterestx <- ZoneXPeriods_tax_melt[i,"RangInterest"]
        Supergroupx <- ZoneXPeriods_tax_melt[i,"Supergroup"]
        Divisionx <- ZoneXPeriods_tax_melt[i,"Division"]
        ZoneXPeriods_tax_melt[i,"Sum"] <- sum(ZoneXPeriods_tax_melt %>% filter(Periods == Periodsx) %>% filter(Zone == Zonex) %>% filter(RangInterest == RangInterestx) %>% dplyr::select(value))/nrow(ZoneXPeriods_tax_melt %>% filter(Periods == Periodsx) %>% filter(Zone == Zonex) %>% filter(RangInterest == RangInterestx))
        ZoneXPeriods_tax_melt[i,"Relative"] <- ZoneXPeriods_tax_melt[i,"value"] *100 / sum(ZoneXPeriods_tax_melt %>% filter(Periods == Periodsx) %>% filter(Zone == Zonex) %>% dplyr::select(value))
        ZoneXPeriods_tax_melt[i,"Proportion"] <- ZoneXPeriods_tax_melt[i,"value"]*100/ZoneXPeriods_tax_melt[i,"Sum"]}
      ZoneXPeriods_tax_melt <- merge(ZoneXPeriods_tax_melt,color_tax %>% select(Label,palet), by.x = "Division", by.y = "Label")
      ZoneXPeriods_tax_melt <- ZoneXPeriods_tax_melt %>% arrange(Supergroup, Division, RangInterest)
      assign(paste("ZoneXPeriods_tax_melt",type,sep = "_"),ZoneXPeriods_tax_melt)
    }
    color_tax_table_order <- unique(ZoneXPeriods_tax_melt_Sequence$Division)
    color_palet_table_order <- unique(ZoneXPeriods_tax_melt_Sequence$palet)
    
    ### Sequence   
    y <- ZoneXPeriods_tax_melt_Sequence %>% group_by(palet, Division, Supergroup, RangInterest, Zone, Periods) %>% summarise_all(sum)
    y$mixte <- paste(y$RangInterest,y$Division, sep=" X ")
    x <- ZoneXPeriods_tax_melt_Sequence %>% select(-Division, -Supergroup, -palet) %>% group_by(RangInterest, Zone, Periods) %>% summarise_all(sum)
    svglite("Functional-Analyse/Hist_Function/ZoneXPeriods-Sequence-Fun-area.svg",width = 11.00,height = 9.00)
    Fun_seq_PeriodsXZone_fig <- ggplot(Rowv = NA, col = colMain, scale = "column") + 
      geom_area(data = y, mapping = aes(y= Relative, x = Periods, fill = factor(Division, level = color_tax_table_order), group = mixte),stat="identity",alpha=0.95, color = NA) +
      geom_area(data = x, mapping = aes(y= Relative, x = Periods, color = RangInterest, group = RangInterest),stat="identity", color = "black", fill = NA, size = 1) + 
      geom_text(data = x %>% arrange(desc(RangInterest)), mapping = aes(y= Relative, x = Periods, label = ifelse(Periods=="06", RangInterest,"")), 
                size = 4, position = position_stack(vjust = 0.5),
                fontface = "bold") +
      facet_grid(.~Zone,scales="free") + 
      theme_unique_art() +
      scale_fill_manual(values = rep(color_palet_table_order,8), na.value = NA) +
      labs(x="Periods",y="Sequence %",fill="Division") + 
      guides(fill=guide_legend(ncol=1))
    print(Fun_seq_PeriodsXZone_fig)
    dev.off()
    ### ASVs
    y <- ZoneXPeriods_tax_melt_ASV %>% group_by(palet, Division, Supergroup, RangInterest, Zone, Periods) %>% summarise_all(sum)
    y$mixte <- paste(y$RangInterest,y$Division, sep=" X ")
    x <- ZoneXPeriods_tax_melt_ASV %>% select(-Division, -Supergroup, -palet) %>% group_by(RangInterest, Zone, Periods) %>% summarise_all(sum)
    svglite("Functional-Analyse/Hist_Function/ZoneXPeriods-ASVs-Fun-area.svg",width = 11.00,height = 9.00)
    Fun_asv_PeriodsXZone_fig <- ggplot(Rowv = NA, col = colMain, scale = "column") + 
      geom_area(data = y, mapping = aes(y= Relative, x = Periods, fill = factor(Division, level = color_tax_table_order), group = mixte),stat="identity",alpha=0.95, color = NA) +
      geom_area(data = x, mapping = aes(y= Relative, x = Periods, color = RangInterest, group = RangInterest),stat="identity", color = "black", fill = NA, size = 1) + 
      geom_text(data = x %>% arrange(desc(RangInterest)), mapping = aes(y= Relative, x = Periods, label = ifelse(Periods=="06", RangInterest,"")), 
                size = 4, position = position_stack(vjust = 0.5),
                fontface = "bold") +
      facet_grid(.~Zone,scales="free") + 
      theme_unique_art() +
      scale_fill_manual(values = rep(color_palet_table_order,8), na.value = NA) +
      labs(x="Periods",y="ASVs %",fill="Division") + 
      guides(fill=guide_legend(ncol=1))
    print(Fun_asv_PeriodsXZone_fig)
    dev.off()
    ## Coplot
    legend <- get_legend(Fun_asv_PeriodsXZone_fig)
    Fun_asv_PeriodsXZone_fig <- Fun_asv_PeriodsXZone_fig + theme(legend.position = "none", axis.text.x = element_blank(),axis.title.x = element_blank()) + scale_x_discrete(position="top")
    Fun_seq_PeriodsXZone_fig <- Fun_seq_PeriodsXZone_fig + theme(legend.position = "none",strip.text = element_blank())

    svglite("Functional-Analyse/Hist_Function/ZoneXPeriods-Full.svg",width = 12.00,height = 12.00)
    ZoneXPeriods_Full_fig <- plot_grid(Fun_asv_PeriodsXZone_fig,NA,NA,NA,NA,legend,Fun_seq_PeriodsXZone_fig, NA,NA, ncol = 3, nrow = 3, rel_widths = c(4,1,1),rel_heights = c(2,0,2))
    print(ZoneXPeriods_Full_fig)
    dev.off()
#
# HIST  --------------------------------------------------------------------
    palette <- c("#D62768CC","#ADB6B6CC","#42B540CC","#357EBDCC","#1B1919CC","#ED0000CC","#FF7F0ECC","#925E9FCC","#D43F3ACC","#9632B8CC","#46B8DACC","#5CB85CCC","#EFC008CC")
  
  # Table_Hist ---------------------------------------------------------------------
    ## Cycle, Zone et Fraction
    for (vt in c("ASV","Sequence")) {
      if (vt=="ASV") { 
        input_data_table <- data_asv_tax }
      if (vt=="Sequence") { 
        input_data_table <- data_seq_tax }
      for (type in c("Total","Only")) {
        for (cdt in c("Cycle","Zone","Fraction")) {
          for (ctg in unique(input_data_table[,cdt])) {
            if (ctg != "Shared") {
              if (type == "Total") {
                interest_data<-paste("Total",ctg,vt,sep="_")
                interest<-paste0("Total",ctg)
                tablecourant <- input_data_table %>% select(ASV_Id,all_of(interest),RangInterest)}
              if (type == "Only") {
                interest_data<-paste("Only",ctg,vt,sep="_")
                interest<-paste0("Total",ctg)
                tablecourant <- input_data_table %>% filter(get(cdt) == ctg) %>% select(ASV_Id,all_of(interest),RangInterest)}
              row.names(tablecourant)<-tablecourant$ASV_Id ; tablecourant <- tablecourant %>% select(-ASV_Id)
              tablecourant <- tablecourant %>% group_by(RangInterest) %>% summarise_all(sum)
              tablecourant[,interest] <- tablecourant[,interest]*100/sum(tablecourant[,interest])
              tablecourant <- tablecourant %>% arrange(desc(RangInterest)) 
              tablecourant$labypos <- cumsum(tablecourant[,interest]) - 0.5*tablecourant[,interest]
              tablecourant$labypos <- as.matrix(tablecourant$labypos)
              tablecourant$label<-NA
              for (i in rownames(tablecourant)) { 
                tablecourant[i,"label"] <- paste0(round(tablecourant[i,interest],1), "%")
                if (tablecourant[i,"label"] == "0%") { tablecourant[i,"label"] <- NA}
                if (is.na(tablecourant[i,"label"]) == FALSE) { tablecourant[i,"label"] <- paste(tablecourant[i,"RangInterest"]," : ",tablecourant[i,"label"], sep = "")}
              }
              tablecourant[,cdt] <- ctg
              colnames(tablecourant)[2]<-"value"
              tablecourant$Sum <- 0
              tablecourant$Count <- 0
              if ((type == "Total")&&(vt=="ASV")) {
                for (i in tablecourant$RangInterest) { tablecourant$Count[which(tablecourant$RangInterest == i)] <- sum(input_data_table  %>% filter(RangInterest == i) %>% select(all_of(interest)))}
                tablecourant$Sum <- sum(tablecourant$Count)}
              if ((type == "Only")&&(vt=="ASV")) {
                for (i in tablecourant$RangInterest) { tablecourant$Count[which(tablecourant$RangInterest == i)] <- nrow(input_data_table %>% filter(get(cdt) == ctg) %>% filter(RangInterest == i))}
                tablecourant$Sum <- sum(tablecourant$Count)}
              if ((type == "Only")&&(vt=="Sequence")) {
                for (i in tablecourant$RangInterest) { tablecourant$Count[which(tablecourant$RangInterest == i)] <- sum(input_data_table %>% filter(get(cdt) == ctg) %>% filter(RangInterest == i) %>% select(all_of(interest)))}
                tablecourant$Sum <- sum(tablecourant$Count)}
              if ((type == "Total")&&(vt=="Sequence")) {
                for (i in tablecourant$RangInterest) { tablecourant$Count[which(tablecourant$RangInterest == i)] <- sum(input_data_table %>% filter(RangInterest == i) %>% select(all_of(interest)))}
                tablecourant$Sum <- sum(tablecourant$Count)}
              assign(interest_data,tablecourant)
            }
          }
          sel<-unique(input_data_table[,cdt]) ; sel<-sel[!sel %in% "Shared"]
          type_bind_data<-paste(type,cdt,vt,sep="_")
          type_bind_data1<-paste(type,sel[1],vt,sep="_")
          type_bind_data2<-paste(type,sel[2],vt,sep="_")
          assign(type_bind_data,as.data.frame(rbind(get(type_bind_data1),get(type_bind_data2))))
          if ((vt=="Sequence")&&(type=="Total")) { 
            assign(type_bind_data,get(type_bind_data) %>% mutate(percent= paste("(",round(Sum*100/(colSums(statRarefy %>% select(`apRarefy-Sequence`))/2),1)," %)",sep ="")))
          }
        }
      }
      ## Period
      for (Period in c("04","06","09","11")) {
        interest_data<-paste("Total",Period,vt,sep="_")
        interest<-paste0("Total",Period)
        tablecourant <- input_data_table %>% select(ASV_Id,all_of(interest),RangInterest)
        row.names(tablecourant)<-tablecourant$ASV_Id ; tablecourant <- tablecourant %>% select(-ASV_Id)
        tablecourant <- tablecourant %>% group_by(RangInterest) %>% summarise_all(sum)
        tablecourant[,interest] <- tablecourant[,interest]*100/sum(tablecourant[,interest])
        tablecourant <- tablecourant %>% arrange(desc(RangInterest)) 
        tablecourant$labypos <- cumsum(tablecourant[,interest]) - 0.5*tablecourant[,interest]
        tablecourant$labypos <- as.matrix(tablecourant$labypos)
        tablecourant$label<-NA
        for (i in rownames(tablecourant)) { 
          tablecourant[i,"label"] <- paste0(round(tablecourant[i,interest],1), "%")
          if (tablecourant[i,"label"] == "0%") { tablecourant[i,"label"] <- NA}
          if (is.na(tablecourant[i,"label"]) == FALSE) { tablecourant[i,"label"] <- paste(tablecourant[i,"RangInterest"]," : ",tablecourant[i,"label"], sep = "")}
        }
        tablecourant[,"Period"] <- Period
        colnames(tablecourant)[2]<-"value"
        tablecourant$Sum <- 0
        tablecourant$Count <- 0
        for (i in tablecourant$RangInterest) { tablecourant$Count[which(tablecourant$RangInterest == i)] <- sum(input_data_table  %>% filter(RangInterest == i) %>% select(all_of(interest)))}
        tablecourant$Sum <- sum(tablecourant$Count)
        assign(interest_data,tablecourant)
      }
      sel<-c("04","06","09","11")
      type_bind_data<-paste("Total","Period",vt,sep="_")
      type_bind_data1<-paste("Total",sel[1],vt,sep="_")
      type_bind_data2<-paste("Total",sel[2],vt,sep="_")
      type_bind_data3<-paste("Total",sel[3],vt,sep="_")
      type_bind_data4<-paste("Total",sel[4],vt,sep="_")
      assign(type_bind_data,as.data.frame(rbind(get(type_bind_data1),get(type_bind_data2),get(type_bind_data3),get(type_bind_data4))))
      if (vt=="Sequence") { 
        assign(type_bind_data,get(type_bind_data) %>% mutate(percent= paste("(",round(Sum*100/(colSums(statRarefy %>% select(`apRarefy-Sequence`))/2),1)," %)",sep ="")))
      }
    }
    ## Only_Sequence_90
    input_data_table=data_seq_tax
    for (cdt in c("Cycle90","Zone90","Fraction90")){
      for (ctg in unique(input_data_table[,cdt])) {
        if (ctg != "Shared") {
          interest_data<-paste("Only",ctg,"Sequence90",sep="_")
          interest<-paste0("Total",ctg)
          tablecourant <- input_data_table %>% filter(get(cdt) == ctg) %>% select(ASV_Id,all_of(interest),RangInterest)
          row.names(tablecourant)<-tablecourant$ASV_Id ; tablecourant <- tablecourant %>% select(-ASV_Id)
          tablecourant <- tablecourant %>% group_by(RangInterest) %>% summarise_all(sum)
          tablecourant[,interest] <- tablecourant[,interest]*100/sum(tablecourant[,interest])
          tablecourant <- tablecourant %>% arrange(desc(RangInterest)) 
          tablecourant$labypos <- cumsum(tablecourant[,interest]) - 0.5*tablecourant[,interest]
          tablecourant$labypos <- as.matrix(tablecourant$labypos)
          tablecourant$label<-NA
          for (i in rownames(tablecourant)) {
            tablecourant[i,"label"] <- paste0(round(tablecourant[i,interest],1), "%")
            if (tablecourant[i,"label"] == "0%") { tablecourant[i,"label"] <- NA}
            if (is.na(tablecourant[i,"label"]) == FALSE) { tablecourant[i,"label"] <- paste(tablecourant[i,"RangInterest"]," : ",tablecourant[i,"label"], sep = "")}
          }
          tablecourant[,cdt] <- ctg
          colnames(tablecourant)[2]<-"value"
          tablecourant$Sum <- 0
          tablecourant$Count <- 0
          for (i in tablecourant$RangInterest) { tablecourant$Count[which(tablecourant$RangInterest == i)] <- sum(input_data_table %>% filter(get(cdt) == ctg) %>% filter(RangInterest == i) %>% select(all_of(interest)))}
          tablecourant$Sum <- sum(tablecourant$Count)}
        assign(interest_data,tablecourant)
      }
      sel<-unique(input_data_table[,cdt]) ; sel<-sel[!sel %in% "Shared"]
      type_bind_data<-paste("Only",cdt,"Sequence90",sep="_")
      type_bind_data1<-paste("Only",sel[1],"Sequence90",sep="_")
      type_bind_data2<-paste("Only",sel[2],"Sequence90",sep="_")
      assign(type_bind_data,as.data.frame(rbind(get(type_bind_data1),get(type_bind_data2))))
    }
    #
  # Figure_Hist ------------------------------------------------------------------
    # Total_ASV -------------------------------------------------------------------
    ## Cycle
    Total_Cycle_ASV_fig <- ggplot(Total_Cycle_ASV, mapping = aes(y= value, x = Cycle, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + scale_fill_paletteer_d("ggthemes::manyeys") + 
      geom_label(aes(y = 106,label = paste(Sum," ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") + labs(x="Cycles",y="ASVs (%)") + theme(legend.position = "none")
    print(Total_Cycle_ASV_fig)
    ## Zone
    Total_Zone_ASV_fig <- ggplot(Total_Zone_ASV, mapping = aes(y= value, x = Zone, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + scale_fill_paletteer_d("ggthemes::manyeys") + 
      geom_label(aes(y = 106,label = paste(Sum," ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") + labs(x="Zones",y="ASVs (%)") + theme(legend.position = "none")
    print(Total_Zone_ASV_fig)
    ## Fraction
    Total_Fraction_ASV_fig <- ggplot(Total_Fraction_ASV, mapping = aes(y= value, x = Fraction, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + scale_fill_paletteer_d("ggthemes::manyeys")+ 
      geom_label(aes(y = 106,label = paste(Sum," ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") + labs(x="Fractions",y="ASVs (%)") + theme(legend.position = "none")
    print(Total_Fraction_ASV_fig) 
    ## Period
    Total_Period_ASV_fig <- ggplot(Total_Period_ASV, mapping = aes(y= value, x = Period, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + scale_fill_paletteer_d("ggthemes::manyeys")+ 
      geom_label(aes(y = 106,label = paste(Sum," ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") + labs(x="Periods",y="ASVs (%)")
    print(Total_Period_ASV_fig) 
    ## Coplote
    svglite("Functional-Analyse/Hist_Function/ASV-Total.svg",width = 12.00,height = 6.00)
    Total_ASV_fig <- plot_grid(Total_Cycle_ASV_fig,Total_Zone_ASV_fig,Total_Fraction_ASV_fig,Total_Period_ASV_fig, ncol = 4, nrow = 1, rel_widths = c(3,3,3,6),rel_heights = c(3))
    print(Total_ASV_fig)
    dev.off()
    #
    # Only_ASV -------------------------------------------------------------------
    ## Cycle
    Only_Cycle_ASV_fig <- ggplot(Only_Cycle_ASV, mapping = aes(y= value, x = Cycle, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + scale_fill_paletteer_d("ggthemes::manyeys") + 
      geom_label(aes(y = 106,label = paste(Sum," ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") + labs(x="Cycles",y="ASVs (%)")
    print(Only_Cycle_ASV_fig)
    ## Zone
    Only_Zone_ASV_fig <- ggplot(Only_Zone_ASV, mapping = aes(y= value, x = Zone, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + scale_fill_paletteer_d("ggthemes::manyeys") + 
      geom_label(aes(y = 106,label = paste(Sum," ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") + labs(x="Zones",y="ASVs (%)")
    print(Only_Zone_ASV_fig)
    ## Fraction
    Only_Fraction_ASV_fig <- ggplot(Only_Fraction_ASV, mapping = aes(y= value, x = Fraction, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + scale_fill_paletteer_d("ggthemes::manyeys")+ 
      geom_label(aes(y = 106,label = paste(Sum," ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") + labs(x="Fractions",y="ASVs (%)")
    print(Only_Fraction_ASV_fig) 
    ## Coplote
    svglite("Functional-Analyse/Hist_Function/ASV-Only.svg",width = 12.00,height = 6.00)
    Only_ASV_fig <- plot_grid(Only_Cycle_ASV_fig,Only_Zone_ASV_fig,Only_Fraction_ASV_fig, ncol = 3, nrow = 1, rel_widths = c(3,3,3),rel_heights = c(3))
    print(Only_ASV_fig)
    dev.off()
    #
    # Total_Sequence -------------------------------------------------------------------
    ## Cycle
    Total_Cycle_Sequence_fig <- ggplot(Total_Cycle_Sequence, mapping = aes(y= value, x = Cycle, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
      scale_fill_paletteer_d("ggthemes::manyeys") + geom_label(aes(y = 108,label = paste(Sum,"sequences",percent,sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
      labs(x="Cycles",y="Sequences (%)") + theme(legend.position = "none")
    print(Total_Cycle_Sequence_fig)
    ## Zone
    Total_Zone_Sequence_fig <- ggplot(Total_Zone_Sequence, mapping = aes(y= value, x = Zone, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
      scale_fill_paletteer_d("ggthemes::manyeys") + geom_label(aes(y = 108,label = paste(Sum,"sequences",percent,sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
      labs(x="Zones",y="Sequences (%)") + theme(legend.position = "none")
    print(Total_Zone_Sequence_fig)
    ## Fraction
    Total_Fraction_Sequence_fig <- ggplot(Total_Fraction_Sequence, mapping = aes(y= value, x = Fraction, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
      scale_fill_paletteer_d("ggthemes::manyeys") + geom_label(aes(y = 108,label = paste(Sum,"sequences",percent,sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
      labs(x="Fractions",y="Sequences (%)") + theme(legend.position = "none")
    print(Total_Fraction_Sequence_fig)
    ## Period
    Total_Period_Sequence_fig <- ggplot(Total_Period_Sequence, mapping = aes(y= value, x = Period, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
      scale_fill_paletteer_d("ggthemes::manyeys") + geom_label(aes(y = 108,label = paste(Sum,"sequences",percent,sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
      labs(x="Periods",y="Sequences (%)")
    print(Total_Period_Sequence_fig)
    ## Coplote
    svglite("Functional-Analyse/Hist_Function/Sequence-Total.svg",width = 12.00,height = 6.00)
    Total_Sequence_fig <- plot_grid(Total_Cycle_Sequence_fig,Total_Zone_Sequence_fig,Total_Fraction_Sequence_fig,Total_Period_Sequence_fig, ncol = 4, nrow = 1, rel_widths = c(3,3,3,6),rel_heights = c(3))
    print(Total_Sequence_fig)
    dev.off()
    #
    # Only_Sequence -------------------------------------------------------------------
      # Sequence100 -------------------------------------------------------------------
    ## Cycle
    Only_Cycle_Sequence_fig <- ggplot(Only_Cycle_Sequence, mapping = aes(y= value, x = Cycle, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
      scale_fill_paletteer_d("ggthemes::manyeys") + geom_label(aes(y = 106,label = paste(Sum,"sequences",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
      labs(x="Cycles",y="Sequences (%)")
    print(Only_Cycle_Sequence_fig)
    ## Zone
    Only_Zone_Sequence_fig <- ggplot(Only_Zone_Sequence, mapping = aes(y= value, x = Zone, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
      scale_fill_paletteer_d("ggthemes::manyeys") + geom_label(aes(y = 106,label = paste(Sum,"sequences",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
      labs(x="Zones",y="Sequences (%)")
    print(Only_Zone_Sequence_fig)
    ## Fraction
    Only_Fraction_Sequence_fig <- ggplot(Only_Fraction_Sequence, mapping = aes(y= value, x = Fraction, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
      scale_fill_paletteer_d("ggthemes::manyeys") + geom_label(aes(y = 106,label = paste(Sum,"sequences",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
      labs(x="Fractions",y="Sequences (%)")
    print(Only_Fraction_Sequence_fig)
    ## Coplote
    svglite("Functional-Analyse/Hist_Function/Sequence100-Only.svg",width = 12.00,height = 6.00)
    Only_Sequence100_fig <- plot_grid(Only_Cycle_Sequence_fig,Only_Zone_Sequence_fig,Only_Fraction_Sequence_fig, ncol = 3, nrow = 1, rel_widths = c(3,3,3),rel_heights = c(3))
    print(Only_Sequence100_fig)
    dev.off()
    #
      # Sequence90 -------------------------------------------------------------------
    ## Cycle
    Only_Cycle90_Sequence90_fig <- ggplot(Only_Cycle90_Sequence90, mapping = aes(y= value, x = Cycle90, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
      scale_fill_paletteer_d("ggthemes::manyeys") + geom_label(aes(y = 106,label = paste(Sum,"sequences",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
      labs(x="Cycles",y="Sequences (%)")
    print(Only_Cycle90_Sequence90_fig)
    ## Zone
    Only_Zone90_Sequence90_fig <- ggplot(Only_Zone90_Sequence90, mapping = aes(y= value, x = Zone90, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
      scale_fill_paletteer_d("ggthemes::manyeys") + geom_label(aes(y = 106,label = paste(Sum,"sequences",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
      labs(x="Zones",y="Sequences (%)")
    print(Only_Zone90_Sequence90_fig)
    ## Fraction
    Only_Fraction90_Sequence90_fig <- ggplot(Only_Fraction90_Sequence90, mapping = aes(y= value, x = Fraction90, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
      scale_fill_paletteer_d("ggthemes::manyeys") + geom_label(aes(y = 106,label = paste(Sum,"sequences",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
      labs(x="Fractions",y="Sequences (%)")
    print(Only_Fraction90_Sequence90_fig)
    ## Coplote
    svglite("Functional-Analyse/Hist_Function/Sequence90-Only.svg",width = 12.00,height = 6.00)
    Only_Sequence90_fig <- plot_grid(Only_Cycle90_Sequence90_fig,Only_Zone90_Sequence90_fig,Only_Fraction90_Sequence90_fig, ncol = 3, nrow = 1, rel_widths = c(3,3,3),rel_heights = c(3))
    print(Only_Sequence90_fig)
    dev.off()
    #
# Tree --------------------------------------------------------------------
    palet_tree_grp <- c("#607848FF","#ad0635FF","#d84040FF","#4080c0FF","#60d0a0FF","#70a830FF","#98d048FF","#d7c11cFF","#40a078FF","#333b43FF","#2c7890FF","#b4d2d2")
    for (type in c("Only","Total")) {
      for (tv in c("ASV","Sequence")) {
        for (input_tree in c("Cycle","Fraction","Zone")) {
          table_tree <- get(paste(type,input_tree,tv,sep="_"))
          table_tree$labelx<-NA
          for (i in row.names(table_tree)) {
            if (is.na(table_tree[i,"label"])==TRUE) {table_tree[i,"labelx"]<-NA}
            else {
              if (str_split(table_tree[i,"label"]," ")[[1]][1]=="Group") {table_tree[i,"labelx"] <- str_split(table_tree[i,"label"],": ")[[1]][3]}
              if (str_split(table_tree[i,"label"]," ")[[1]][1]!="Group") {table_tree[i,"labelx"] <- str_split(table_tree[i,"label"],": ")[[1]][2]}
            }}
          q <-ggplot(table_tree, aes(area = value, fill = RangInterest, label = "")) +
            geom_treemap(alpha=0.80,color="white",size=2) +
            scale_fill_manual(values = palet_tree_grp) + 
            facet_grid(.~get(input_tree)) +
            theme_unique_darktris() +
            theme(legend.position="bottom") +
            geom_treemap_text(fontface = "bold", colour = "black", place = "centre",grow = FALSE,size=12)
          svglite(paste(paste0("Functional-Analyse/Tree_Function/",type),input_tree,tv,"Tree.svg",sep="_"),width = 12.00,height = 6.00)
          print(q)
          dev.off()
          ## with label
          r <-ggplot(table_tree, aes(area = value, fill = RangInterest, label = paste0(round(value,1),"%"))) +
            geom_treemap(alpha=0.80,color="white",size=2) +
            scale_fill_manual(values = palet_tree_grp) + 
            facet_grid(.~get(input_tree)) +
            theme_unique_darktris() +
            theme(legend.position="bottom") +
            geom_treemap_text(fontface = "bold", colour = "black", place = "centre",grow = FALSE,size=12)
          svglite(paste(paste0("Functional-Analyse/Tree_Function/",type),input_tree,tv,"Tree_Lab.svg",sep="_"),width = 12.00,height = 6.00)
          print(r)
          dev.off()
        }
      } 
    }
    ## Only_Sequence90
    for (input_tree in c("Cycle90","Fraction90","Zone90")) {
      table_tree <- get(paste("Only",input_tree,"Sequence90",sep="_"))
      table_tree$labelx<-NA
      for (i in row.names(table_tree)) {
        if (is.na(table_tree[i,"label"])==TRUE) {table_tree[i,"labelx"]<-NA}
        else {
          if (str_split(table_tree[i,"label"]," ")[[1]][1]=="Group") {table_tree[i,"labelx"] <- str_split(table_tree[i,"label"],": ")[[1]][3]}
          if (str_split(table_tree[i,"label"]," ")[[1]][1]!="Group") {table_tree[i,"labelx"] <- str_split(table_tree[i,"label"],": ")[[1]][2]}
        }}
      q <-ggplot(table_tree, aes(area = value, fill = RangInterest, label = "")) +
        geom_treemap(alpha=0.80,color="white",size=2) +
        scale_fill_manual(values = palet_tree_grp) + 
        facet_grid(.~get(input_tree)) +
        theme_unique_darktris() +
        theme(legend.position="bottom") +
        geom_treemap_text(fontface = "bold", colour = "black", place = "centre",grow = FALSE,size=12)
      svglite(paste(paste0("Functional-Analyse/Tree_Function/","Only"),input_tree,"Sequence90","Tree.svg",sep="_"),width = 12.00,height = 6.00)
      print(q)
      dev.off()
      ## with label
      r <-ggplot(table_tree, aes(area = value, fill = RangInterest, label = paste0(round(value,1),"%"))) +
        geom_treemap(alpha=0.80,color="white",size=2) +
        scale_fill_manual(values = palet_tree_grp) + 
        facet_grid(.~get(input_tree)) +
        theme_unique_darktris() +
        theme(legend.position="bottom") +
        geom_treemap_text(fontface = "bold", colour = "black", place = "centre",grow = FALSE,size=12)
      svglite(paste(paste0("Functional-Analyse/Tree_Function/","Only"),input_tree,"Sequence90","Tree_Lab.svg",sep="_"),width = 12.00,height = 6.00)
      print(r)
      dev.off()
    }
    #
# Polar ------------------------------------------------------
    #palet_tree_grp_x <- c("#ad0635FF","#d84040FF","#4080c0FF","#60d0a0FF","#70a830FF","#98d048FF","#d7c11cFF","#b4d2d2")
    palet_tree_grp_x <- c("#ad0635FF","#70a830FF","#98d048FF","#d84040FF","#60d0a0FF","#4080c0FF","#d7c11cFF","#b4d2d2")
    
  # Sequence ----------------------------------------------------------------
    ## Total column
    Polar_seq <- data_seq_tax
    Polar_seq <- Polar_seq %>% filter(nbGroup %in% c("Unassigned",unique(data$nbGroup)))
    Polar_seq$TotalAnnee <- 0
    for (i in rownames(Polar_seq)) { Polar_seq[i,"TotalAnnee"] <- Polar_seq[i,"Total04"] + Polar_seq[i,"Total06"] + Polar_seq[i,"Total09"] + Polar_seq[i,"Total11"]}
    Polar_seq <- Polar_seq %>% select(TotalAnnee,RangInterest)
    ## Total Figure
    Polar_seq <- Polar_seq %>% group_by(RangInterest) %>% summarise_all(sum)
    for (i in rownames(Polar_seq)) { Polar_seq[i,"Total"] <- (Polar_seq[i,"TotalAnnee"] * 100) / sum(Polar_seq$TotalAnnee)}
    ## Label
    Polar_seq <- Polar_seq %>%
      arrange(desc(RangInterest)) %>%
      mutate(lab.ypos = cumsum(Total) - 0.5*Total)
    Polar_seq$label <- paste(round(Polar_seq$Total,1), "%", sep = "")
    for (i in rownames(Polar_seq)) {
      if (Polar_seq[i,"label"] == "0%") { Polar_seq[i,"label"] <- NA}}
    for (i in rownames(Polar_seq)) {
      if (is.na(Polar_seq[i,"label"]) == FALSE) { Polar_seq[i,"label"] <- paste(Polar_seq[i,"RangInterest"]," : ",Polar_seq[i,"label"], sep = "")}}
    ## Figure
    svglite("Functional-Analyse/Composition_Function/Polar-Total-seq.svg",width = 4.50,height = 4.50)
    axim <- ggplot(Polar_seq, mapping = aes(y= Total, x = 2, fill = RangInterest), Rowv = NA, col = colMain, scale = "column") +
      geom_bar(stat="identity", color = "black", width = 1,alpha=0.8) + coord_polar("y") + 
      #geom_label_repel(aes(y = lab.ypos,label = label),color = "black",size = 3,segment.color = "black",show.legend = FALSE, nudge_x = 0.5) +
      geom_text_repel(data = Polar_seq, aes(x=2,y = lab.ypos,label = label),fontface= "bold",color = "black",size = 3,segment.color = "black",show.legend = FALSE, nudge_x = 0.9, force_pull = 0) +
      scale_y_continuous(limits=c(0,sum(Polar_seq %>% select(Total)+0.01)),position="right")+
      xlim(1,3) +
      theme_unique_darkbis() + 
      #facet_wrap( ~ variable , nrow = 4) +
      labs(y = "Sequences",x="") + scale_fill_manual(values=palet_tree_grp_x) #scale_fill_paletteer_d("ggthemes::manyeys")
    print(axim)
    dev.off()
    ## table
    write.table(Polar_seq,file="Functional-Analyse/Composition_Function/Table-Total-seq.csv",row.names = F, quote = F, sep = "\t")
    #
  # Sequence-Grp ----------------------------------------------------------------
    k=0
    for (l in unique(data$nbGroup)) { print(l)
      k<-k+1
      ## Total column
      Polar_seq <- data_seq_tax
      Polar_seq <- Polar_seq %>% filter(nbGroup==l)
      Polar_seq$TotalAnnee <- 0
      for (i in rownames(Polar_seq)) { Polar_seq[i,"TotalAnnee"] <- Polar_seq[i,"Total04"] + Polar_seq[i,"Total06"] + Polar_seq[i,"Total09"] + Polar_seq[i,"Total11"]}
      Polar_seq <- Polar_seq %>% select(TotalAnnee,Division)
      ## Total Figure
      Polar_seq <- Polar_seq %>% group_by(Division) %>% summarise_all(sum)
      for (i in rownames(Polar_seq)) { Polar_seq[i,"Total"] <- (Polar_seq[i,"TotalAnnee"] * 100) / sum(Polar_seq$TotalAnnee)}
      ## Label
      Polar_seq <- Polar_seq %>%
        arrange(desc(Division)) %>%
        mutate(lab.ypos = cumsum(Total) - 0.5*Total)
      Polar_seq$label <- paste(round(Polar_seq$Total,1), "%", sep = "")
      for (i in rownames(Polar_seq)) {
        if (as.numeric(str_split(Polar_seq[i,"label"],"%")[[1]][1]) <= 0) { Polar_seq[i,"label"] <- NA}}
      for (i in rownames(Polar_seq)) {
        if (is.na(Polar_seq[i,"label"]) == FALSE) { Polar_seq[i,"label"] <- paste(Polar_seq[i,"Division"]," : ",Polar_seq[i,"label"], sep = "")}}
      Polar_seq$Group<-l
      if (k==1) {Polar_seq_GRP<-Polar_seq}
      if (k>1) {Polar_seq_GRP <- rbind(Polar_seq_GRP,Polar_seq)}
    }
    paletxi <- rev(c("#c04848CC","#f27435CC","#707470CC","#948c75CC","#424254CC","#e04644CC","#706acfCC","#6a4a3cCC","#2a044aCC","#e84a5fCC","#f02475CC","#e8bf56CC","#1b676bCC","#3b2d38CC","#cc2a41CC","#f02475CC","#cfbe27CC","#b8af03CC","#d27f48CC","#2cb0e0CC","#3fb8afCC","#3B3B3B99","#519548CC","#c2d1fcCC","#555152CC","#8fbe00CC","#64908aCC","#0e2430CC","#bcbdacCC"))
    ## Figure
    svglite("Functional-Analyse/Composition_Function/Polar-Total-seq-grp.svg",width = 10,height = 10)
    ax <- ggplot(Polar_seq_GRP, mapping = aes(y= Total, x = 2, fill = Division), Rowv = NA, col = colMain, scale = "column") +
      geom_bar(stat="identity", color = "white", width = 1) + coord_polar("y") + 
      geom_label_repel(aes(y = lab.ypos,label = label),color = "white",size = 3,segment.color = "black",show.legend = FALSE, nudge_x = 0.5,nudge_y = 0.5,max.overlaps = 20) +
      scale_y_continuous(limits=c(0,sum(Polar_seq_GRP %>% select(Total)+0.01)/length(unique(data$nbGroup))))+
      xlim(1,2.5) +
      theme_unique_darkbis() + 
      facet_wrap( ~ Group , nrow = 3) +
      labs(y = "Sequences",x="") + scale_fill_manual(values = paletxi) #+ scale_fill_paletteer_d("ggthemes::manyeys")
    print(ax)
    dev.off()
    ## Table
    write.table(Polar_seq_GRP,file="Functional-Analyse/Composition_Function/Table-Total-seq-grp.csv",row.names = F, quote = F, sep = "\t")
    ## Table Class
    seq_table_class_grp <- data_seq_tax %>% select(Total04,Total06,Total09,Total11,TotalMonimolimnion,TotalMixolimnion,Domain,Supergroup,Division,Class,Order,RangInterest)
    seq_table_class_grp <- seq_table_class_grp %>% group_by(Domain,Supergroup,Division,Class,Order,RangInterest) %>% summarise_all(sum) %>% arrange(RangInterest)
    write.table(seq_table_class_grp, file="Functional-Analyse/Composition_Function/Table-seq-grp.csv",row.names = F, quote = F,sep="\t")
#
  # ASV ----------------------------------------------------------------
    ## Total column
    Polar_asv <- data_asv_tax
    Polar_asv <- Polar_asv %>% filter(nbGroup %in% c("Unassigned",unique(data$nbGroup)))
    Polar_asv$TotalAnnee <- 0
    for (i in rownames(Polar_asv)) { if (Polar_asv[i,"Total04"] + Polar_asv[i,"Total06"] + Polar_asv[i,"Total09"] + Polar_asv[i,"Total11"] > 0) {Polar_asv[i,"TotalAnnee"] <- 1}}
    Polar_asv <- Polar_asv %>% select(TotalAnnee,RangInterest)
    ## Total Figure
    Polar_asv <- Polar_asv %>% group_by(RangInterest) %>% summarise_all(sum)
    for (i in rownames(Polar_asv)) { Polar_asv[i,"Total"] <- (Polar_asv[i,"TotalAnnee"] * 100) / sum(Polar_asv$TotalAnnee)}
    ## Label
    Polar_asv <- Polar_asv %>%
      arrange(desc(RangInterest)) %>%
      mutate(lab.ypos = cumsum(Total) - 0.5*Total)
    Polar_asv$label <- paste(round(Polar_asv$Total,1), "%", sep = "")
    for (i in rownames(Polar_asv)) {
      if (Polar_asv[i,"label"] == "0%") { Polar_asv[i,"label"] <- NA}}
    for (i in rownames(Polar_asv)) {
      if (is.na(Polar_asv[i,"label"]) == FALSE) { Polar_asv[i,"label"] <- paste(Polar_asv[i,"RangInterest"]," : ",Polar_asv[i,"label"], sep = "")}}
    ## Figure
    svglite("Functional-Analyse/Composition_Function/Polar-Total-asv.svg",width = 4.50,height = 4.50)
    bxim <- ggplot(Polar_asv, mapping = aes(y= Total, x = 2, fill = RangInterest), Rowv = NA, col = colMain, scale = "column") +
      geom_bar(stat="identity", color = "black", width = 1,alpha=0.8) + coord_polar("y") + 
      #geom_label_repel(aes(y = lab.ypos,label = label),color = "white",size = 3,segment.color = "black",show.legend = FALSE, nudge_x = 0.5) +
      geom_text_repel(data = Polar_asv, aes(x=2,y = lab.ypos,label = label),fontface= "bold",color = "black",size = 3,segment.color = "black",show.legend = FALSE, nudge_x = 0.9, force_pull = 0) +
      scale_y_continuous(limits=c(0,sum(Polar_asv %>% select(Total)+0.01)),position="right")+
      xlim(1,3) +
      theme_unique_darkbis() + 
      #facet_wrap( ~ variable , nrow = 4) +
      labs(y = "ASVs",x="") + scale_fill_manual(values=palet_tree_grp_x) #scale_fill_paletteer_d("ggthemes::manyeys")
    print(bxim)
    dev.off()
    ## Table
    write.table(Polar_asv,file="Functional-Analyse/Composition_Function/Table-Total-asv.csv",row.names = F, quote = F, sep = "\t")
    #
  # ASV-Grp ----------------------------------------------------------------
    k=0
    for (l in unique(data$nbGroup)) { print(l)
      k<-k+1
      ## Total column
      Polar_asv <- data_asv_tax
      Polar_asv <- Polar_asv %>% filter(nbGroup==l)
      Polar_asv$TotalAnnee <- 0
      for (i in rownames(Polar_asv)) { if (Polar_asv[i,"Total04"] + Polar_asv[i,"Total06"] + Polar_asv[i,"Total09"] + Polar_asv[i,"Total11"] > 0) {Polar_asv[i,"TotalAnnee"] <- 1}}
      Polar_asv <- Polar_asv %>% select(TotalAnnee,Division)
      ## Total Figure
      Polar_asv <- Polar_asv %>% group_by(Division) %>% summarise_all(sum)
      for (i in rownames(Polar_asv)) { Polar_asv[i,"Total"] <- (Polar_asv[i,"TotalAnnee"] * 100) / sum(Polar_asv$TotalAnnee)}
      ## Label
      Polar_asv <- Polar_asv %>%
        arrange(desc(Division)) %>%
        mutate(lab.ypos = cumsum(Total) - 0.5*Total)
      Polar_asv$label <- paste(round(Polar_asv$Total,1), "%", sep = "")
      for (i in rownames(Polar_asv)) {
        if (Polar_asv[i,"label"] == "0%") { Polar_asv[i,"label"] <- NA}}
      for (i in rownames(Polar_asv)) {
        if (is.na(Polar_asv[i,"label"]) == FALSE) { Polar_asv[i,"label"] <- paste(Polar_asv[i,"Division"]," : ",Polar_asv[i,"label"], sep = "")}}
      Polar_asv$Group<-l
      if (k==1) {Polar_asv_GRP<-Polar_asv}
      if (k>1) {Polar_asv_GRP <- rbind(Polar_asv_GRP,Polar_asv)}
    }
    ## Figure
    svglite("Functional-Analyse/Composition_Function/Polar-Total-asv-grp.svg",width = 10,height = 10)
    bx <- ggplot(Polar_asv_GRP, mapping = aes(y= Total, x = 2, fill = Division), Rowv = NA, col = colMain, scale = "column") +
      geom_bar(stat="identity", color = "white", width = 1) + coord_polar("y") + 
      geom_label_repel(aes(y = lab.ypos,label = label),color = "white",size = 3,segment.color = "black",show.legend = FALSE, nudge_x = 0.5,nudge_y = 0.5,max.overlaps = 20) +
      scale_y_continuous(limits=c(0,sum(Polar_asv_GRP %>% select(Total)+0.01)/length(unique(data$nbGroup))))+
      xlim(1,2.5) +
      theme_unique_darkbis() + 
      facet_wrap( ~ Group , nrow = 3) +
      labs(y = "ASVs",x="")  + scale_fill_manual(values = paletxi) #+ scale_fill_paletteer_d("ggthemes::manyeys")
    print(bx)
    dev.off()
    ## Table
    write.table(Polar_asv_GRP,file="Functional-Analyse/Composition_Function/Table-Total-asv-grp.csv",row.names = F, quote = F, sep = "\t")
    ## Table Class
    asv_table_class_grp <- data_asv_tax %>% select(Total04,Total06,Total09,Total11,TotalMonimolimnion,TotalMixolimnion,Domain,Supergroup,Division,Class,Order,RangInterest)
    asv_table_class_grp <- asv_table_class_grp %>% group_by(Domain,Supergroup,Division,Class,Order,RangInterest) %>% summarise_all(sum) %>% arrange(RangInterest)
    write.table(asv_table_class_grp, file="Functional-Analyse/Composition_Function/Table-asv-grp.csv",row.names = F, quote = F,sep="\t")
    #
  # Coplot ------------------------------------------------------------------
    svglite("Functional-Analyse/Composition_Function/Totalcompo.svg",width = 8,height = 5)
    Total_fig <- plot_grid(bxim,axim, ncol = 2, nrow = 1, labels =c("A","B"), rel_widths = c(3,3),rel_heights = c(3))
    print(Total_fig)
    dev.off()
    #
# Small vs Large ----------------------------------------------------------
    data_seq_tax_Small <- data_seq_tax %>% filter(Small==1)
    data_seq_tax_Large <- data_seq_tax %>% filter(Large==2)
    data_asv_tax_Small <- data_asv_tax %>% filter(Small==1)
    data_asv_tax_Large <- data_asv_tax %>% filter(Large==2)
    #
    palet_tree_grp <- c("#ad0635FF","#70a830FF","#98d048FF","#d84040FF","#60d0a0FF","#4080c0FF","#d7c11cFF","#b4d2d2")
    
  # Hist -------------------------------------------------------------------
    for (type in c("ASV","Sequence")) {
      if (type == "ASV") { input_data_table <- asv_mat_rare }
      if (type == "Sequence") { input_data_table <- seq_mat_rare }
      MixolimnionXSmall <- samples_df %>% filter(Zone == "Mixolimnion") %>% filter(Fraction == "Small") %>% filter(Replicat == 1)
      MixolimnionXSmall<-MixolimnionXSmall$Amplicon
      MixolimnionXSmall_seq <- as.data.frame(rowSums(input_data_table %>% select(all_of(MixolimnionXSmall))))
      colnames(MixolimnionXSmall_seq) <- "MixolimnionXSmall"
      MixolimnionXSmall_seq$ASV_Id <- row.names(MixolimnionXSmall_seq)
      MixolimnionXLarge <- samples_df %>% filter(Zone == "Mixolimnion") %>% filter(Fraction == "Large") %>% filter(Replicat == 1)
      MixolimnionXLarge<-MixolimnionXLarge$Amplicon
      MixolimnionXLarge_seq <- as.data.frame(rowSums(input_data_table %>% select(all_of(MixolimnionXLarge))))
      colnames(MixolimnionXLarge_seq) <- "MixolimnionXLarge"
      MixolimnionXLarge_seq$ASV_Id <- row.names(MixolimnionXLarge_seq)
      MonimolimnionXSmall <- samples_df %>% filter(Zone == "Monimolimnion") %>% filter(Fraction == "Small") %>% filter(Replicat == 1)
      MonimolimnionXSmall<-MonimolimnionXSmall$Amplicon
      MonimolimnionXSmall_seq <- as.data.frame(rowSums(input_data_table %>% select(all_of(MonimolimnionXSmall))))
      colnames(MonimolimnionXSmall_seq) <- "MonimolimnionXSmall"
      MonimolimnionXSmall_seq$ASV_Id <- row.names(MonimolimnionXSmall_seq)
      MonimolimnionXLarge <- samples_df %>% filter(Zone == "Monimolimnion") %>% filter(Fraction == "Large") %>% filter(Replicat == 1)
      MonimolimnionXLarge<-MonimolimnionXLarge$Amplicon
      MonimolimnionXLarge_seq <- as.data.frame(rowSums(input_data_table %>% select(all_of(MonimolimnionXLarge))))
      colnames(MonimolimnionXLarge_seq) <- "MonimolimnionXLarge"
      MonimolimnionXLarge_seq$ASV_Id <- row.names(MonimolimnionXLarge_seq)
      MixolimnionXFraction_seq <- merge(x = MixolimnionXSmall_seq,y = MixolimnionXLarge_seq, by = "ASV_Id")
      AnxiqueXFraction_seq <- merge(x = MonimolimnionXSmall_seq,y = MonimolimnionXLarge_seq, by = "ASV_Id")
      ZoneXFraction_seq <- merge(x = MixolimnionXFraction_seq,y = AnxiqueXFraction_seq, by = "ASV_Id")
      if (type == "ASV") { ZoneXFraction_seq[,2:5][ZoneXFraction_seq[,2:5]>1]<-1}
      ZoneXFraction_seq_tax <- merge(data_seq_tax %>% select(ASV_Id,RangInterest), ZoneXFraction_seq, by = "ASV_Id")
      ZoneXFraction_seq_tax <- ZoneXFraction_seq_tax %>% filter(!RangInterest %in% c("Lichen","Other multicellular Fungi", "Metazoa","Embryophyceae")) # Modif AM21102022
      ZoneXFraction_seq_tax <- ZoneXFraction_seq_tax %>% select(-ASV_Id)
      ZoneXFraction_seq_tax <- ZoneXFraction_seq_tax %>% group_by(RangInterest) %>% summarise_all(sum)
      ZoneXFraction_seq_tax_melt <-melt(ZoneXFraction_seq_tax,id.vars = c("RangInterest"))
      ZoneXFraction_seq_tax_melt$Sum <- 0
      ZoneXFraction_seq_tax_melt$Prop <- 0
      for (i in row.names(ZoneXFraction_seq_tax_melt)) { var <- ZoneXFraction_seq_tax_melt[i,"variable"]
      ZoneXFraction_seq_tax_melt[i,"Sum"] <- sum(ZoneXFraction_seq_tax_melt %>% filter(variable == var) %>% select(value))
      ZoneXFraction_seq_tax_melt[i,"Proportion"] <- ZoneXFraction_seq_tax_melt[i,"value"]*100/ZoneXFraction_seq_tax_melt[i,"Sum"]}
      ZoneXFraction_seq_tax_melt <- ZoneXFraction_seq_tax_melt %>% mutate()
      ZoneXFraction_seq_tax_melt <- separate(ZoneXFraction_seq_tax_melt,"variable",c("Zone","Fraction"),sep="X")
      assign(paste("ZoneXFraction_seq_tax_melt",type,sep = "_"),ZoneXFraction_seq_tax_melt)
      if (type == "ASV") { ZoneXFraction_asv <- ZoneXFraction_seq}
    }
    ## Plot Sequence
    svglite("Functional-Analyse/Hist_Function/ZoneXFraction-Sequence-Total.svg",width = 8.00,height = 6.00)
    Total_ZoneXFraction_Sequence_fig <- ggplot(ZoneXFraction_seq_tax_melt_Sequence, mapping = aes(y= Proportion, x = Zone, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
      scale_fill_paletteer_d("ggthemes::manyeys") + facet_grid(.~Fraction,scales="free") +
      labs(x="Zones",y="Sequences %") + geom_label(aes(y = 106,label = paste(Sum,"sequences",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
    print(Total_ZoneXFraction_Sequence_fig)
    dev.off()
    ## Plot ASV
    svglite("Functional-Analyse/Hist_Function/ZoneXFraction-ASV-Total.svg",width = 8.00,height = 6.00)
    Total_ZoneXFraction_ASV_fig <- ggplot(ZoneXFraction_seq_tax_melt_ASV, mapping = aes(y= Proportion, x = Zone, fill = RangInterest, group = RangInterest), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") +
      scale_fill_paletteer_d("ggthemes::manyeys") + facet_grid(.~Fraction,scales="free") +
      labs(x="Zones",y="ASV %") + geom_label(aes(y = 106,label = paste(Sum,"ASVs",sep ="\n")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
    print(Total_ZoneXFraction_ASV_fig)
    dev.off()
    #
  # Tree --------------------------------------------------------------------
    for (tv in c("ASV","Sequence")) {
      table_tree <- get(paste("ZoneXFraction_seq_tax_melt",tv,sep="_"))
      table_tree$labelx<-NA
      for (i in row.names(table_tree)) {
        if (is.na(table_tree[i,"Proportion"])==TRUE) {table_tree[i,"labelx"]<-NA}
        else { table_tree[i,"labelx"] <- paste0(round(table_tree[i,"Proportion"],2),"%")}}
      q <-ggplot(table_tree, aes(area = Proportion, fill = RangInterest, label = paste0(round(Proportion,1),"%"))) +
        geom_treemap(alpha=0.80,color="white",size=2) +
        scale_fill_manual(values = palet_tree_grp) + 
        facet_grid(Zone~Fraction) +
        theme_unique_darktris() +
        theme(legend.position="bottom") +
        geom_treemap_text(fontface = "bold", colour = "black", place = "centre",grow = FALSE,size=12)
      svglite(paste(paste0("Functional-Analyse/Tree_Function/","ZoneXFraction"),tv,"Tree_Lab.svg",sep="_"),width = 12.00,height = 6.00)
      print(q)
      dev.off()
      r <-ggplot(table_tree, aes(area = Proportion, fill = RangInterest, label = "")) +
        geom_treemap(alpha=0.80,color="white",size=2) +
        scale_fill_manual(values = palet_tree_grp) + 
        facet_grid(Zone~Fraction) +
        theme_unique_darktris() +
        theme(legend.position="bottom") +
        geom_treemap_text(fontface = "bold", colour = "black", place = "centre",grow = FALSE,size=12)
      svglite(paste(paste0("Functional-Analyse/Tree_Function/","ZoneXFraction"),tv,"Tree.svg",sep="_"),width = 12.00,height = 6.00)
      print(r)
      dev.off()
    }
    #
# Tree-2 ------------------------------------------------------------------
  # Sequence ----------------------------------------------------------------
    Group_FraXZo_seq_tax <- merge(data_seq_tax %>% select(ASV_Id,RangInterest,Supergroup,Division,Class,Order,Family,Genus,Last), ZoneXFraction_seq, by = "ASV_Id")
    Group_FraXZo_seq_tax <- Group_FraXZo_seq_tax %>% filter(!RangInterest %in% c("Lichen","Other multicellular Fungi", "Metazoa","Embryophyceae")) # Modif AM21102022
    Group_FraXZo_seq_tax <- Group_FraXZo_seq_tax %>% select(-ASV_Id) %>% 
      group_by(RangInterest,Supergroup,Division,Class,Order,Family,Genus,Last) %>%
      summarise_all(sum)
    Group_FraXZo_seq_tax_melt <- melt(Group_FraXZo_seq_tax, id.vars = c("RangInterest","Supergroup","Division","Class","Order","Family","Genus","Last"))
    Group_FraXZo_seq_tax_melt <- separate(Group_FraXZo_seq_tax_melt,variable,c("Zone","Fraction"),"X")
    color_data_X <- as.data.frame(palet_tree_grp)
    color_data_X$RangInterest<-(unique(Group_FraXZo_seq_tax_melt$RangInterest))
    Group_FraXZo_seq_tax_melt <- merge(Group_FraXZo_seq_tax_melt,color_data_X,by="RangInterest")
    for (i in c("Supergroup","Division","Class","Order","Family","Genus","Last")) {
      Group_FraXZo_seq_tax_melt[,i][Group_FraXZo_seq_tax_melt[,"RangInterest"]=="Unassigned"]<-""}
    ## Table
    write.table(Group_FraXZo_seq_tax_melt, file = "Functional-Analyse/Table/ZoneXFractionXGenus_Sequence_avgroup_Tree.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    ## Family
    #Group_FraXZo_seq_tax_melt <- Group_FraXZo_seq_tax_melt %>% select(-Last, -Genus) %>% group_by(RangInterest,Supergroup,Division,Class,Order,Family,Fraction,Zone,palet_tree_grp) %>% summarise_all(sum)
    ## Order
    #Group_FraXZo_seq_tax_melt <- Group_FraXZo_seq_tax_melt %>% select(-Last, -Genus, -Family) %>% group_by(RangInterest,Supergroup,Division,Class,Order,Fraction,Zone,palet_tree_grp) %>% summarise_all(sum)
    ## Class
    Group_FraXZo_seq_tax_melt <- Group_FraXZo_seq_tax_melt %>% select(-Last, -Genus, -Family, -Order) %>% group_by(RangInterest,Supergroup,Division,Class,Fraction,Zone,palet_tree_grp) %>% summarise_all(sum)
    ### Plot
    legendXi<-get_legend(r)
    s <-ggplot(Group_FraXZo_seq_tax_melt, aes(area = value, fill = palet_tree_grp, subgroup = RangInterest, label = Class)) +
      geom_treemap(alpha=0.80,color="white",size=2,fill=Group_FraXZo_seq_tax_melt$palet_tree_grp) +
      #scale_fill_manual(values = palet_tree_grp) + 
      facet_grid(Zone~Fraction) +
      geom_treemap_subgroup_border(size=2,color="black") +
      #geom_treemap_subgroup_text(place = "centre", grow = T, alpha = 0.5, colour = "white",fontface = "italic", min.size = 0) +
      geom_treemap_text(fontface = "bold", colour = "black", place = "centre",grow = FALSE,size=12,min.size=4)+
      theme_unique_darktris() + theme(panel.border = element_rect(fill=NA,size=2),legend.position="none")
    t<-plot_grid(s,legendXi, ncol = 1, nrow = 2, rel_heights = c(5,1))
    svglite(paste(paste0("Functional-Analyse/Tree2_Function/","ZoneXFractionXGenus"),"Sequence","Tree.svg",sep="_"),width = 14.00,height = 8.00)
    print(t)
    dev.off()
    write.table(Group_FraXZo_seq_tax_melt, file = "Functional-Analyse/Table/ZoneXFractionXGenus_Sequence_apgroup_Tree.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    ## Poster _ Monimolimnion vs Mixo
    Group_FraXZo_seq_tax_melt <- as.data.frame(Group_FraXZo_seq_tax_melt)
    Group_FraXZo_seq_tax_melt_Poster <- Group_FraXZo_seq_tax_melt %>% select(-Fraction) %>% group_by(RangInterest,Supergroup,Division,Class,Zone,palet_tree_grp) %>% summarise_all(sum)
    ### Plot
    s <-ggplot(Group_FraXZo_seq_tax_melt_Poster, aes(area = value, fill = palet_tree_grp, subgroup = RangInterest, label = Class)) +
      geom_treemap(alpha=0.80,color="white",size=2,fill=Group_FraXZo_seq_tax_melt_Poster$palet_tree_grp) +
      #scale_fill_manual(values = palet_tree_grp) + 
      facet_grid(.~Zone) +
      geom_treemap_subgroup_border(size=2,color="black") +
      #geom_treemap_subgroup_text(place = "centre", grow = T, alpha = 0.5, colour = "white",fontface = "italic", min.size = 0) +
      geom_treemap_text(fontface = "bold", colour = "black", place = "centre",grow = FALSE,size=12,min.size=4)+
      theme_unique_darktris() + theme(panel.border = element_rect(fill=NA,size=2),legend.position="none")
    t<-plot_grid(s,legendXi, ncol = 1, nrow = 2, rel_heights = c(5,1))
    svglite(paste(paste0("Functional-Analyse/Tree2_Function/","ZoneXFractionXGenus"),"SEQ","Tree_POSTER.svg",sep="_"),width = 14.00,height = 6.00)
    print(t)
    dev.off()
    #
  # ASV ---------------------------------------------------------------------
    Group_FraXZo_asv_tax <- merge(data_seq_tax %>% select(ASV_Id,RangInterest,Supergroup,Division,Class,Order,Family,Genus,Last), ZoneXFraction_asv, by = "ASV_Id")
    Group_FraXZo_asv_tax <- Group_FraXZo_asv_tax %>% filter(!RangInterest %in% c("Lichen","Other multicellular Fungi", "Metazoa","Embryophyceae")) # Modif AM21102022
    
    Group_FraXZo_asv_tax <- Group_FraXZo_asv_tax %>% select(-ASV_Id) %>% 
      group_by(RangInterest,Supergroup,Division,Class,Order,Family,Genus,Last) %>%
      summarise_all(sum)
    Group_FraXZo_asv_tax_melt <- melt(Group_FraXZo_asv_tax, id.vars = c("RangInterest","Supergroup","Division","Class","Order","Family","Genus","Last"))
    Group_FraXZo_asv_tax_melt <- separate(Group_FraXZo_asv_tax_melt,variable,c("Zone","Fraction"),"X")
    color_data_X <- as.data.frame(palet_tree_grp)
    color_data_X$RangInterest<-(unique(Group_FraXZo_asv_tax_melt$RangInterest))
    Group_FraXZo_asv_tax_melt <- merge(Group_FraXZo_asv_tax_melt,color_data_X,by="RangInterest")
    
    for (i in c("Supergroup","Division","Class","Order","Family","Genus","Last")) {
      Group_FraXZo_asv_tax_melt[,i][Group_FraXZo_asv_tax_melt[,"RangInterest"]=="Unassigned"]<-""}
    write.table(Group_FraXZo_asv_tax_melt, file = "Functional-Analyse/Table/ZoneXFractionXGenus_ASV_avgroup_Tree.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    ## Family
    #Group_FraXZo_asv_tax_melt <- Group_FraXZo_asv_tax_melt %>% select(-Last, -Genus) %>% group_by(RangInterest,Supergroup,Division,Class,Order,Family,Fraction,Zone,palet_tree_grp) %>% summarise_all(sum)
    ## Order
    #Group_FraXZo_asv_tax_melt <- Group_FraXZo_asv_tax_melt %>% select(-Last, -Genus, -Family) %>% group_by(RangInterest,Supergroup,Division,Class,Order,Fraction,Zone,palet_tree_grp) %>% summarise_all(sum)
    ## Class
    Group_FraXZo_asv_tax_melt <- Group_FraXZo_asv_tax_melt %>% select(-Last, -Genus, -Family, -Order) %>% group_by(RangInterest,Supergroup,Division,Class,Fraction,Zone,palet_tree_grp) %>% summarise_all(sum)
    ### Plot
    legendXi<-get_legend(r)
    s <-ggplot(Group_FraXZo_asv_tax_melt, aes(area = value, fill = palet_tree_grp, subgroup = RangInterest, label = Class)) +
      geom_treemap(alpha=0.80,color="white",size=2,fill=Group_FraXZo_asv_tax_melt$palet_tree_grp) +
      facet_grid(Zone~Fraction) +
      geom_treemap_subgroup_border(size=2,color="black") +
      geom_treemap_text(fontface = "bold", colour = "black", place = "centre",grow = FALSE,size=12,min.size=4)+
      theme_unique_darktris() + theme(panel.border = element_rect(fill=NA,size=2),legend.position="none")
    t<-plot_grid(s,legendXi, ncol = 1, nrow = 2, rel_heights = c(5,1))
    svglite(paste(paste0("Functional-Analyse/Tree2_Function/","ZoneXFractionXGenus"),"ASV","Tree.svg",sep="_"),width = 14.00,height = 8.00)
    print(t)
    dev.off()
    write.table(Group_FraXZo_asv_tax_melt, file = "Functional-Analyse/Table/ZoneXFractionXGenus_ASV_apgroup_Tree.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    #
    ## Poster ASV _ Monimolimnion vs Mixo
    Group_FraXZo_asv_tax_melt <- as.data.frame(Group_FraXZo_asv_tax_melt)
    Group_FraXZo_asv_tax_melt_Poster <- Group_FraXZo_asv_tax_melt %>% select(-Fraction) %>% group_by(RangInterest,Supergroup,Division,Class,Zone,palet_tree_grp) %>% summarise_all(sum)
    ### Plot
    s <-ggplot(Group_FraXZo_asv_tax_melt_Poster, aes(area = value, fill = palet_tree_grp, subgroup = RangInterest, label = Class)) +
      geom_treemap(alpha=0.80,color="white",size=2,fill=Group_FraXZo_asv_tax_melt_Poster$palet_tree_grp) +
      #scale_fill_manual(values = palet_tree_grp) + 
      facet_grid(.~Zone) +
      geom_treemap_subgroup_border(size=2,color="black") +
      #geom_treemap_subgroup_text(place = "centre", grow = T, alpha = 0.5, colour = "white",fontface = "italic", min.size = 0) +
      geom_treemap_text(fontface = "bold", colour = "black", place = "centre",grow = FALSE,size=12,min.size=4)+
      theme_unique_darktris() + theme(panel.border = element_rect(fill=NA,size=2),legend.position="none")
    t<-plot_grid(s,legendXi, ncol = 1, nrow = 2, rel_heights = c(5,1))
    svglite(paste(paste0("Functional-Analyse/Tree2_Function/","ZoneXFractionXGenus"),"ASV","Tree_POSTER.svg",sep="_"),width = 14.00,height = 6.00)
    print(t)
    dev.off()
# Periods Tree_4 ----------------------------------------------------------
  # Sequence ----------------------------------------------------------------
    Group_Periods_seq_tax <- data_seq_tax %>% select(ASV_Id,RangInterest,Supergroup,Division,Class,Order,Family,Genus,Last, Total04,Total06,Total09,Total11)
    Group_Periods_seq_tax <- Group_Periods_seq_tax %>% filter(!RangInterest %in% c("Lichen","Other multicellular Fungi", "Metazoa","Embryophyceae")) # Modif AM21102022
    Group_Periods_seq_tax <- Group_Periods_seq_tax %>% select(-ASV_Id) %>% 
      group_by(RangInterest,Supergroup,Division,Class,Order,Family,Genus,Last) %>%
      summarise_all(sum)
    Group_Periods_seq_tax_melt <- melt(Group_Periods_seq_tax, id.vars = c("RangInterest","Supergroup","Division","Class","Order","Family","Genus","Last"))
    color_data_X <- as.data.frame(palet_tree_grp)
    color_data_X$RangInterest <- unique(Group_Periods_seq_tax_melt$RangInterest)
    Group_Periods_seq_tax_melt <- merge(Group_Periods_seq_tax_melt,color_data_X,by="RangInterest")
    for (i in c("Supergroup","Division","Class","Order","Family","Genus","Last")) {
      Group_Periods_seq_tax_melt[,i][Group_Periods_seq_tax_melt[,"RangInterest"]=="Unassigned"]<-""}
    ## Table
    write.table(Group_Periods_seq_tax_melt, file = "Functional-Analyse/Table/PeriodsXGenus_Sequence_avgroup_Tree.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    ## Family
    #Group_FraXZo_seq_tax_melt <- Group_FraXZo_seq_tax_melt %>% select(-Last, -Genus) %>% group_by(RangInterest,Supergroup,Division,Class,Order,Family,Fraction,Zone,palet_tree_grp) %>% summarise_all(sum)
    ## Order
    #Group_FraXZo_seq_tax_melt <- Group_FraXZo_seq_tax_melt %>% select(-Last, -Genus, -Family) %>% group_by(RangInterest,Supergroup,Division,Class,Order,Fraction,Zone,palet_tree_grp) %>% summarise_all(sum)
    ## Class
    Group_Periods_seq_tax_melt <- Group_Periods_seq_tax_melt %>% select(-Last, -Genus, -Family, -Order) %>% group_by(RangInterest,Supergroup,Division,Class,variable,palet_tree_grp) %>% summarise_all(sum)
    ### Plot
    legendXi<-get_legend(r)
    s <-ggplot(Group_Periods_seq_tax_melt, aes(area = value, fill = palet_tree_grp, subgroup = RangInterest, label = Class)) +
      geom_treemap(alpha=0.80,color="white",size=2,fill=Group_Periods_seq_tax_melt$palet_tree_grp) +
      #scale_fill_manual(values = palet_tree_grp) + 
      facet_grid(variable~.) +
      geom_treemap_subgroup_border(size=2,color="black") +
      #geom_treemap_subgroup_text(place = "centre", grow = T, alpha = 0.5, colour = "white",fontface = "italic", min.size = 0) +
      geom_treemap_text(fontface = "bold", colour = "black", place = "centre",grow = FALSE,size=12,min.size=4)+
      theme_unique_darktris() + theme(panel.border = element_rect(fill=NA,size=2),legend.position="none") +
      labs(title="Sequences abundance")
    t<-plot_grid(s,legendXi, ncol = 1, nrow = 2, rel_heights = c(5,1))
    svglite(paste(paste0("Functional-Analyse/Tree2_Function/","PeriodsXGenus"),"Sequence","Tree.svg",sep="_"),width = 14.00,height = 8.00)
    print(t)
    dev.off()
    write.table(Group_Periods_seq_tax_melt, file = "Functional-Analyse/Table/PeriodsXGenus_Sequence_apgroup_Tree.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    #
  # ASV ---------------------------------------------------------------------
    Group_Periods_asv_tax <- data_asv_tax %>% select(ASV_Id,RangInterest,Supergroup,Division,Class,Order,Family,Genus,Last,Total04,Total06,Total09,Total11)
    Group_Periods_asv_tax <- Group_Periods_asv_tax %>% filter(!RangInterest %in% c("Lichen","Other multicellular Fungi", "Metazoa","Embryophyceae")) # Modif AM21102022
    colnames(Group_Periods_asv_tax) <- c("ASV_Id","RangInterest","Supergroup","Division","Class","Order","Family","Genus","Last","April","June","September","November")
    Group_Periods_asv_tax <- Group_Periods_asv_tax %>% select(-ASV_Id) %>% 
      group_by(RangInterest,Supergroup,Division,Class,Order,Family,Genus,Last) %>%
      summarise_all(sum)
    Group_Periods_asv_tax_melt <- melt(Group_Periods_asv_tax, id.vars = c("RangInterest","Supergroup","Division","Class","Order","Family","Genus","Last"))
    color_data_X <- as.data.frame(palet_tree_grp)
    color_data_X$RangInterest<-(unique(Group_Periods_asv_tax_melt$RangInterest))
    Group_Periods_asv_tax_melt <- merge(Group_Periods_asv_tax_melt,color_data_X,by="RangInterest")
    
    for (i in c("Supergroup","Division","Class","Order","Family","Genus","Last")) {
      Group_Periods_asv_tax_melt[,i][Group_Periods_asv_tax_melt[,"RangInterest"]=="Unassigned"]<-""}
    write.table(Group_Periods_asv_tax_melt, file = "Functional-Analyse/Table/PeriodsXGenus_ASV_avgroup_Tree.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    # Class
    Group_Periods_asv_tax_melt <- Group_Periods_asv_tax_melt %>% select(-Last, -Genus, -Family, -Order) %>% group_by(RangInterest,Supergroup,Division,Class,variable,palet_tree_grp) %>% summarise_all(sum)
    Group_Periods_asv_tax_melt <- as.data.frame(Group_Periods_asv_tax_melt)
    ### Plot
    legendXi<-get_legend(r)
    u <-ggplot(Group_Periods_asv_tax_melt, aes(area = value, fill = palet_tree_grp, subgroup = RangInterest, label = Class)) +
      geom_treemap(alpha=0.80,color="white",size=2,fill=Group_Periods_asv_tax_melt$palet_tree_grp) +
      facet_grid(variable~.) +
      geom_treemap_subgroup_border(size=2,color="black") +
      geom_treemap_text(fontface = "bold", colour = "black", place = "centre",grow = FALSE,size=12,min.size=4)+
      theme_unique_darktris() + theme(panel.border = element_rect(fill=NA,size=2),legend.position="none") +
      labs(title="ASVs richness")
    v<-plot_grid(u,legendXi, ncol = 1, nrow = 2, rel_heights = c(5,1))
    svglite(paste(paste0("Functional-Analyse/Tree2_Function/","PeriodsXGenus"),"ASV","Tree.svg",sep="_"),width = 14.00,height = 8.00)
    print(v)
    dev.off()
    write.table(Group_Periods_asv_tax_melt, file = "Functional-Analyse/Table/PeriodsXGenus_ASV_apgroup_Tree.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    # Coplot
    svglite(paste(paste0("Functional-Analyse/Tree2_Function/","PeriodsXGenus"),"All","Tree.svg",sep="_"),width = 10.00,height = 8.00)
    plot_grid(s + 
                theme(strip.text = element_blank(), plot.title = element_text(color = "black", face = "bold", hjust = 0.5,size = 18)),
              NULL,u +
                theme(plot.title = element_text(color = "black", face = "bold", hjust = 0.5,size = 18)),
              NULL,legendXi,NULL, ncol = 3, nrow = 2, rel_heights = c(5,0.5), rel_widths = c(5,0.01,5))
    dev.off()
    
    
    
# Tree ON LAST ------------------------------------------------------------
  # Sequence ZoneXFraction ----------------------------------------------------------------
    Group_FraXZo_seq_tax <- merge(data_seq_tax %>% select(ASV_Id,RangInterest,Supergroup,Division,Class,Order,Family,Genus,Last), ZoneXFraction_seq, by = "ASV_Id")
    Group_FraXZo_seq_tax <- Group_FraXZo_seq_tax %>% filter(!RangInterest %in% c("Lichen","Other multicellular Fungi", "Metazoa","Embryophyceae")) # Modif AM21102022
    Group_FraXZo_seq_tax <- Group_FraXZo_seq_tax %>% select(-ASV_Id) %>% 
      group_by(RangInterest,Supergroup,Division,Class,Order,Family,Genus,Last) %>%
      summarise_all(sum)
    Group_FraXZo_seq_tax_melt <- melt(Group_FraXZo_seq_tax, id.vars = c("RangInterest","Supergroup","Division","Class","Order","Family","Genus","Last"))
    Group_FraXZo_seq_tax_melt <- separate(Group_FraXZo_seq_tax_melt,variable,c("Zone","Fraction"),"X")
    color_data_X <- as.data.frame(palet_tree_grp)
    color_data_X$RangInterest<-(unique(Group_FraXZo_seq_tax_melt$RangInterest))
    Group_FraXZo_seq_tax_melt <- merge(Group_FraXZo_seq_tax_melt,color_data_X,by="RangInterest")
    for (i in c("Supergroup","Division","Class","Order","Family","Genus","Last")) {
      Group_FraXZo_seq_tax_melt[,i][Group_FraXZo_seq_tax_melt[,"RangInterest"]=="Unassigned"]<-""}
    ## LAST
    Group_FraXZo_seq_tax_melt <- Group_FraXZo_seq_tax_melt %>% group_by(RangInterest,Supergroup,Division,Class,Order,Family,Genus,Last,Fraction,Zone,palet_tree_grp) %>% summarise_all(sum)
    ### Plot
    legendXi<-get_legend(r)
    s <-ggplot(Group_FraXZo_seq_tax_melt, aes(area = value, fill = palet_tree_grp, subgroup = RangInterest, label = Last)) +
      geom_treemap(alpha=0.80,color="white",size=2,fill=Group_FraXZo_seq_tax_melt$palet_tree_grp) +
      #scale_fill_manual(values = palet_tree_grp) + 
      facet_grid(Zone~Fraction) +
      geom_treemap_subgroup_border(size=2,color="black") +
      #geom_treemap_subgroup_text(place = "centre", grow = T, alpha = 0.5, colour = "white",fontface = "italic", min.size = 0) +
      geom_treemap_text(fontface = "bold", colour = "black", place = "centre",grow = FALSE,size=12,min.size=4)+
      theme_unique_darktris() + theme(panel.border = element_rect(fill=NA,size=2),legend.position="none")
    t<-plot_grid(s,legendXi, ncol = 1, nrow = 2, rel_heights = c(5,1))
    svglite(paste(paste0("Functional-Analyse/Tree2_Function/","ZoneXFractionXLast"),"Sequence","Tree.svg",sep="_"),width = 18.00,height = 10.00)
    print(t)
    dev.off()
    write.table(Group_FraXZo_seq_tax_melt, file = "Functional-Analyse/Table/ZoneXFractionXLast_Sequence_apgroup_Tree.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    
    
    
  # ASV ZoneXFraction ---------------------------------------------------------------------
    Group_FraXZo_asv_tax <- merge(data_seq_tax %>% select(ASV_Id,RangInterest,Supergroup,Division,Class,Order,Family,Genus,Last), ZoneXFraction_asv, by = "ASV_Id")
    Group_FraXZo_asv_tax <- Group_FraXZo_asv_tax %>% filter(!RangInterest %in% c("Lichen","Other multicellular Fungi", "Metazoa","Embryophyceae")) # Modif AM21102022
    
    Group_FraXZo_asv_tax <- Group_FraXZo_asv_tax %>% select(-ASV_Id) %>% 
      group_by(RangInterest,Supergroup,Division,Class,Order,Family,Genus,Last) %>%
      summarise_all(sum)
    Group_FraXZo_asv_tax_melt <- melt(Group_FraXZo_asv_tax, id.vars = c("RangInterest","Supergroup","Division","Class","Order","Family","Genus","Last"))
    Group_FraXZo_asv_tax_melt <- separate(Group_FraXZo_asv_tax_melt,variable,c("Zone","Fraction"),"X")
    color_data_X <- as.data.frame(palet_tree_grp)
    color_data_X$RangInterest<-(unique(Group_FraXZo_asv_tax_melt$RangInterest))
    Group_FraXZo_asv_tax_melt <- merge(Group_FraXZo_asv_tax_melt,color_data_X,by="RangInterest")
    for (i in c("Supergroup","Division","Class","Order","Family","Genus","Last")) {
      Group_FraXZo_asv_tax_melt[,i][Group_FraXZo_asv_tax_melt[,"RangInterest"]=="Unassigned"]<-""}
    ## Last
    Group_FraXZo_asv_tax_melt <- Group_FraXZo_asv_tax_melt %>% group_by(RangInterest,Supergroup,Division,Class,Last,Genus,Family,Order,Fraction,Zone,palet_tree_grp) %>% summarise_all(sum)
    ### Plot
    legendXi<-get_legend(r)
    s <-ggplot(Group_FraXZo_asv_tax_melt, aes(area = value, fill = palet_tree_grp, subgroup = RangInterest, label = Last)) +
      geom_treemap(alpha=0.80,color="white",size=2,fill=Group_FraXZo_asv_tax_melt$palet_tree_grp) +
      facet_grid(Zone~Fraction) +
      geom_treemap_subgroup_border(size=2,color="black") +
      geom_treemap_text(fontface = "bold", colour = "black", place = "centre",grow = FALSE,size=12,min.size=4)+
      theme_unique_darktris() + theme(panel.border = element_rect(fill=NA,size=2),legend.position="none")
    t<-plot_grid(s,legendXi, ncol = 1, nrow = 2, rel_heights = c(5,1))
    svglite(paste(paste0("Functional-Analyse/Tree2_Function/","ZoneXFractionXLast"),"ASV","Tree.svg",sep="_"),width = 18.00,height = 10.00)
    print(t)
    dev.off()
    write.table(Group_FraXZo_asv_tax_melt, file = "Functional-Analyse/Table/ZoneXFractionXLast_ASV_apgroup_Tree.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    #
    
  # Sequence Periods ----------------------------------------------------------------
    Group_Periods_seq_tax <- data_seq_tax %>% select(ASV_Id,RangInterest,Supergroup,Division,Class,Order,Family,Genus,Last, Total04,Total06,Total09,Total11)
    Group_Periods_seq_tax <- Group_Periods_seq_tax %>% filter(!RangInterest %in% c("Lichen","Other multicellular Fungi", "Metazoa","Embryophyceae")) # Modif AM21102022
    Group_Periods_seq_tax <- Group_Periods_seq_tax %>% select(-ASV_Id) %>% 
      group_by(RangInterest,Supergroup,Division,Class,Order,Family,Genus,Last) %>%
      summarise_all(sum)
    Group_Periods_seq_tax_melt <- melt(Group_Periods_seq_tax, id.vars = c("RangInterest","Supergroup","Division","Class","Order","Family","Genus","Last"))
    color_data_X <- as.data.frame(palet_tree_grp)
    color_data_X$RangInterest <- unique(Group_Periods_seq_tax_melt$RangInterest)
    Group_Periods_seq_tax_melt <- merge(Group_Periods_seq_tax_melt,color_data_X,by="RangInterest")
    for (i in c("Supergroup","Division","Class","Order","Family","Genus","Last")) {
      Group_Periods_seq_tax_melt[,i][Group_Periods_seq_tax_melt[,"RangInterest"]=="Unassigned"]<-""}
    ## Last
    Group_Periods_seq_tax_melt <- Group_Periods_seq_tax_melt  %>% group_by(RangInterest,Supergroup,Division,Class,Last,Genus,Family,Order,variable,palet_tree_grp) %>% summarise_all(sum)
    ### Plot
    legendXi<-get_legend(r)
    s <-ggplot(Group_Periods_seq_tax_melt, aes(area = value, fill = palet_tree_grp, subgroup = RangInterest, label = Last)) +
      geom_treemap(alpha=0.80,color="white",size=2,fill=Group_Periods_seq_tax_melt$palet_tree_grp) +
      #scale_fill_manual(values = palet_tree_grp) + 
      facet_grid(variable~.) +
      geom_treemap_subgroup_border(size=2,color="black") +
      #geom_treemap_subgroup_text(place = "centre", grow = T, alpha = 0.5, colour = "white",fontface = "italic", min.size = 0) +
      geom_treemap_text(fontface = "bold", colour = "black", place = "centre",grow = FALSE,size=12,min.size=4)+
      theme_unique_darktris() + theme(panel.border = element_rect(fill=NA,size=2),legend.position="none") +
      labs(title="Sequences abundance")
    t<-plot_grid(s,legendXi, ncol = 1, nrow = 2, rel_heights = c(5,1))
    svglite(paste(paste0("Functional-Analyse/Tree2_Function/","PeriodsXLast"),"Sequence","Tree.svg",sep="_"),width = 14.00,height = 8.00)
    print(t)
    dev.off()
    write.table(Group_Periods_seq_tax_melt, file = "Functional-Analyse/Table/PeriodsXLast_Sequence_apgroup_Tree.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    #
  # ASV Periods ---------------------------------------------------------------------
    Group_Periods_asv_tax <- data_asv_tax %>% select(ASV_Id,RangInterest,Supergroup,Division,Class,Order,Family,Genus,Last,Total04,Total06,Total09,Total11)
    Group_Periods_asv_tax <- Group_Periods_asv_tax %>% filter(!RangInterest %in% c("Lichen","Other multicellular Fungi", "Metazoa","Embryophyceae")) # Modif AM21102022
    colnames(Group_Periods_asv_tax) <- c("ASV_Id","RangInterest","Supergroup","Division","Class","Order","Family","Genus","Last","April","June","September","November")
    Group_Periods_asv_tax <- Group_Periods_asv_tax %>% select(-ASV_Id) %>% 
      group_by(RangInterest,Supergroup,Division,Class,Order,Family,Genus,Last) %>%
      summarise_all(sum)
    Group_Periods_asv_tax_melt <- melt(Group_Periods_asv_tax, id.vars = c("RangInterest","Supergroup","Division","Class","Order","Family","Genus","Last"))
    color_data_X <- as.data.frame(palet_tree_grp)
    color_data_X$RangInterest<-(unique(Group_Periods_asv_tax_melt$RangInterest))
    Group_Periods_asv_tax_melt <- merge(Group_Periods_asv_tax_melt,color_data_X,by="RangInterest")
    for (i in c("Supergroup","Division","Class","Order","Family","Genus","Last")) {
      Group_Periods_asv_tax_melt[,i][Group_Periods_asv_tax_melt[,"RangInterest"]=="Unassigned"]<-""}
    # Class
    Group_Periods_asv_tax_melt <- Group_Periods_asv_tax_melt %>% group_by(RangInterest,Supergroup,Division,Class,Last,Genus,Family,Order,variable,palet_tree_grp) %>% summarise_all(sum)
    #Group_Periods_asv_tax_melt <- as.data.frame(Group_Periods_asv_tax_melt)
    ### Plot
    legendXi<-get_legend(r)
    u <-ggplot(Group_Periods_asv_tax_melt, aes(area = value, fill = palet_tree_grp, subgroup = RangInterest, label = Last)) +
      geom_treemap(alpha=0.80,color="white",size=2,fill=Group_Periods_asv_tax_melt$palet_tree_grp) +
      facet_grid(variable~.) +
      geom_treemap_subgroup_border(size=2,color="black") +
      geom_treemap_text(fontface = "bold", colour = "black", place = "centre",grow = FALSE,size=12,min.size=4)+
      theme_unique_darktris() + theme(panel.border = element_rect(fill=NA,size=2),legend.position="none") +
      labs(title="ASVs richness")
    v<-plot_grid(u,legendXi, ncol = 1, nrow = 2, rel_heights = c(5,1))
    svglite(paste(paste0("Functional-Analyse/Tree2_Function/","PeriodsXLast"),"ASV","Tree.svg",sep="_"),width = 14.00,height = 8.00)
    print(v)
    dev.off()
    write.table(Group_Periods_asv_tax_melt, file = "Functional-Analyse/Table/PeriodsXLast_ASV_apgroup_Tree.csv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    # Coplot
    svglite(paste(paste0("Functional-Analyse/Tree2_Function/","PeriodsXLast"),"All","Tree.svg",sep="_"),width = 18.00,height = 10.00)
    plot_grid(s + 
                theme(strip.text = element_blank(), plot.title = element_text(color = "black", face = "bold", hjust = 0.5,size = 18)),
              NULL,u +
                theme(plot.title = element_text(color = "black", face = "bold", hjust = 0.5,size = 18)),
              NULL,legendXi,NULL, ncol = 3, nrow = 2, rel_heights = c(5,0.5), rel_widths = c(5,0.01,5))
    dev.off()
    
    
    
# Table ASVs majoritaires -------------------------------------------------
  # Table ONLY---------------------------------------------------------------------
    # Cycle-100 -------------------------------------------------------------------
    ## Day
    onlyDayTable <- data_seq_tax %>% filter(Cycle == "Day") %>% select(ASV_Id,TotalDay,RangInterest,Cycle)
    onlyDayTable <- onlyDayTable %>% filter(TotalDay > 0.001*sum(onlyDayTable$TotalDay))
    colnames(onlyDayTable)[2]  <- "value"
    ## Night
    onlyNightTable <- data_seq_tax %>% filter(Cycle == "Night") %>% select(ASV_Id,TotalNight,RangInterest,Cycle)
    onlyNightTable <- onlyNightTable %>% filter(TotalNight > 0.001*sum(onlyNightTable$TotalNight))
    colnames(onlyNightTable)[2]  <- "value"
    ## Bind
    onlyCycleTable <- rbind(onlyDayTable,onlyNightTable)
    onlyCycleTable <- merge(onlyCycleTable,profile_tax_table %>% select(-Last), by = "ASV_Id")
    ## Table
    write.table(onlyCycleTable, file = "Functional-Analyse/Table/Only_Cycles.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(Only_Cycle_ASV, file = "Functional-Analyse/Table/Only_Cycles_ASV.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(Only_Cycle_Sequence, file = "Functional-Analyse/Table/Only_Cycles_Sequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    #
    # Zone-100 -------------------------------------------------------------------
    ## Mixolimnion
    onlyMixolimnionTable <- data_seq_tax %>% filter(Zone == "Mixolimnion") %>% select(ASV_Id,TotalMixolimnion,RangInterest,Zone)
    onlyMixolimnionTable <- onlyMixolimnionTable %>% filter(TotalMixolimnion > 0.001*sum(onlyMixolimnionTable$TotalMixolimnion))
    colnames(onlyMixolimnionTable)[2]  <- "value"
    ## Monimolimnion
    onlyMonimolimnionTable <- data_seq_tax %>% filter(Zone == "Monimolimnion") %>% select(ASV_Id,TotalMonimolimnion,RangInterest,Zone)
    onlyMonimolimnionTable <- onlyMonimolimnionTable %>% filter(TotalMonimolimnion > 0.001*sum(onlyMonimolimnionTable$TotalMonimolimnion))
    colnames(onlyMonimolimnionTable)[2]  <- "value"
    ## Bind
    onlyZoneTable <- rbind(onlyMixolimnionTable,onlyMonimolimnionTable)
    onlyZoneTable <- merge(onlyZoneTable,profile_tax_table %>% select(-Last), by = "ASV_Id")
    ## Table
    write.table(onlyZoneTable, file = "Functional-Analyse/Table/Only_Zone.txt", sep = "\t", col.names = TRUE, row.names = FALSE,quote = FALSE)
    write.table(Only_Zone_ASV, file = "Functional-Analyse/Table/Only_Zone_ASV.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(Only_Zone_Sequence, file = "Functional-Analyse/Table/Only_Zone_Sequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    #
    # Fraction-100 -------------------------------------------------------------------
    ## Small
    onlySmallTable <- data_seq_tax %>% filter(Fraction == "Small") %>% select(ASV_Id,TotalSmall,RangInterest,Fraction)
    onlySmallTable <- onlySmallTable %>% filter(TotalSmall > 0.001*sum(onlySmallTable$TotalSmall))
    colnames(onlySmallTable)[2]  <- "value"
    ## Large
    onlyLargeTable <- data_seq_tax %>% filter(Fraction == "Large") %>% select(ASV_Id,TotalLarge,RangInterest,Fraction)
    onlyLargeTable <- onlyLargeTable %>% filter(TotalLarge > 0.001*sum(onlyLargeTable$TotalLarge))
    colnames(onlyLargeTable)[2]  <- "value"
    ## Bind
    onlyFractionTable <- rbind(onlySmallTable,onlyLargeTable)
    onlyFractionTable <- merge(onlyFractionTable,profile_tax_table %>% select(-Last), by = "ASV_Id")
    ## Table
    write.table(onlyFractionTable, file = "Functional-Analyse/Table/Only_Fraction.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(Only_Fraction_ASV, file = "Functional-Analyse/Table/Only_Fraction_ASV.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(Only_Fraction_Sequence, file = "Functional-Analyse/Table/Only_Fraction_Sequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    #
    # Cycle-90 -------------------------------------------------------------------
    ## Day
    onlyDayTable90 <- data_seq_tax %>% filter(Cycle90 == "Day") %>% select(ASV_Id,TotalDay,RangInterest,Cycle90)
    onlyDayTable90 <- onlyDayTable90 %>% filter(TotalDay > 0.001*sum(onlyDayTable90$TotalDay))
    colnames(onlyDayTable90)[2]  <- "value"
    ## Night
    onlyNightTable90 <- data_seq_tax %>% filter(Cycle90 == "Night") %>% select(ASV_Id,TotalNight,RangInterest,Cycle90)
    onlyNightTable90 <- onlyNightTable90 %>% filter(TotalNight > 0.001*sum(onlyNightTable90$TotalNight))
    colnames(onlyNightTable90)[2]  <- "value"
    ## Bind
    onlyCycleTable90 <- rbind(onlyDayTable90,onlyNightTable90)
    onlyCycleTable90 <- merge(onlyCycleTable90,profile_tax_table %>% select(-Last), by = "ASV_Id")
    ## Table
    write.table(onlyCycleTable90, file = "Functional-Analyse/Table/Only_Cycles90.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(Only_Cycle90_Sequence90, file = "Functional-Analyse/Table/Only_Cycles_Sequence90.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    #
    # Zone-90 -------------------------------------------------------------------
    ## Mixolimnion
    onlyMixolimnionTable90 <- data_seq_tax %>% filter(Zone90 == "Mixolimnion") %>% select(ASV_Id,TotalMixolimnion,RangInterest,Zone90)
    onlyMixolimnionTable90 <- onlyMixolimnionTable90 %>% filter(TotalMixolimnion > 0.001*sum(onlyMixolimnionTable90$TotalMixolimnion))
    colnames(onlyMixolimnionTable90)[2]  <- "value"
    ## Monimolimnion
    onlyMonimolimnionTable90 <- data_seq_tax %>% filter(Zone90 == "Monimolimnion") %>% select(ASV_Id,TotalMonimolimnion,RangInterest,Zone90)
    onlyMonimolimnionTable90 <- onlyMonimolimnionTable90 %>% filter(TotalMonimolimnion > 0.001*sum(onlyMonimolimnionTable90$TotalMonimolimnion))
    colnames(onlyMonimolimnionTable90)[2]  <- "value"
    ## Bind
    onlyZoneTable90 <- rbind(onlyMixolimnionTable90,onlyMonimolimnionTable90)
    onlyZoneTable90 <- merge(onlyZoneTable90,profile_tax_table %>% select(-Last), by = "ASV_Id")
    ## Table
    write.table(onlyZoneTable90, file = "Functional-Analyse/Table/Only_Zone90.txt", sep = "\t", col.names = TRUE, row.names = FALSE,quote = FALSE)
    write.table(Only_Zone90_Sequence90, file = "Functional-Analyse/Table/Only_Zone_Sequence90.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    #
    # Fraction-90 -------------------------------------------------------------------
    ## Small
    onlySmallTable90 <- data_seq_tax %>% filter(Fraction90 == "Small") %>% select(ASV_Id,TotalSmall,RangInterest,Fraction90)
    onlySmallTable90 <- onlySmallTable90 %>% filter(TotalSmall > 0.001*sum(onlySmallTable90$TotalSmall))
    colnames(onlySmallTable90)[2]  <- "value"
    ## Large
    onlyLargeTable90 <- data_seq_tax %>% filter(Fraction90 == "Large") %>% select(ASV_Id,TotalLarge,RangInterest,Fraction90)
    onlyLargeTable90 <- onlyLargeTable90 %>% filter(TotalLarge > 0.001*sum(onlyLargeTable90$TotalLarge))
    colnames(onlyLargeTable90)[2]  <- "value"
    ## Bind
    onlyFractionTable90 <- rbind(onlySmallTable90,onlyLargeTable90)
    onlyFractionTable90 <- merge(onlyFractionTable90,profile_tax_table %>% select(-Last), by = "ASV_Id")
    ## Table
    write.table(onlyFractionTable90, file = "Functional-Analyse/Table/Only_Fraction90.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(Only_Fraction90_Sequence90, file = "Functional-Analyse/Table/Only_Fraction_Sequence90.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    #
  # Table Total---------------------------------------------------------------------
    # Cycle Total -------------------------------------------------------------------
    ## Day
    totalDayTable <- data_seq_tax %>% select(ASV_Id,TotalDay,RangInterest,Cycle)
    totalDayTable$Cycle <- "Day"
    totalDayTable <- totalDayTable %>% filter(TotalDay > 0.05*sum(totalDayTable$TotalDay))
    colnames(totalDayTable)[2]  <- "value"
    ## Night
    totalNightTable <- data_seq_tax %>% select(ASV_Id,TotalNight,RangInterest,Cycle)
    totalNightTable$Cycle <- "Night"
    totalNightTable <- totalNightTable %>% filter(TotalNight > 0.05*sum(totalNightTable$TotalNight))
    colnames(totalNightTable)[2]  <- "value"
    ## Bind
    totalCycleTable <- rbind(totalDayTable,totalNightTable)
    totalCycleTable <- merge(totalCycleTable,profile_tax_table %>% select(-Last), by = "ASV_Id")
    ## Table
    write.table(totalCycleTable, file = "Functional-Analyse/Table/Total_Cycles.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(Total_Cycle_ASV, file = "Functional-Analyse/Table/Total_Cycles_ASV.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(Total_Cycle_Sequence, file = "Functional-Analyse/Table/Total_Cycles_Sequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    #
    # Zone Total -------------------------------------------------------------------
    ## Mixolimnion
    totalMixolimnionTable <- data_seq_tax %>% select(ASV_Id,TotalMixolimnion,RangInterest,Zone)
    totalMixolimnionTable$Zone <- rep("Mixolimnion",each = nrow(totalMixolimnionTable))
    totalMixolimnionTable <- totalMixolimnionTable %>% filter(TotalMixolimnion > 0.05*sum(totalMixolimnionTable$TotalMixolimnion))
    colnames(totalMixolimnionTable)[2]  <- "value"
    ## Monimolimnion
    totalMonimolimnionTable <- data_seq_tax %>% select(ASV_Id,TotalMonimolimnion,RangInterest,Zone)
    totalMonimolimnionTable$Zone <- rep("Monimolimnion",each = nrow(totalMonimolimnionTable))
    totalMonimolimnionTable <- totalMonimolimnionTable %>% filter(TotalMonimolimnion > 0.05*sum(totalMonimolimnionTable$TotalMonimolimnion))
    colnames(totalMonimolimnionTable)[2]  <- "value"
    ## Bind
    totalZoneTable <- rbind(totalMixolimnionTable,totalMonimolimnionTable)
    totalZoneTable <- merge(totalZoneTable,profile_tax_table %>% select(-Last), by = "ASV_Id")
    ## Table
    write.table(totalZoneTable, file = "Functional-Analyse/Table/Total_Zone.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(Total_Zone_ASV, file = "Functional-Analyse/Table/Total_Zone_ASV.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(Total_Zone_Sequence, file = "Functional-Analyse/Table/Total_Zone_Sequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    #
    # Fraction Total -------------------------------------------------------------------
    ## Small
    totalSmallTable <- data_seq_tax %>% select(ASV_Id,TotalSmall,RangInterest,Fraction)
    totalSmallTable$Fraction <- rep("Small",each = nrow(totalSmallTable))
    totalSmallTable <- totalSmallTable %>% filter(TotalSmall > 0.05*sum(totalSmallTable$TotalSmall))
    colnames(totalSmallTable)[2]  <- "value"
    ## Large
    totalLargeTable <- data_seq_tax %>% select(ASV_Id,TotalLarge,RangInterest,Fraction)
    totalLargeTable$Fraction <- rep("Large",each = nrow(totalLargeTable))
    totalLargeTable <- totalLargeTable %>% filter(TotalLarge > 0.05*sum(totalLargeTable$TotalLarge))
    colnames(totalLargeTable)[2]  <- "value"
    ## Bind
    totalFractionTable <- rbind(totalSmallTable,totalLargeTable)
    totalFractionTable <- merge(totalFractionTable,profile_tax_table %>% select(-Last), by = "ASV_Id")
    ## Table
    write.table(totalFractionTable, file = "Functional-Analyse/Table/Total_Fraction.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(Total_Fraction_ASV, file = "Functional-Analyse/Table/Total_Fraction_ASV.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(Total_Fraction_Sequence, file = "Functional-Analyse/Table/Total_Fraction_Sequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    #
    # Period Total -------------------------------------------------------------------
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
    totalPeriodTable <- merge(totalPeriodTable,profile_tax_table %>% select(-Last), by = "ASV_Id")
    ##
    write.table(totalPeriodTable, file = "Functional-Analyse/Table/Total_Periods.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(Total_Period_ASV, file = "Functional-Analyse/Table/Total_Periods_ASV.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(Total_Period_Sequence, file = "Functional-Analyse/Table/Total_Periods_Sequence.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    #