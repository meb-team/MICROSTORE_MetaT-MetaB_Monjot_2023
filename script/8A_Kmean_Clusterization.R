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
# Create result files
  system("mkdir Functional-Analyse/Group_FUNCTION")
  system("mkdir Functional-Analyse/Hist_Function")
  system("mkdir Functional-Analyse/Tree_Function")
  system("mkdir Functional-Analyse/Composition_Function")
  system("mkdir Functional-Analyse/Table")
  system("mkdir Functional-Analyse/Tree2_Function")

#
# Input trait Table ---------------------------------------------------------
  tableVinput <- read.csv(file = paste("../../../rawdata",input,sep = "/"), sep = "\t")
#
# Prepare table -----------------------------------------------------------
  tableTraits <- tableVinput %>% select(Last,Size.min,Size.max,Cover,Shape,Spicule,Symmetry,Polarity,Colony,Motility,Plast,Ingestion.mode,Biotic.interaction,Resting.stage,Suspected.trophism,IN.OUT)
  tableTraits <- tableTraits %>% distinct(Last,.keep_all=TRUE)
  row.names(tableTraits)<-tableTraits$Last ; tableTraits <- tableTraits %>% select(-Last)
  tableTraits <- as.data.frame(tableTraits)
# Catagorize Trait --------------------------------------------------------
## Size
  tableTraits[,"Size.min"]<- round(log10(tableTraits[,"Size.min"]))
  tableTraits[,"Size.max"]<- round(log10(tableTraits[,"Size.max"]))
## As numeric
  facteur=c("Size.min","Size.max","Cover","Shape","Spicule","Symmetry","Polarity","Colony","Motility","Plast","Ingestion.mode","Biotic.interaction","Resting.stage","Suspected.trophism")
  tableTraits[,all_of(facteur)] <- lapply(tableTraits[,all_of(facteur)], as.character)
  tableTraits <- tableTraits  %>% filter(IN.OUT=="Yes") %>% select(-IN.OUT)
## Remove unnatoted
  tableTraits$Count<-0
  for (i in row.names(tableTraits)) {
    if (sum(is.na(tableTraits[i,])) > 9) {tableTraits[i,"Count"]<-1}}
  tableTraits<- tableTraits %>% filter(Count ==0) %>% select(-Count)

# Gower distance ----------------------------------------------------------
  tableTraits[,all_of(facteur)] <- lapply(tableTraits[,all_of(facteur)], as.factor)
  gower.dissimilarity.mtrx <- daisy(tableTraits,metric = c("gower"))
  dissimilarity.mtrx.csv.content = as.matrix(gower.dissimilarity.mtrx,na.rm=TRUE)
  dissimilarity.mtrx.csv.content[is.na(dissimilarity.mtrx.csv.content)==TRUE]<-0
#
  # PCoA --------------------------------------------------------------------
    set.seed(123)
    res.pcoa<-cmdscale(dissimilarity.mtrx.csv.content,eig = T,k=7)
    p <- res.pcoa$points[,1:7]
    coord <- p ; colnames(coord) <- c("Dim.1","Dim.2","Dim.3","Dim.4","Dim.5","Dim.6","Dim.7")
    coord <- as.data.frame(coord)
    eigtable <- res.pcoa$eig[1:7]
    eigtable <- as.data.frame(eigtable) ; eigtable$row.names <- c("Dim.1","Dim.2","Dim.3","Dim.4","Dim.5","Dim.6","Dim.7")
    Dim1seq <- paste("Dim 1 [",round(eigtable %>% filter(row.names == "Dim.1") %>% select(eigtable),1),"%]", sep = "")
    Dim2seq <- paste("Dim 2 [",round(eigtable %>% filter(row.names == "Dim.2") %>% select(eigtable),1),"%]", sep = "")
    Dim3seq <- paste("Dim 3 [",round(eigtable %>% filter(row.names == "Dim.3") %>% select(eigtable),1),"%]", sep = "")
    Dim4seq <- paste("Dim 4 [",round(eigtable %>% filter(row.names == "Dim.4") %>% select(eigtable),1),"%]", sep = "")
    Dim5seq <- paste("Dim 5 [",round(eigtable %>% filter(row.names == "Dim.5") %>% select(eigtable),1),"%]", sep = "")
    Dim6seq <- paste("Dim 6 [",round(eigtable %>% filter(row.names == "Dim.6") %>% select(eigtable),1),"%]", sep = "")
    Dim7seq <- paste("Dim 7 [",round(eigtable %>% filter(row.names == "Dim.7") %>% select(eigtable),1),"%]", sep = "")
    ## Fusion with trait
    data_seq_trait <- merge(x = coord, y = tableTraits, by = 'row.names')
    row.names(data_seq_trait) <- data_seq_trait$Row.names ; data_seq_trait <- data_seq_trait %>% select(-Row.names)
    ## Spearman ranked analysis
    tableTraits_numic <- as.data.frame(lapply(tableTraits,as.numeric)) ; row.names(tableTraits_numic)<-row.names(coord)
    correlation_table <-corr.test(coord,tableTraits_numic,method = "spearman")
    correlation_table_r <- correlation_table[["r"]]
    correlation_table_p <- correlation_table[["p"]]
    correlation_table <- melt(correlation_table_r)
    correlation_table$pvalue<-melt(correlation_table_p) %>% select(value)
    correlation_table[,"value"][correlation_table[,"pvalue"] > 0.05] <- 0
    correlation_table[,"value"][correlation_table[,"value"] < 0.3 & correlation_table[,"value"] > -0.3] <- 0
    correlation_table[,1:2] <- data.frame(lapply(correlation_table[,1:2], function(x) {gsub("\\.", " ", x)}))
    ggplot(data = correlation_table, aes(x=Var1, y=Var2, fill=abs(value))) + 
      scale_fill_gradient2(low="#245d72", mid="#edededFF",high="#a50326", space ="Lab" ) +
      geom_tile(color="white") + labs(x="",y="",fill="Absolute value of the\nspearman coefficient") + theme_unique_darkcris()
    ggsave(paste0("Functional-Analyse/Group_FUNCTION/","Correlation",".svg"), device = "svg", width = 5, height = 5)
#
  # Plot --------------------------------------------------------------------
    myplots <- list()
    i<-0
    for (factor in all_of(facteur)) {
      # Dim1 and 2 -------------------------------------------------------------------
      i<-i+1
      print(noquote(factor))
      data <- data_seq_trait %>% filter(is.na(get(factor))==FALSE)
      a <- ggplot(data, aes(y = `Dim.2`, x = `Dim.1`, color = data[,factor])) + geom_point(size = 2) +
        geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
        geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
        scale_linetype_manual(values = c(0,2,2)) +
        theme(axis.title = element_text(face="bold", size=12), 
              axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
              title = element_text(face="bold", size=14),
              legend.title = element_text(face="bold"),
              legend.position = "right",
              legend.text = element_text(size=10)) +
        labs(x=Dim1seq,y=Dim2seq,color = factor, alpha = factor, linetype = factor) + theme(legend.position = "none")
      #
      # Dim3 and 4 -------------------------------------------------------------------
      i<-i+1
      print(noquote(factor))
      data <- data_seq_trait %>% filter(is.na(get(factor))==FALSE)
      b <- ggplot(data, aes(y = `Dim.4`, x = `Dim.3`, color = data[,factor])) + geom_point(size = 2) +
        geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
        geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
        scale_linetype_manual(values = c(0,2,2)) +
        theme(axis.title = element_text(face="bold", size=12), 
              axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
              title = element_text(face="bold", size=14),
              legend.title = element_text(face="bold"),
              legend.position = "right",
              legend.text = element_text(size=10)) +
        labs(x=Dim3seq,y=Dim4seq,color = factor, alpha = factor, linetype = factor) + theme(legend.position = "none")
      #
      # Dim5 and 6 -------------------------------------------------------------------
      i<-i+1
      print(noquote(factor))
      data <- data_seq_trait %>% filter(is.na(get(factor))==FALSE)
      c <- ggplot(data, aes(y = `Dim.6`, x = `Dim.5`, color = data[,factor])) + geom_point(size = 2) +
        geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
        geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
        scale_linetype_manual(values = c(0,2,2)) +
        theme(axis.title = element_text(face="bold", size=12), 
              axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
              title = element_text(face="bold", size=14),
              legend.title = element_text(face="bold"),
              legend.position = "right",
              legend.text = element_text(size=10)) +
        labs(x=Dim5seq,y=Dim6seq,color = factor, alpha = factor, linetype = factor)
      legend<-get_legend(c)
      c<-c + theme(legend.position = "none")
      #
      # Dim7 and 1 -------------------------------------------------------------------
      i<-i+1
      print(noquote(factor))
      data <- data_seq_trait %>% filter(is.na(get(factor))==FALSE)
      d <- ggplot(data, aes(y = `Dim.1`, x = `Dim.7`, color = data[,factor])) + geom_point(size = 2) +
        geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
        geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
        scale_linetype_manual(values = c(0,2,2)) +
        theme(axis.title = element_text(face="bold", size=12), 
              axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
              title = element_text(face="bold", size=14),
              legend.title = element_text(face="bold"),
              legend.position = "right",
              legend.text = element_text(size=10)) +
        labs(x=Dim7seq,y=Dim1seq,color = factor, alpha = factor, linetype = factor) + theme(legend.position = "none")
      d<-d + theme(legend.position = "none")
      
      #svglite(paste0("Group_FUNCTION/","PCoA_",factor,".svg"),width = 12.00,height = 4.00)
      assign(paste0("All_plot_",factor),plot_grid(plotlist=list(a,b,c,d,legend),ncol = 5, nrow = 1))
      #All_plot <- plot_grid(plotlist=list(a,b,c,d,legend),ncol = 5, nrow = 1)
      #print(All_plot)
      #dev.off()
    }
    list_data_Grp_all_ <- ls()
    list_data_Grp_all_<- grep(list_data_Grp_all_,pattern="All_plot_",value = TRUE)
    list_data_Grp_all_<- lapply(list_data_Grp_all_,get)
    All_plot <- plot_grid(plotlist=list_data_Grp_all_,ncol = 2, nrow = 7, align = "hv")
    svglite(paste0("Functional-Analyse/Group_FUNCTION/","PCoA_All.svg"),width = 24.00,height = 20.00)
    print(All_plot)
    dev.off()
#
# Kmean Cluster ----------------------------------------------
    set.seed(123) # 123
    spe.KM.cascade_bis <- cascadeKM(coord %>% select(Dim.1,Dim.2), inf.gr=3, sup.gr=12, iter=1000, criterion="ssi")
    KM.cascade.data <- spe.KM.cascade_bis
    svglite(paste0("Functional-Analyse/Group_FUNCTION/","Kmeans-ssi",".svg"),width = 10.00,height = 6.00)
    plot(KM.cascade.data, sortg=TRUE)
    dev.off()
    data_Group<-as.data.frame(KM.cascade.data$partition)
    nbofGroup<-as.data.frame(KM.cascade.data[["results"]]["ssi",]) ; colnames(nbofGroup) <- "ssi" ; nbofGroup$Group <- row.names(nbofGroup)
    nbofGroupi<- as.character(nbofGroup %>% filter(ssi==max(nbofGroup$ssi)) %>% select(Group))
    data_Group_seq <- merge(x = data_Group %>% select(all_of(nbofGroupi)), y = data_seq_trait, by = 'row.names')
    data_Group_seq$nbGroup <- data_Group_seq[,nbofGroupi]
    ###
    Ngrp <- as.numeric(length(unique(data_Group[,nbofGroupi])))
#
  # Plot Kmean-ssi PCoA --------------------------------------------------------------------
    data <- data_Group_seq
    for (i in row.names(data)) {data[i,"nbGroup"] <- paste("Cluster",data[i,"nbGroup"])}
    #colnames(data)<-c("Row.names","nbGroup","Dim.1","Dim.2","Dim.3","Dim.4","Dim.5","Dim.6","Dim.7","Size.min","Size.max","Cover","Shape","Spicule","Symmetry","Polarity","Colony","Motility","Plast.origin","Ingestion.mode","Symbiontic","Resting.stage","Suspected.trophism")
    #data$nbGroup <- as.character(data$nbGroup)
    a <- ggplot(data, aes(y = `Dim.2`, x = `Dim.1`, color = nbGroup)) + geom_point(size = 2) +
      geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
      geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
      stat_ellipse(aes(linetype = nbGroup, alpha = nbGroup),geom = "polygon",type = "norm") +
      scale_linetype_manual(values = rep(2,Ngrp)) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      scale_alpha_manual(values = rep(0,Ngrp)) +
      scale_color_manual(values = c("#F8766D","#0099B4CC","#00A5FF","#AD002ACC","#EEA236CC","#00468BCC","#c2d1fcCC","#64908aCC",'lightgrey',"lightblue","#d7c11c","darkblue")) +
      labs(x=Dim1seq,y=Dim2seq,color = "Clusters", alpha = "Clusters", linetype = "Clusters") + theme(legend.position = "none")
    print(a)
    b <- ggplot(data, aes(y = `Dim.4`, x = `Dim.3`, color = nbGroup)) + geom_point(size = 2) +
      geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
      geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
      stat_ellipse(aes(linetype = nbGroup, alpha = nbGroup),geom = "polygon",type = "norm") +
      scale_linetype_manual(values = rep(2,Ngrp)) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      scale_alpha_manual(values = rep(0,Ngrp)) +
      scale_color_manual(values = c("#F8766D","#0099B4CC","#00A5FF","#AD002ACC","#EEA236CC","#00468BCC","#c2d1fcCC","#64908aCC",'lightgrey',"lightblue","#d7c11c","darkblue")) +
      labs(x=Dim3seq,y=Dim4seq,color = "Clusters", alpha = "Clusters", linetype = "Clusters") + theme(legend.position = "none")
    print(b)
    c <- ggplot(data, aes(y = `Dim.6`, x = `Dim.5`, color = nbGroup)) + geom_point(size = 2) +
      geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
      geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
      stat_ellipse(aes(linetype = nbGroup, alpha = nbGroup),geom = "polygon",type = "norm") +
      scale_linetype_manual(values = rep(2,Ngrp)) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      scale_alpha_manual(values = rep(0,Ngrp)) +
      scale_color_manual(values = c("#F8766D","#0099B4CC","#00A5FF","#AD002ACC","#EEA236CC","#00468BCC","#c2d1fcCC","#64908aCC",'lightgrey',"lightblue","#d7c11c","darkblue")) +
      labs(x=Dim5seq,y=Dim6seq,color = "Clusters", alpha = "Clusters", linetype = "Clusters")
    print(c)
    legend<-get_legend(c)
    d <- ggplot(data, aes(y = `Dim.1`, x = `Dim.7`, color = nbGroup)) + geom_point(size = 2) +
      geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
      geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
      stat_ellipse(aes(linetype = nbGroup, alpha = nbGroup),geom = "polygon",type = "norm") +
      scale_linetype_manual(values = rep(2,Ngrp)) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      scale_alpha_manual(values = rep(0,Ngrp)) +
      scale_color_manual(values = c("#F8766D","#0099B4CC","#00A5FF","#AD002ACC","#EEA236CC","#00468BCC","#c2d1fcCC","#64908aCC",'lightgrey',"lightblue","#d7c11c","darkblue")) +
      labs(x=Dim7seq,y=Dim1seq,color = "Clusters", alpha = "Clusters", linetype = "Clusters")
    print(d)
    d <- d + theme(legend.position = "none")
    c <- c + theme(legend.position = "none")
    ## Plot  
    svglite(paste0("Functional-Analyse/Group_FUNCTION/","PCoA_Kmean-ssi",".svg"),width = 12.00,height = 4.00)
    All_plot <- plot_grid(plotlist=list(a,b,c,d,legend),ncol = 5, nrow = 1)
    print(All_plot)
    dev.off()
    #
    svglite(paste0("Functional-Analyse/Group_FUNCTION/","PCoA_Kmean-ssi-DIM12",".svg"),width = 5.00,height = 4.00)
    All_plot <- plot_grid(plotlist=list(a,legend),ncol = 2, nrow = 1, rel_widths = c(4,2))
    print(All_plot)
    dev.off()
    #
    # GB16/03/2023
    for (i in colnames(data_Group)) {
      data_Group_seq_ijk <- merge(x = data_Group %>% select(all_of(i)), y = data_seq_trait, by = 'row.names')
      # exchange the name of the groups to get the same annotation as before ()
      data_Group_seq_ijk$nbGroup <- data_Group_seq_ijk[,i]
      ###
      Ngrp <- as.numeric(length(unique(data_Group[,i])))
      #
      dataijk <- data_Group_seq_ijk
      for (j in row.names(dataijk)) {dataijk[j,"nbGroup"] <- paste("Cluster",dataijk[j,"nbGroup"])}
      a <- ggplot(dataijk, aes(y = `Dim.2`, x = `Dim.1`, color = nbGroup)) + geom_point(size = 2) +
        geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
        geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
        stat_ellipse(aes(linetype = nbGroup, alpha = nbGroup),geom = "polygon",type = "norm") +
        scale_linetype_manual(values = rep(2,Ngrp)) +
        theme(axis.title = element_text(face="bold", size=12), 
              axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
              title = element_text(face="bold", size=14),
              legend.title = element_text(face="bold"),
              legend.position = "right",
              legend.text = element_text(size=10)) +
        scale_alpha_manual(values = rep(0,Ngrp)) +
        scale_color_manual(values = c("#F8766D","#0099B4CC","#00A5FF","#AD002ACC","#EEA236CC","#00468BCC","#c2d1fcCC","#64908aCC",'lightgrey',"lightblue","#d7c11c","darkblue","lightgreen")) +
        labs(x=Dim1seq,y=Dim2seq,color = "Functional groups", alpha = "Functional groups", linetype = "Functional groups") + theme(legend.position = "none")
      svglite(paste0("Functional-Analyse/Group_FUNCTION/","PCoA_Kmean-ssi-",i,"cluster.svg"),width = 4.00,height = 4.00)
      print(a)
      dev.off()
    }
#
  # Group definition --------------------------------------------------------
    z<-as.numeric(str_split(nbofGroupi," ")[[1]][1])
    for (i in 1:z) { print(i)
      data_Grp <- data %>% filter(nbGroup==paste("Cluster",i))
      data_Grp_Size.min<- data_Grp %>% group_by(Value=Size.min) %>% summarise(n = n()) ; data_Grp_Size.min$Factor="Size min"
      data_Grp_Size.max<- data_Grp %>% group_by(Value=Size.max) %>% summarise(n = n()) ; data_Grp_Size.max$Factor="Size max"
      data_Grp_Cover<- data_Grp %>% group_by(Value=Cover) %>% summarise(n = n()) ; data_Grp_Cover$Factor="Cover"
      data_Grp_Shape<- data_Grp %>% group_by(Value=Shape) %>% summarise(n = n()) ; data_Grp_Shape$Factor="Shape"
      data_Grp_Colony<- data_Grp %>% group_by(Value=Colony) %>% summarise(n = n()) ; data_Grp_Colony$Factor="Colony"
      data_Grp_Spicule<- data_Grp %>% group_by(Value=Spicule) %>% summarise(n = n()) ; data_Grp_Spicule$Factor="Spicule"
      data_Grp_Symmetry<- data_Grp %>% group_by(Value=Symmetry) %>% summarise(n = n()) ; data_Grp_Symmetry$Factor="Symmetry"
      data_Grp_Polarity<- data_Grp %>% group_by(Value=Polarity) %>% summarise(n = n()) ; data_Grp_Polarity$Factor="Polarity"
      data_Grp_Plast.origin<- data_Grp %>% group_by(Value=Plast) %>% summarise(n = n()) ; data_Grp_Plast.origin$Factor="Plast"
      data_Grp_Symbiotic<- data_Grp %>% group_by(Value=Biotic.interaction) %>% summarise(n = n()) ; data_Grp_Symbiotic$Factor="Biotic interaction"
      data_Grp_Resting.stage<- data_Grp %>% group_by(Value=Resting.stage) %>% summarise(n = n()) ; data_Grp_Resting.stage$Factor="Resting stage"
      data_Grp_Suspected.trophism<- data_Grp %>% group_by(Value=Suspected.trophism) %>% summarise(n = n()) ; data_Grp_Suspected.trophism$Factor="Suspected trophism"
      data_Grp_Motility<- data_Grp %>% group_by(Value=Motility) %>% summarise(n = n()) ; data_Grp_Motility$Factor="Motility"
      data_Grp_Ingestion.mode<- data_Grp %>% group_by(Value=Ingestion.mode) %>% summarise(n = n()) ; data_Grp_Ingestion.mode$Factor="Ingestion mode"
      assign(paste0("data_Grp_all_",i),rbind(data_Grp_Size.min,data_Grp_Size.max,data_Grp_Cover,data_Grp_Shape,data_Grp_Colony,data_Grp_Spicule,data_Grp_Symmetry,data_Grp_Polarity,data_Grp_Plast.origin,data_Grp_Symbiotic,data_Grp_Resting.stage,data_Grp_Suspected.trophism,data_Grp_Motility,data_Grp_Ingestion.mode))
    }
    #
  # PlotMultiplGRP ----------------------------------------------------------
    list_data_Grp_all_ <- ls()
    list_data_Grp_all_<- grep(list_data_Grp_all_,pattern="^data_Grp_all_",value = TRUE)
    k<-0
    for (i in list_data_Grp_all_){
      k<-k+1
      z<-as.numeric(str_split(i,"_")[[1]][4])
      assign(i,get(i) %>% filter(is.na(Value)!=TRUE) %>% add_column(Group=as.character(z)))
      if (k==1) {data_Grp_all<-get(i)}
      if (k>1) {data_Grp_all <- rbind(data_Grp_all,get(i))}
    }
    
    ## Plot
    for (i in row.names(data_Grp_all)) {data_Grp_all[i,"Group"] <- paste("Cluster",data_Grp_all[i,"Group"])}
    svglite(paste0("Functional-Analyse/Group_FUNCTION/","Def_group_by_factor",".svg"),width = 10.00,height = 10.00)
    histGrpall <- ggplot(data_Grp_all, mapping = aes(y=Value, x = n, fill = Factor, group = Factor), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity") + 
      scale_fill_paletteer_d("ggthemes::manyeys") + scale_y_discrete(position = "right") +
      facet_grid(Factor~Group,scales="free",space="free_y",switch = "y") + theme_unique_darkcris() + theme(strip.background = element_rect(color="black", fill="white", size=1, linetype="blank"))
    legendGRP <- get_legend(histGrpall)
    histGrpall <- histGrpall + theme(legend.position = "none",strip.text.y.left = element_text(angle = 0, hjust = 1)) #+ labs(x="Cycles",y="ASVs (%)") 
    print(histGrpall)
    dev.off()
    save.image("Functional-Analyse/Group_FUNCTION/RHistory-Traitbased.RData")
    #
# Code-Table_Last.Group ---------------------------------------------------
  ## Code Table
    cluster_name <- data %>% select(nbGroup) %>% distinct(nbGroup) %>% arrange(nbGroup)
    cluster_name$name <- ""
  ## Write table and print fig
    svglite(paste0("../../../rawdata/Def_group_by_factor.svg"),width = 10.00,height = 10.00)
    print(histGrpall)
    dev.off()
    #
    write.table(cluster_name,"../../../rawdata/cluster_name.tsv",row.names = FALSE,sep="\t",quote = FALSE)
  ## Save Rdata
    save.image("Functional-Analyse/8A_Kmean_Clusterization.RData")
#