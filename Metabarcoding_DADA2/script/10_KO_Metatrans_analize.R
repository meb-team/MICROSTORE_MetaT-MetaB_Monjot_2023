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
input <- "main_table.mapping.unique.read_per_kb.noHuman.noConta.noMetazoa.annot.tsv"
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
#
# Import package -----------------------------------------------------------
pkg <- c("ggplot2","dplyr","tidyr","cowplot","FactoMineR","factoextra","reshape2","varhandle","ggrepel","ggpubr","ggsci","scales","hrbrthemes","data.table","svglite","treemap", "VennDiagram","stringr","paletteer","elementalist","gtools","vegan","treemapify","gplots","rgl","car","scatterplot3d")
lapply(pkg, require, character.only = TRUE)

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
dir.create("Table")
dir.create("Stat")
dir.create("Tree_Function")
dir.create("Polar")
dir.create("Hist")
dir.create("Balance")
dir.create("Hist_Meta")
dir.create("Quaterneryplot")
dir.create("Venn")
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
get_colors <- function(groups, group.col = palette()){
  groups <- as.factor(groups)
  ngrps <- length(levels(groups))
  if(ngrps > length(group.col)) 
    group.col <- rep(group.col, ngrps)
  color <- group.col[as.numeric(groups)]
  names(color) <- as.vector(groups)
  return(color)
}
rgl_init <- function(new.device = FALSE, bg = "white", width = 640) { 
  if( new.device | rgl.cur() == 0 ) {
    rgl.open()
    par3d(windowRect = 50 + c( 0, 0, width, width ) )
    rgl.bg(color = bg )
  }
  rgl.clear(type = c("shapes", "bboxdeco"))
  rgl.viewpoint(theta = 15, phi = 20, zoom = 0.7)
}
html_to_pdf <- function(html_file, pdf_file) {
  cmd <- sprintf("pandoc %s -t html -o %s", html_file, pdf_file)
  system(cmd)
}
# Input Tables ---------------------------------------------------------
tableVinput <- fread(file = paste("../../../../rawdata/",input,sep = "/"), sep = "\t")
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
infdataini <- read.table(file = "../../../../rawdata/data-inf-metaT.txt", sep = ";", header = TRUE)
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
samples_df <- infdataini[,all_of(col)]
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
# Normalization TPM -------------------------------------------------------
  tableVinput_abundance_normalize <- setDF(tableVinput_abundance,rownames = tableVinput_abundance$unigene) %>% dplyr::select(-unigene, -PfamDom, -ko, -ko_BestHitNotSignif, -GOterms)
  tableVinput_abundance_normalize["scalingFactor",] <- colSums(tableVinput_abundance_normalize) / 1000000
  for (i in colnames(tableVinput_abundance_normalize)) {
    tableVinput_abundance_normalize[1:nrow(tableVinput_abundance_normalize)-1,i] <- tableVinput_abundance_normalize[1:nrow(tableVinput_abundance_normalize)-1,i] / tableVinput_abundance_normalize[nrow(tableVinput_abundance_normalize),i] }
  tableVinput_abundance_normalize <- tableVinput_abundance_normalize[1:nrow(tableVinput_abundance_normalize)-1,]
#
# Stat 2 and first filter Somme = 0 --------------------------------------------------------------------
  tableVinput_abundance_normalize$Somme <- rowSums(tableVinput_abundance_normalize)
  nrow(tableVinput_abundance_normalize)
  nrow(tableVinput_abundance_normalize %>% filter(Somme != 0))
  tableVinput_abundance_normalize0 <- tableVinput_abundance_normalize %>% dplyr::filter(Somme == 0)
  tableVinput_abundance_normalize <- tableVinput_abundance_normalize %>% dplyr::filter(Somme != 0) %>% dplyr::select(-Somme)
  
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
    oldcol <- colnames(tableVinput_abundance_normalize)
    w <- 0
    for (i in oldcol) { f <- as.character(MetaT_Summary %>% dplyr::filter(Ref..collab == i) %>% dplyr::select(Echantillon))
      if (f != "character(0)") { 
        if (w == 0) { newcol <- f 
          selectcol <- i }
        else { newcol <- c(all_of(newcol),f)
          selectcol <- c(all_of(selectcol),i)}}
    w <- w+1 }
    tableVinput_abundance_normalize <- tableVinput_abundance_normalize %>% dplyr::select(all_of(selectcol))
    colnames(tableVinput_abundance_normalize) <- newcol
#
# Stat 3 and second filter Somme = 0 --------------------------------------------------------------------
    tableVinput_abundance_normalize$Somme <- rowSums(tableVinput_abundance_normalize)
    nrow(tableVinput_abundance_normalize)
    nrow(tableVinput_abundance_normalize %>% filter(Somme != 0))
    tableVinput_abundance_normalize0 <- tableVinput_abundance_normalize %>% dplyr::filter(Somme == 0)
    tableVinput_abundance_normalize <- tableVinput_abundance_normalize %>% dplyr::filter(Somme != 0) %>% dplyr::select(-Somme)
#    
# Merge Tables ------------------------------------------------------------
    tableVinput_abundance_normalize_annots <- merge(tableVinput_abundance_normalize,tableVinput_abundance %>% dplyr::select(unigene,PfamDom,ko,ko_BestHitNotSignif,GOterms),by.x="row.names",by.y="unigene")
    tableVinput_abundance_normalize_annots_tax <- merge(tableVinput_abundance_normalize_annots,Tax_table, by.x = "Row.names", by.y = "Unigene")
    tableVinput_abundance_normalize_annots_tax_funct <- 
    rm(tableVinput_abundance_normalize_annots)
#    
# Taxonomic Annotation ---------------------------------------------------------
  # Process Taxonomic data ------------------------------------------------------------
    # Prepare Taxonomy
    uniq_tax <- as.data.frame(unique(tableVinput_abundance_normalize_annots_tax$Values)) ; colnames(uniq_tax) <- "uniq"
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
    # Merge in tableVinput_abundance_normalize_annots_tax
    tableVinput_abundance_normalize_annots_tax <- merge(tableVinput_abundance_normalize_annots_tax, uniq_tax %>% dplyr::select("uniq", "PR2_taxonomy","PR2_taxonomy_complete"), by.x = "Values", by.y = "uniq")    
#
  # General Tax Stat --------------------------------------------------------------------
    data_seq_summarize <- tableVinput_abundance_normalize_annots_tax %>% dplyr::select(Values,Row.names,PR2_taxonomy_complete)
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
    if (inputmode == FALSE) {save.image("R_TaxonomyOK_image.RData")}
    
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
    tableVinput_abundance_normalize_annots_tax_Funct <- dplyr::left_join(tableVinput_abundance_normalize_annots_tax,data_seq_funct %>% dplyr::select(uniq,nbGroup,Last),by = c("Values"= "uniq"))
    tableVinput_abundance_normalize_annots_tax_Funct <- separate(tableVinput_abundance_normalize_annots_tax_Funct, PR2_taxonomy_complete, c("Domain","Supergroup","Division","Class","Order","Family","Genus","Species"), sep =";")
    interest <- tableVinput_abundance_normalize_annots_tax_Funct %>% filter(grepl(ko_BestHitNotSignif,pattern=",")) %>% dplyr::select(Row.names)
    tableVinput_abundance_normalize_annots_tax_Funct[,"ko_BestHitNotSignif"][tableVinput_abundance_normalize_annots_tax_Funct[,"Row.names"] %in% interest$Row.names] <- NA
    tableVinput_abundance_normalize_annots_tax_Funct[,"ko_BestHitNotSignif"][is.na(tableVinput_abundance_normalize_annots_tax_Funct[,"ko"]) == FALSE] <- NA
    tableVinput_abundance_normalize_annots_tax_Funct <- tableVinput_abundance_normalize_annots_tax_Funct %>% tidyr::unite(., col = "Kegg.Onthology",  ko, ko_BestHitNotSignif, na.rm=TRUE, sep = ",",remove=FALSE)
    #tableVinput_abundance_normalize_annots_tax_Funct <- separate_rows(tableVinput_abundance_normalize_annots_tax_Funct,"Kegg.Onthology",sep = ",") # Duplicate row by KO 
    tableVinput_abundance_normalize_annots_tax_Funct[tableVinput_abundance_normalize_annots_tax_Funct==""] <- NA
    tableVinput_abundance_normalize_annots_tax_Funct[,"nbGroup"][tableVinput_abundance_normalize_annots_tax_Funct[,"nbGroup"]=="Unassigned"] <- NA
  #
    nrow(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::filter(nchar(Kegg.Onthology)>6))
  ## General Grp Stat
    data_seq_summarize <- tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(Values,Row.names,ko,PfamDom,ko_BestHitNotSignif,Kegg.Onthology,GOterms,nbGroup)
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
    Fun_table_TF <- tableVinput_abundance_normalize_annots_tax_Funct
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
      geom_label(aes(y = value + 0.08*max(value),label = paste(value,"\n","RPKM: ", round(value*100/sum(Fun_table_TF$Somme),1)," %",sep =""), color = variable, alpha = variable),size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
      scale_color_manual(values = c("#0000FF00","black")) +
      scale_alpha_manual(values = c(0,0.5)) +
      scale_linetype_manual(values = c("blank","dashed")) +
      guides(fill=guide_legend("Annotation"), color = "none", linetype = "none") +
      labs(x="Annotation type",y="Abundance (RPKM%)") +
      theme(legend.title = element_text(face="bold"),
            axis.title = element_text(color = "black", face = "bold"))
    print(fun_plot)
    ggsave("Stat/Stat-Functional-Abundance.svg", device = "svg", width = 9, height = 4)
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


# Stat through conditions -------------------------------------------------
  ## General Grp Stat
  for (i in c("SmallXMonimolimnion","LargeXMonimolimnion","SmallXMixolimnion","LargeXMixolimnion")) {
    data_seq_summarize <- tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(all_of(get(i)),Values,Row.names,ko,PfamDom,ko_BestHitNotSignif,Kegg.Onthology,GOterms,nbGroup)
    data_seq_summarize$Somme <- rowSums(data_seq_summarize %>% dplyr::select(all_of(get(i))))
    data_seq_summarize <- data_seq_summarize %>% filter(Somme != 0) %>% dplyr::select(Values,Row.names,ko,PfamDom,ko_BestHitNotSignif,Kegg.Onthology,GOterms,nbGroup)
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
    ggsave(paste0("Stat/Stat-Functional-",i,".svg"), device = "svg", width = 9, height = 4)
    ## General Grp Stat Abundance
    #
    Fun_table_TF <- tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(all_of(get(i)),Values,Row.names,ko,PfamDom,ko_BestHitNotSignif,Kegg.Onthology,GOterms,nbGroup)
    Fun_table_TF$Somme <- rowSums(Fun_table_TF %>% dplyr::select(all_of(get(i))))
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
      geom_label(aes(y = value + 0.08*max(value),label = paste(value,"\n","RPKM: ", round(value*100/sum(Fun_table_TF$Somme),1)," %",sep =""), color = variable, alpha = variable),size = 3,show.legend = FALSE, fill = "#3B3B3B99") +
      scale_color_manual(values = c("#0000FF00","black")) +
      scale_alpha_manual(values = c(0,0.5)) +
      scale_linetype_manual(values = c("blank","dashed")) +
      guides(fill=guide_legend("Annotation"), color = "none", linetype = "none") +
      labs(x="Annotation type",y="Abundance (RPKM%)") +
      theme(legend.title = element_text(face="bold"),
            axis.title = element_text(color = "black", face = "bold"))
    print(fun_plot)
    ggsave(paste0("Stat/Stat-Functional-Abundance-",i,".svg"), device = "svg", width = 9, height = 4)
  }
# dataframe condition -----------------------------------------------------
    row.names(tableVinput_abundance_normalize_annots_tax_Funct) <- tableVinput_abundance_normalize_annots_tax_Funct$Row.names
  #Cycle
    ## Day
    Day_seq <- as.data.frame(rowSums(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(all_of(CycleDay))))
    colnames(Day_seq) <- "TotalDay"
    ## Night
    Night_seq <- as.data.frame(rowSums(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(all_of(CycleNight))))
    colnames(Night_seq) <- "TotalNight"
    ## Merge
    Cycle_seq <- merge(x = Day_seq,y = Night_seq, by = "row.names")
  # Zone
    ## Mixolimnion
    Mixolimnion_seq <- as.data.frame(rowSums(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(all_of(ZoneMixolimnion))))
    colnames(Mixolimnion_seq) <- "TotalMixolimnion"
    ## Monimolimnion
    Monimolimnion_seq <- as.data.frame(rowSums(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(all_of(ZoneMonimolimnion))))
    colnames(Monimolimnion_seq) <- "TotalMonimolimnion"
    ## Merge
    Zone_seq <- merge(x = Mixolimnion_seq,y = Monimolimnion_seq, by = "row.names")
  # Fraction
    ## Small
    Small_seq <- as.data.frame(rowSums(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(all_of(FractionSmall))))
    colnames(Small_seq) <- "TotalSmall"
    ## Large
    Large_seq <- as.data.frame(rowSums(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(all_of(FractionLarge))))
    colnames(Large_seq) <- "TotalLarge"
    ## Merge
    Fraction_seq <- merge(x = Small_seq,y = Large_seq, by = "row.names")
  # FractionXZone
  ## Abundance
    ## SmallXMonimolimnion
    SmallXMonimolimnion_seq <- as.data.frame(rowSums(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(all_of(SmallXMonimolimnion))))
    colnames(SmallXMonimolimnion_seq) <- "SmallXMonimolimnion"
    ## SmallXMixolimnion
    SmallXMixolimnion_seq <- as.data.frame(rowSums(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(all_of(SmallXMixolimnion))))
    colnames(SmallXMixolimnion_seq) <- "SmallXMixolimnion"
    ## LargeXMonimolimnion
    LargeXMonimolimnion_seq <- as.data.frame(rowSums(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(all_of(LargeXMonimolimnion))))
    colnames(LargeXMonimolimnion_seq) <- "LargeXMonimolimnion"
    ## LargeXMixolimnion
    LargeXMixolimnion_seq <- as.data.frame(rowSums(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(all_of(LargeXMixolimnion))))
    colnames(LargeXMixolimnion_seq) <- "LargeXMixolimnion"
    ## Merge
    FractionXZone_seq <- merge(x = SmallXMonimolimnion_seq,y = SmallXMixolimnion_seq, by = "row.names")
    FractionXZone_seq <- merge(x = FractionXZone_seq,y = LargeXMonimolimnion_seq, by.x = "Row.names", by.y = "row.names")
    FractionXZone_seq <- merge(x = FractionXZone_seq,y = LargeXMixolimnion_seq, by.x = "Row.names", by.y = "row.names")
  # FractionXZoneXCycle
    ## Abundance
    ## SmallXMonimolimnionXDay
    SmallXMonimolimnionXDay_seq <- as.data.frame(rowSums(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(all_of(SmallXMonimolimnionXDay))))
    colnames(SmallXMonimolimnionXDay_seq) <- "SmallXMonimolimnionXDay"
    ## SmallXMixolimnionXDay
    SmallXMixolimnionXDay_seq <- as.data.frame(rowSums(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(all_of(SmallXMixolimnionXDay))))
    colnames(SmallXMixolimnionXDay_seq) <- "SmallXMixolimnionXDay"
    ## LargeXMonimolimnionXDay
    LargeXMonimolimnionXDay_seq <- as.data.frame(rowSums(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(all_of(LargeXMonimolimnionXDay))))
    colnames(LargeXMonimolimnionXDay_seq) <- "LargeXMonimolimnionXDay"
    ## LargeXMixolimnionXDay
    LargeXMixolimnionXDay_seq <- as.data.frame(rowSums(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(all_of(LargeXMixolimnionXDay))))
    colnames(LargeXMixolimnionXDay_seq) <- "LargeXMixolimnionXDay"
    
    ## SmallXMonimolimnionXNight
    SmallXMonimolimnionXNight_seq <- as.data.frame(rowSums(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(all_of(SmallXMonimolimnionXNight))))
    colnames(SmallXMonimolimnionXNight_seq) <- "SmallXMonimolimnionXNight"
    ## SmallXMixolimnionXNight
    SmallXMixolimnionXNight_seq <- as.data.frame(rowSums(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(all_of(SmallXMixolimnionXNight))))
    colnames(SmallXMixolimnionXNight_seq) <- "SmallXMixolimnionXNight"
    ## LargeXMonimolimnionXNight
    LargeXMonimolimnionXNight_seq <- as.data.frame(rowSums(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(all_of(LargeXMonimolimnionXNight))))
    colnames(LargeXMonimolimnionXNight_seq) <- "LargeXMonimolimnionXNight"
    ## LargeXMixolimnionXNight
    LargeXMixolimnionXNight_seq <- as.data.frame(rowSums(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(all_of(LargeXMixolimnionXNight))))
    colnames(LargeXMixolimnionXNight_seq) <- "LargeXMixolimnionXNight"
    
    ## Merge
    DayXSmallXZone_seq <- merge(x = SmallXMonimolimnionXDay_seq,y = SmallXMixolimnionXDay_seq, by = "row.names")
    DayXLargeXZone_seq <- merge(x = LargeXMonimolimnionXDay_seq,y = LargeXMixolimnionXDay_seq, by = "row.names")
    DayXFractionXZone_seq <- merge(x = DayXSmallXZone_seq,y = DayXLargeXZone_seq, by.x = "Row.names", by.y = "Row.names")
    NightXSmallXZone_seq <- merge(x = SmallXMonimolimnionXNight_seq,y = SmallXMixolimnionXNight_seq, by = "row.names")
    NightXLargeXZone_seq <- merge(x = LargeXMonimolimnionXNight_seq,y = LargeXMixolimnionXNight_seq, by = "row.names")
    NightXFractionXZone_seq <- merge(x = NightXSmallXZone_seq,y = NightXLargeXZone_seq, by.x = "Row.names", by.y = "Row.names")
    #
    CycleXFractionXZone_seq <- merge(x = DayXFractionXZone_seq,y = NightXFractionXZone_seq, by.x = "Row.names", by.y = "Row.names")
  #
  ## unigene P/A
    # FractionXZone P/A
    tableVinput_richness_normalize_annots_tax_Funct <- tableVinput_abundance_normalize_annots_tax_Funct
    tableVinput_richness_normalize_annots_tax_Funct[,all_of(MetaT_Summary$Echantillon)][tableVinput_richness_normalize_annots_tax_Funct[,all_of(MetaT_Summary$Echantillon)]>0] <- 1
    ## SmallXMonimolimnion
    SmallXMonimolimnion_pa <- as.data.frame(rowSums(tableVinput_richness_normalize_annots_tax_Funct %>% dplyr::select(all_of(SmallXMonimolimnion))))
    colnames(SmallXMonimolimnion_pa) <- "SmallXMonimolimnion"
    ## SmallXMixolimnion
    SmallXMixolimnion_pa <- as.data.frame(rowSums(tableVinput_richness_normalize_annots_tax_Funct %>% dplyr::select(all_of(SmallXMixolimnion))))
    colnames(SmallXMixolimnion_pa) <- "SmallXMixolimnion"
    ## LargeXMonimolimnion
    LargeXMonimolimnion_pa <- as.data.frame(rowSums(tableVinput_richness_normalize_annots_tax_Funct %>% dplyr::select(all_of(LargeXMonimolimnion))))
    colnames(LargeXMonimolimnion_pa) <- "LargeXMonimolimnion"
    ## LargeXMixolimnion
    LargeXMixolimnion_pa <- as.data.frame(rowSums(tableVinput_richness_normalize_annots_tax_Funct %>% dplyr::select(all_of(LargeXMixolimnion))))
    colnames(LargeXMixolimnion_pa) <- "LargeXMixolimnion"
    ## Merge
    FractionXZone_pa <- merge(x = SmallXMonimolimnion_pa,y = SmallXMixolimnion_pa, by = "row.names")
    FractionXZone_pa <- merge(x = FractionXZone_pa,y = LargeXMonimolimnion_pa, by.x = "Row.names", by.y = "row.names")
    FractionXZone_pa <- merge(x = FractionXZone_pa,y = LargeXMixolimnion_pa, by.x = "Row.names", by.y = "row.names")
  #Periods
    ## 04
    `04_seq` <- as.data.frame(rowSums(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(all_of(Periods04))))
    colnames(`04_seq`) <- "Total04"
    ## 06
    `06_seq` <- as.data.frame(rowSums(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(all_of(Periods06))))
    colnames(`06_seq`) <- "Total06"
    ## 09
    `09_seq` <- as.data.frame(rowSums(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(all_of(Periods09))))
    colnames(`09_seq`) <- "Total09"
    ## 11
    `11_seq` <- as.data.frame(rowSums(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select(all_of(Periods11))))
    colnames(`11_seq`) <- "Total11"
    ## Merge
    Periods_seq <- merge(x = `04_seq`,y = `06_seq`, by = "row.names")
    Periods_seq <- merge(x = Periods_seq,y = `09_seq`, by.x = "Row.names", by.y = "row.names")
    Periods_seq <- merge(x = Periods_seq,y = `11_seq`, by.x = "Row.names", by.y = "row.names")
#
# Venn --------------------------------------------------------------------
    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
    # Cycle
    CycleVenn <- venn.diagram(x = list(Cycle_seq %>% filter(TotalDay > 0) %>% dplyr::select(Row.names) %>% unlist(), 
                                       Cycle_seq %>% filter(TotalNight > 0) %>% dplyr::select(Row.names) %>% unlist()), 
                              category.names = c("Day" , "Night"),
                              filename = NULL,
                              compression = "lzw",
                              lwd = 1,
                              main = "Cycles",
                              col=c("#a14647", '#6d6ea8'),
                              fill = c(alpha("#a14647",0.3), alpha('#6d6ea8',0.3)),
                              cex = 0.8,
                              main.fontfamily = "sans",
                              main.fontface = "bold",
                              fontfamily = "sans",
                              cat.cex = 1,
                              cat.default.pos = "outer",
                              cat.pos = c(-35,35),
                              cat.fontfamily = "sans",
                              cat.fontface = "bold",
                              cat.col = c("#a14647", '#6d6ea8'))
    ggsave(CycleVenn, file="Venn/CycleVenn.svg", device = "svg", width = 4, height = 4)
    # Fraction
    FractionVenn <- venn.diagram(x = list(Fraction_seq %>% filter(TotalSmall > 0) %>% dplyr::select(Row.names) %>% unlist(), 
                                          Fraction_seq %>% filter(TotalLarge > 0) %>% dplyr::select(Row.names) %>% unlist()), 
                                 category.names = c("Small" , "Large"),
                                 filename = NULL,
                                 compression = "lzw",
                                 lwd = 1,
                                 main = "Fractions",
                                 col=c('#6d6ea8',"#a14647"),
                                 fill = c(alpha("#6d6ea8",0.3), alpha('#a14647',0.3)),
                                 cex = 0.8,
                                 main.fontfamily = "sans",
                                 main.fontface = "bold",
                                 fontfamily = "sans",
                                 cat.cex = 1,
                                 cat.default.pos = "outer",
                                 cat.pos = c(-35,35),
                                 cat.fontfamily = "sans",
                                 cat.fontface = "bold",
                                 cat.col = c("#6d6ea8", '#a14647'))
    ggsave(FractionVenn, file="Venn/FractionVenn.svg", device = "svg", width = 4, height = 4)
    # Zone
    ZoneVenn <- venn.diagram(x = list(Zone_seq %>% filter(TotalMixolimnion > 0) %>% dplyr::select(Row.names) %>% unlist(), 
                                      Zone_seq %>% filter(TotalMonimolimnion > 0) %>% dplyr::select(Row.names) %>% unlist()), 
                             category.names = c("Mixolimnion" , "Monimolimnion"),
                             filename = NULL,
                             compression = "lzw",
                             lwd = 1,
                             main = "Zones",
                             col=c("#a14647", '#6d6ea8'),
                             fill = c(alpha("#a14647",0.3), alpha('#6d6ea8',0.3)),
                             cex = 0.8,
                             main.fontfamily = "sans",
                             main.fontface = "bold",
                             fontfamily = "sans",
                             cat.cex = 1,
                             cat.default.pos = "outer",
                             cat.pos = c(-30,30),
                             cat.fontfamily = "sans",
                             cat.fontface = "bold",
                             cat.col = c("#a14647", '#6d6ea8'))
    ggsave(ZoneVenn, file="Venn/ZoneVenn.svg", device = "svg", width = 4, height = 4)
    # Periods
    PeriodsVenn <- venn.diagram(x = list(Periods_seq %>% filter(Total04 > 0) %>% dplyr::select(Row.names) %>% unlist(),
                                         Periods_seq %>% filter(Total06 > 0) %>% dplyr::select(Row.names) %>% unlist(),
                                         Periods_seq %>% filter(Total09 > 0) %>% dplyr::select(Row.names) %>% unlist(),
                                         Periods_seq %>% filter(Total11 > 0) %>% dplyr::select(Row.names) %>% unlist()), 
                             category.names = c("April" , "June", "September", "November"),
                             filename = NULL,
                             compression = "lzw",
                             lwd = 1,
                             main = "Periods",
                             col=c("#6d6ea8", '#7b8e4b', '#e3a235', '#a14647'),
                             fill = c(alpha("#6d6ea8",0.3), alpha('#7b8e4b',0.3), alpha('#e3a235',0.3), alpha('#a14647',0.3)),
                             cex = 0.8,
                             main.fontfamily = "sans",
                             main.fontface = "bold",
                             fontfamily = "sans",
                             cat.cex = 1,
                             cat.default.pos = "outer",
                             #cat.pos = c(-35,35),
                             cat.fontfamily = "sans",
                             cat.fontface = "bold",
                             cat.col = c("#6d6ea8", '#7b8e4b', '#e3a235', '#a14647'))
    ggsave(PeriodsVenn, file="Venn/PeriodsVenn.svg", device = "svg", width = 4, height = 4)
    
    # Coplot
    svglite("Venn/VENCDT.svg",width = 8.00,height = 8.00)
    b_plot <- plot_grid(CycleVenn,NULL,ZoneVenn,NULL,NULL,NULL,FractionVenn,NULL,PeriodsVenn, labels = c('A', '','B','','','', 'C','', 'D'), ncol = 3, nrow = 3, rel_widths = c(5,1,5),rel_heights = c(5,1,5))
    print(b_plot)
    dev.off()
    
    #KO_seq_FractionXZoneXCycle <- separate_rows(KO_seq_FractionXZoneXCycle,"Kegg.Onthology",sep = ",") # Duplicate row by KO 
    #KO_seq_FractionXZoneXCycle <- KO_seq_FractionXZoneXCycle %>% dplyr::select(-"Row.names") %>% group_by(Kegg.Onthology,nbGroup) %>% summarise_all(sum)
    #
    #colSums(KO_seq_FractionXZoneXCycle %>% filter(is.na(nbGroup)==FALSE) %>% filter(!nbGroup %in% c("Unassigned","Lichen","Other multicellular Fungi", "Embryophyceae","Unaffiliated")) %>% dplyr::select(-Row.names,-nbGroup,-Kegg.Onthology))
    
    
# FractionXZone Taxonomy analysis --------------------------------------------------
  paletspe <- rep(c("#e53a35","#f44336","#ef5350","#e57373","#ef9a9a","#f4d4d4",
                    "#7b1fa2","#8e24aa","#8e24aaFF",
                    "#303f9f","#303f9fFF",
                    "#3065c0","#3f88e5","#049be5FF","#64b5f6","#b8d9fd","#64b5f680",
                    "#2e7d32","#689f38",#"#7cb342",
                    "#fbc02c","#fdd835","#ffeb3a",#"#ffee58","#fff176",
                    "#c6c6c6",#"#f49800","#f57c0099","#f4980080",
                    "#616161","#757575","#9e9e9e","#bdbdbd",#"#e0e0e0",
                    "#f57c00",
                    "#1876d2","#1876d295",#"#4cc6da","#80deea","#b2ebf2",
                    "#3897a7","#41acc1","#4cc6da","#80deea","#b2ebf2"),
                  4)
    ko2Met_Table_Final_i <- ko2Met_Table_Final %>% dplyr::select('Pathway/Gene',KO_id) %>% distinct()
    ko2Met_Table_Final_x <- ko2Met_Table_Final %>% dplyr::select(-lvl_C_val,-lvl_A_id,-lvl_B_id,-lvl_C_id,-lvl_B_val) %>% distinct()
    write.table(ko2Met_Table_Final_x,"Table/KO2Met.tsv",row.names = FALSE,col.names = TRUE,quote = FALSE, sep ="\t")
    
  # Taxonomy ----------------------------------------------------------------
  data_seq_FractionXZone <- dplyr::left_join(FractionXZone_seq,tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select("Row.names","Supergroup","Division"),by="Row.names")
  data_seq_FractionXZone <- data_seq_FractionXZone %>% dplyr::select(-"Row.names") %>% group_by(Supergroup,Division) %>% summarise_all(sum)
  # Other
  data_seq_FractionXZone$Somme <- rowSums(data_seq_FractionXZone[,3:6])
  data_seq_FractionXZone[,"Division"][data_seq_FractionXZone[,"Somme"]<0.00005*colSums(data_seq_FractionXZone[,"Somme"])] <- "Other<0.005%"
  data_seq_FractionXZone <- data_seq_FractionXZone %>% dplyr::select(-Somme)
  data_seq_FractionXZone[,"Supergroup"][data_seq_FractionXZone[,"Division"]=="Other<0.005%"] <- "Other<0.005%"
  data_seq_FractionXZone <- data_seq_FractionXZone %>% group_by(Supergroup,Division) %>% summarise_all(sum)
  #
  data_seq_FractionXZone <- reshape2::melt(data_seq_FractionXZone,id.vars = c("Supergroup","Division"))
  data_seq_FractionXZone <- data_seq_FractionXZone %>% mutate(Division = ifelse(is.na(Division) == TRUE, paste0('Unaffiliated_',Supergroup),Division))
  data_seq_FractionXZone[data_seq_FractionXZone=="Unaffiliated_NA"] <- NA
  data_seq_FractionXZone$Sum <- 0
  data_seq_FractionXZone$Prop <- 0
  for (i in row.names(data_seq_FractionXZone)) { var <- data_seq_FractionXZone[i,"variable"]
  data_seq_FractionXZone[i,"Sum"] <- sum(data_seq_FractionXZone %>% filter(variable == var) %>% dplyr::select(value))
  data_seq_FractionXZone[i,"Proportion"] <- data_seq_FractionXZone[i,"value"]*100/data_seq_FractionXZone[i,"Sum"]}
  data_seq_FractionXZone <- separate(data_seq_FractionXZone,"variable",c("Fraction","Zone"),sep="X")
  data_seq_FractionXZone[is.na(data_seq_FractionXZone)==TRUE] <- "Not Affiliated"
  data_seq_FractionXZone <- data_seq_FractionXZone %>% arrange(Supergroup)
  for (i in row.names(data_seq_FractionXZone)) { if (grepl(pattern="_X",data_seq_FractionXZone[i,"Division"])==TRUE) {data_seq_FractionXZone[i,"Division"] <- paste0("Unaffiliated_",data_seq_FractionXZone[i,"Supergroup"])}}
  OrderDiv <- unique(data_seq_FractionXZone$Division)
  #plot
  Total_ZoneXFraction_Sequence_fig <- ggplot(data_seq_FractionXZone, mapping = aes(y= Proportion, x = Zone, fill = factor(Division,level=OrderDiv), group = factor(Division,level=OrderDiv)), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity",color="black") +
    scale_fill_manual(values = paletspe, na.value = NA) + facet_grid(.~Fraction,scales="free") + theme_unique_art() +
    labs(x="Zones",y="Reads abundance (RPKM%)",fill="Division")
  print(Total_ZoneXFraction_Sequence_fig)
  ggsave("Hist/Tax-Div_abundance.svg", device = "svg", width = 10, height = 10)
  
  # Function ----------------------------------------------------------------

    # Hist plot --------------------------------------------------------------------
  ## Metabolisme
  #data_seq_FractionXZone <- dplyr::left_join(FractionXZone_seq,tableVinput_abundance_normalize_annots_tax_Funct %>% filter(Division == "Ochrophyta") %>% dplyr::select("Row.names","Kegg.Onthology"),by="Row.names")
  data_seq_FractionXZone <- dplyr::left_join(FractionXZone_seq,tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select("Row.names","Kegg.Onthology"),by="Row.names")
  data_seq_FractionXZone <- separate_rows(data_seq_FractionXZone,"Kegg.Onthology",sep = ",") # Duplicate row by KO 
  data_seq_FractionXZone <- data_seq_FractionXZone %>% dplyr::select(-"Row.names") %>% group_by(Kegg.Onthology) %>% summarise_all(sum)
  #
  data_seq_FractionXZone_i <- left_join(data_seq_FractionXZone,ko2Met_Table_Final_i,by=c("Kegg.Onthology"="KO_id"))
  colnames(data_seq_FractionXZone_i)[colnames(data_seq_FractionXZone_i) == 'Pathway/Gene'] <- 'Label'
  #
  data_seq_FractionXZone_i[,"Label"] <- sapply(strsplit(as.character(data_seq_FractionXZone_i$Label),"_-_"), `[`, 1)
  data_seq_FractionXZone_i[,"Label"] <- sapply(str_replace_all(as.character(data_seq_FractionXZone_i$Label),"_"," "), `[`, 1)
  data_seq_FractionXZone_i <- data_seq_FractionXZone_i %>% dplyr::select(-"Kegg.Onthology") %>% group_by(Label) %>% summarise_all(sum)
  data_seq_FractionXZone_i[is.na(data_seq_FractionXZone_i)==TRUE] <- "Unaffiliated"
  data_seq_FractionXZone_i <- reshape2::melt(data_seq_FractionXZone_i,id.vars = c("Label"))
  data_seq_FractionXZone_i$Sum <- 0
  for (i in row.names(data_seq_FractionXZone_i)) { var <- data_seq_FractionXZone_i[i,"variable"]
  data_seq_FractionXZone_i[i,"Sum"] <- sum(data_seq_FractionXZone_i %>% filter(variable == var) %>% dplyr::select(value))
  data_seq_FractionXZone_i[i,"Proportion"] <- data_seq_FractionXZone_i[i,"value"]*100/data_seq_FractionXZone_i[i,"Sum"]}
  data_seq_FractionXZone_i <- separate(data_seq_FractionXZone_i,"variable",c("Fraction","Zone"),sep="X")
  data_seq_FractionXZone_i<-data_seq_FractionXZone_i %>% filter(Label!="Unaffiliated")
  #plot
  color_ko_table_i <- color_ko_table %>% filter(Label %in% data_seq_FractionXZone_i$Label)
  Total_ZoneXFraction_Sequence_fig <- ggplot(data_seq_FractionXZone_i, mapping = aes(y= Proportion, x = Zone, fill = factor(Label, levels = color_ko_table_order), group = factor(Label, levels = color_ko_table_order)), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity",color="black") +
    scale_fill_manual(values=color_ko_table_i$palet) + facet_grid(.~Fraction,scales="free") + theme_unique_art() +
    labs(x="Zones",y="Unigenes %",fill="Metabolism")
  print(Total_ZoneXFraction_Sequence_fig)
  ggsave("Hist/KO-Div_abundance.svg", device = "svg", width = 10, height = 10)

    # Balance plot ------------------------------------------------------------
  ## Mixolimnion
  data_seq_FractionXZone_i_mixolimnion <- data_seq_FractionXZone_i %>% filter(Zone == "Mixolimnion") %>% dplyr::select(Label,value,Zone,Fraction,Sum) %>% group_by(Label,Zone,Fraction,Sum) %>% summarise_all(sum)
  ## Monimolimnion
  data_seq_FractionXZone_i_monimolimnion <- data_seq_FractionXZone_i %>% filter(Zone == "Monimolimnion") %>% dplyr::select(Label,value,Zone,Fraction,Sum) %>% group_by(Label,Zone,Fraction,Sum) %>% summarise_all(sum)
  ## Bind
  balance_totalZoneXFractionUnigen<-data_seq_FractionXZone_i_mixolimnion %>% dplyr::select(-Sum)
  for (i in row.names(balance_totalZoneXFractionUnigen)) { balance_totalZoneXFractionUnigen[i,"value"] <- foldchange(data_seq_FractionXZone_i_mixolimnion[i,"value"],data_seq_FractionXZone_i_monimolimnion[i,"value"])
  balance_totalZoneXFractionUnigen[i,"value"] <- foldchange2logratio(balance_totalZoneXFractionUnigen[i,"value"][[1]])
  if (data_seq_FractionXZone_i_mixolimnion[i,"value"] < data_seq_FractionXZone_i_monimolimnion[i,"value"]) {balance_totalZoneXFractionUnigen[i,"Zone"] <- "Monimolimnion"}
  if (data_seq_FractionXZone_i_mixolimnion[i,"value"] > data_seq_FractionXZone_i_monimolimnion[i,"value"]) {balance_totalZoneXFractionUnigen[i,"Zone"] <- "Mixolimnion"}
  balance_totalZoneXFractionUnigen[i,"label"] <- paste(round(max(data_seq_FractionXZone_i_mixolimnion[i,"value"],data_seq_FractionXZone_i_monimolimnion[i,"value"]),1),"VS",round(min(data_seq_FractionXZone_i_mixolimnion[i,"value"],data_seq_FractionXZone_i_monimolimnion[i,"value"]),1),sep=" ")
  balance_totalZoneXFractionUnigen[i,"Abundance"] <- max(data_seq_FractionXZone_i_mixolimnion[i,"value"],data_seq_FractionXZone_i_monimolimnion[i,"value"])}
  #
  #for (i in unique(balance_totalZoneXFractionUnigen$Label)) { if (sum(abs(c(balance_totalZoneXFractionUnigen %>% filter(Label==i))$value))<0.5) { balance_totalZoneXFractionUnigen[,"Label"][balance_totalZoneXFractionUnigen[,"Label"]==i]<-"∑(log2)<0.5"}}
  #balance_totalZoneXFractionUnigen <- balance_totalZoneXFractionUnigen %>% filter(Label != "∑(log2)<0.5")
  #
  ## plot
  balance_palet <- paletteer_d("ggthemes::Nuriel_Stone",n=2, direction=-1)
  balxi <- ggplot(balance_totalZoneXFractionUnigen, aes(x = factor(Label,level=rev(color_ko_table_order)), y = value)) +
    geom_bar(aes(fill = Zone,alpha=log10(Abundance)), stat = 'identity',color="black") +  facet_grid(.~Fraction) +
    coord_flip() +
    geom_text(aes(label = label ,y = 0,vjust = ifelse(value >= 0, 0.5, 0.5), hjust = ifelse(value>= 0, 1.2, -0.2)),color = "black",size = 3,show.legend = FALSE) +
    labs(x="",y="TPM log2-ratio",fill="Zones") +
    theme(axis.text.y = element_blank()) + theme_unique_art() + scale_fill_manual(values=balance_palet) + guides(alpha = "none")
  print(balxi)
  ggsave("Balance/KO-Balance_abundance.svg", device = "svg", width = 16.00, height = 12)
  #
# FractionXZone Funct° grp analysis --------------------------------------------------
  palet_tree_grp <- c("#607848FF","#ad0635FF","#d84040FF","#4080c0FF","#60d0a0FF","#70a830FF","#98d048FF","#d7c11cFF","#40a078FF","#333b43FF","#2c7890FF","#b4d2d2")
  
  # Polar -------------------------------------------------------------------
    # General -----------------------------------------------------------------
      # Richness P/A ------------------------------------------------------------
  data_pa_FractionXZone <- dplyr::left_join(FractionXZone_pa,tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select("Row.names","nbGroup","Supergroup","Division"),by="Row.names")
  data_pa_FractionXZone[,"nbGroup"][is.na(data_pa_FractionXZone[,"nbGroup"])==TRUE] <- "Unassigned"
  Polar_pa <- data_pa_FractionXZone
  Polar_pa$TotalAnnee <- 0
  Polar_pa[,"TotalAnnee"][Polar_pa[,"SmallXMonimolimnion"] + Polar_pa[,"SmallXMixolimnion"] + Polar_pa[,"LargeXMonimolimnion"] + Polar_pa[,"LargeXMixolimnion"] > 0] <- 1
  Polar_pa <- Polar_pa %>% dplyr::select(TotalAnnee,nbGroup)
  ## Total Figure
  Polar_pa <- Polar_pa %>% dplyr::group_by(nbGroup) %>% dplyr::summarise_all(sum) %>% mutate(Total = TotalAnnee*100/sum(TotalAnnee))
  ## Label
  Polar_pa <- Polar_pa %>%
    dplyr::arrange(desc(nbGroup)) %>%
    dplyr::mutate(lab.ypos = cumsum(Total) - 0.5*Total)
  Polar_pa$label <- paste(round(Polar_pa$Total,1), "%", sep = "")
  for (i in rownames(Polar_pa)) {
    if (Polar_pa[i,"label"] == "0%") { Polar_pa[i,"label"] <- NA}}
  for (i in rownames(Polar_pa)) {
    if (is.na(Polar_pa[i,"label"]) == FALSE) { Polar_pa[i,"label"] <- paste(Polar_pa[i,"nbGroup"]," : ",Polar_pa[i,"label"], sep = "")}}
  ## Figure
  bxim <- ggplot(Polar_pa, mapping = aes(y= Total, x = 2, fill = nbGroup), Rowv = NA, col = colMain, scale = "column") +
    geom_bar(stat="identity", color = "white", width = 1,alpha=0.8) + coord_polar("y") + 
    geom_label_repel(aes(y = lab.ypos,label = label),color = "white",size = 3,segment.color = "black",show.legend = FALSE, nudge_x = 0.5) +
    scale_y_continuous(limits=c(0,sum(Polar_pa %>% dplyr::select(Total)+0.01)))+
    xlim(1,2.5) +
    theme_unique_darkbis() + 
    labs(y = "P/A unigene",x="") + scale_fill_manual(values=palet_tree_grp)
  svglite("Polar/Polar-pa-unigene.svg",width = 4.50,height = 4.50)
  print(bxim)
  dev.off()

      # Abundance ---------------------------------------------------------------
  data_seq_FractionXZone <- dplyr::left_join(FractionXZone_seq,tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select("Row.names","nbGroup","Supergroup","Division"),by="Row.names")
  data_seq_FractionXZone[,"nbGroup"][is.na(data_seq_FractionXZone[,"nbGroup"])==TRUE] <- "Unassigned"
  Polar_seq <- data_seq_FractionXZone
  Polar_seq$TotalAnnee <- 0
  Polar_seq[,"TotalAnnee"] <- Polar_seq[,"SmallXMonimolimnion"] + Polar_seq[,"SmallXMixolimnion"] + Polar_seq[,"LargeXMonimolimnion"] + Polar_seq[,"LargeXMixolimnion"]
  Polar_seq <- Polar_seq %>% dplyr::select(TotalAnnee,nbGroup)
  ## Total Figure
  Polar_seq <- Polar_seq %>% dplyr::group_by(nbGroup) %>% dplyr::summarise_all(sum) %>% mutate(Total = TotalAnnee*100/sum(TotalAnnee))
  ## Label
  Polar_seq <- Polar_seq %>%
    dplyr::arrange(desc(nbGroup)) %>%
    dplyr::mutate(lab.ypos = cumsum(Total) - 0.5*Total)
  Polar_seq$label <- paste(round(Polar_seq$Total,1), "%", sep = "")
  for (i in rownames(Polar_seq)) {
    if (Polar_seq[i,"label"] == "0%") { Polar_seq[i,"label"] <- NA}}
  for (i in rownames(Polar_seq)) {
    if (is.na(Polar_seq[i,"label"]) == FALSE) { Polar_seq[i,"label"] <- paste(Polar_seq[i,"nbGroup"]," : ",Polar_seq[i,"label"], sep = "")}}
  ## Figure
  bxim <- ggplot(Polar_seq, mapping = aes(y= Total, x = 2, fill = nbGroup), Rowv = NA, col = colMain, scale = "column") +
    geom_bar(stat="identity", color = "white", width = 1,alpha=0.8) + coord_polar("y") + 
    geom_label_repel(aes(y = lab.ypos,label = label),color = "white",size = 3,segment.color = "black",show.legend = FALSE, nudge_x = 0.5) +
    scale_y_continuous(limits=c(0,sum(Polar_seq %>% dplyr::select(Total)+0.01)))+
    xlim(1,2.5) +
    theme_unique_darkbis() + 
    labs(y = "Abundance unigene",x="") + scale_fill_manual(values=palet_tree_grp)
  svglite("Polar/Polar-seq-unigene.svg",width = 4.50,height = 4.50)
  print(bxim)
  dev.off()
  
  
    # By group ----------------------------------------------------------------

      # Richness P/A ------------------------------------------------------------
  # Processing
  data_pa_FractionXZone <- dplyr::left_join(FractionXZone_pa,tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select("Row.names","nbGroup","Supergroup","Division"),by="Row.names")
  data_pa_FractionXZone[,"nbGroup"][is.na(data_pa_FractionXZone[,"nbGroup"])==TRUE] <- "Unassigned"
  data_pa_FractionXZone <- data_pa_FractionXZone %>% mutate(Division = ifelse(is.na(Division) == TRUE, paste0('Unaffiliated_',Supergroup),Division))
  data_pa_FractionXZone[data_pa_FractionXZone=="Unaffiliated_NA"] <- NA
  k=0
  for (l in unique(data_pa_FractionXZone$nbGroup)) { print(l)
    if ( !(l %in% c("Unassigned","Lichen","Other multicellular Fungi", "Embryophyceae"))) {
    k<-k+1
    ## Total column
    Polar_pa <- data_pa_FractionXZone
    Polar_pa <- Polar_pa %>% dplyr::filter(nbGroup==l)
    Polar_pa$TotalAnnee <- 0
    Polar_pa[,"TotalAnnee"][Polar_pa[,"SmallXMonimolimnion"] + Polar_pa[,"SmallXMixolimnion"] + Polar_pa[,"LargeXMonimolimnion"] + Polar_pa[,"LargeXMixolimnion"] > 0] <- 1
    Polar_pa <- Polar_pa %>% dplyr::select(TotalAnnee,Division)
    ## Total Figure
    Polar_pa <- Polar_pa %>% dplyr::group_by(Division) %>% dplyr::summarise_all(sum) %>% mutate(Total = TotalAnnee*100/sum(TotalAnnee))
    ## Label
    Polar_pa <- Polar_pa %>%
      dplyr::arrange(desc(Division)) %>%
      dplyr::mutate(lab.ypos = cumsum(Total) - 0.5*Total)
    Polar_pa$label <- paste(round(Polar_pa$Total,1), "%", sep = "")
    for (i in rownames(Polar_pa)) {
      if (Polar_pa[i,"label"] == "0%") { Polar_pa[i,"label"] <- NA}}
    for (i in rownames(Polar_pa)) {
      if (is.na(Polar_pa[i,"label"]) == FALSE) { Polar_pa[i,"label"] <- paste(Polar_pa[i,"Division"]," : ",Polar_pa[i,"label"], sep = "")}}
    Polar_pa$Group<-l
    if (k==1) {Polar_pa_GRP<-Polar_pa}
    if (k>1) {Polar_pa_GRP <- rbind(Polar_pa_GRP,Polar_pa)}
    }}
  #plot
  paletxi <- c("#c04848CC","#f27435CC","#707470CC","#948c75CC","#424254CC","#e04644CC","#706acfCC","#6a4a3cCC","#2a044aCC","#e84a5fCC","#f02475CC","#e8bf56CC","#1b676bCC","#3b2d38CC","#cc2a41CC","#f02475CC","#cfbe27CC","#b8af03CC","#d27f48CC","#2cb0e0CC","#3fb8afCC","#3B3B3B99","#519548CC","#c2d1fcCC","#555152CC","#8fbe00CC","#64908aCC","#0e2430CC","#bcbdacCC")
  bx <- ggplot(Polar_pa_GRP, mapping = aes(y= Total, x = 2, fill = Division), Rowv = NA, col = colMain, scale = "column") +
    geom_bar(stat="identity", color = "white", width = 1) + coord_polar("y") + 
    geom_label_repel(aes(y = lab.ypos,label = label),color = "white",size = 3,segment.color = "black",show.legend = FALSE, nudge_x = 0.6,nudge_y = 0.5) +
    scale_y_continuous(limits=c(0,100.1))+
    xlim(1,2.5) +
    theme_unique_darkbis() + 
    facet_wrap( ~ Group , nrow = 3) +
    labs(y = "P/A unigene",x="")  + scale_fill_manual(values = paletxi) #+ scale_fill_paletteer_d("ggthemes::manyeys")
  svglite("Polar/Polar-pa-unigene-bygrp.svg",width = 10,height = 10)
  print(bx)
  dev.off()
      # Abundance --------------------------------------------------------------
  # Processing
  data_seq_FractionXZone <- dplyr::left_join(FractionXZone_seq,tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select("Row.names","nbGroup","Supergroup","Division"),by="Row.names")
  data_seq_FractionXZone[,"nbGroup"][is.na(data_seq_FractionXZone[,"nbGroup"])==TRUE] <- "Unassigned"
  data_seq_FractionXZone <- data_seq_FractionXZone %>% mutate(Division = ifelse(is.na(Division) == TRUE, paste0('Unaffiliated_',Supergroup),Division))
  data_seq_FractionXZone[data_seq_FractionXZone=="Unaffiliated_NA"] <- NA
  k=0
  for (l in unique(data_seq_FractionXZone$nbGroup)) { print(l)
    if ( !(l %in% c("Unassigned","Lichen","Other multicellular Fungi", "Embryophyceae"))) {
      ## Total column
      k<-k+1
      Polar_seq <- data_seq_FractionXZone
      Polar_seq <- Polar_seq %>% dplyr::filter(nbGroup==l)
      Polar_seq$TotalAnnee <- 0
      Polar_seq[,"TotalAnnee"] <- Polar_seq[,"SmallXMonimolimnion"] + Polar_seq[,"SmallXMixolimnion"] + Polar_seq[,"LargeXMonimolimnion"] + Polar_seq[,"LargeXMixolimnion"]
      Polar_seq <- Polar_seq %>% dplyr::select(TotalAnnee,Division)
      ## Total Figure
      Polar_seq <- Polar_seq %>% dplyr::group_by(Division) %>% dplyr::summarise_all(sum) %>% mutate(Total = TotalAnnee*100/sum(TotalAnnee))
      ## Label
      Polar_seq <- Polar_seq %>%
        dplyr::arrange(desc(Division)) %>%
        dplyr::mutate(lab.ypos = cumsum(Total) - 0.5*Total)
      Polar_seq$label <- paste(round(Polar_seq$Total,1), "%", sep = "")
      for (i in rownames(Polar_seq)) {
        if (Polar_seq[i,"label"] == "0%") { Polar_seq[i,"label"] <- NA}}
      for (i in rownames(Polar_seq)) {
        if (is.na(Polar_seq[i,"label"]) == FALSE) { Polar_seq[i,"label"] <- paste(Polar_seq[i,"Division"]," : ",Polar_seq[i,"label"], sep = "")}}
      Polar_seq$Group<-l
      if (k==1) {Polar_seq_GRP<-Polar_seq}
      if (k>1) {Polar_seq_GRP <- rbind(Polar_seq_GRP,Polar_seq)}
    }}
  #plot
  paletxi <- c("#c04848CC","#f27435CC","#707470CC","#948c75CC","#424254CC","#e04644CC","#706acfCC","#6a4a3cCC","#2a044aCC","#e84a5fCC","#f02475CC","#e8bf56CC","#1b676bCC","#3b2d38CC","#cc2a41CC","#f02475CC","#cfbe27CC","#b8af03CC","#d27f48CC","#2cb0e0CC","#3fb8afCC","#3B3B3B99","#519548CC","#c2d1fcCC","#555152CC","#8fbe00CC","#64908aCC","#0e2430CC","#bcbdacCC")
  bx <- ggplot(Polar_seq_GRP, mapping = aes(y= Total, x = 2, fill = Division), Rowv = NA, col = colMain, scale = "column") +
    geom_bar(stat="identity", color = "white", width = 1) + coord_polar("y") + 
    geom_label_repel(aes(y = lab.ypos,label = label),color = "white",size = 3,segment.color = "black",show.legend = FALSE, nudge_x = 0.6,nudge_y = 0.5) +
    scale_y_continuous(limits=c(0,100.1)) +
    xlim(1,2.5) +
    theme_unique_darkbis() + 
    facet_wrap( ~ Group , nrow = 3) +
    labs(y = "P/A unigene",x="")  + scale_fill_manual(values = paletxi) #+ scale_fill_paletteer_d("ggthemes::manyeys")
  svglite("Polar/Polar-seq-unigene-bygrp.svg",width = 10,height = 10)
  print(bx)
  dev.off()
  # Hist general ------------------------------------------------------------
  paletspe_i <- rep(c("#e53a35","#f44336","#ef5350",
                      "#303f9f","#303f9fFF",
                      "#3065c0","#3f88e5",
                      "#2e7d32",
                      "#fbc02c","#fdd835",
                      "#616161","#757575",
                      "#1876d2",
                      "#3897a7","#41acc1"),
                    4)
    # KO Hist ------------------------------------------------------------
  KO_seq_FractionXZone <- dplyr::left_join(FractionXZone_seq,tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select("Row.names","nbGroup","Kegg.Onthology"),by="Row.names")
  KO_seq_FractionXZone <- separate_rows(KO_seq_FractionXZone,"Kegg.Onthology",sep = ",") # Duplicate row by KO 
  KO_seq_FractionXZone <- KO_seq_FractionXZone %>% dplyr::select(-"Row.names") %>% group_by(Kegg.Onthology,nbGroup) %>% summarise_all(sum)
  #
  KO_seq_FractionXZone_i <- left_join(KO_seq_FractionXZone,ko2Met_Table_Final_i,by=c("Kegg.Onthology"="KO_id"))
  colnames(KO_seq_FractionXZone_i)[colnames(KO_seq_FractionXZone_i) == 'Pathway/Gene'] <- 'Label'
  #
  KO_seq_FractionXZone_i[,"Label"] <- sapply(strsplit(as.character(KO_seq_FractionXZone_i$Label),"_-_"), `[`, 1)
  KO_seq_FractionXZone_i[,"Label"] <- sapply(str_replace_all(as.character(KO_seq_FractionXZone_i$Label),"_"," "), `[`, 1)
  KO_seq_FractionXZone_i <- as.data.frame(KO_seq_FractionXZone_i)
  KO_seq_FractionXZone_i <- KO_seq_FractionXZone_i %>% dplyr::select(-"Kegg.Onthology") %>% group_by(nbGroup,Label) %>% summarise_all(sum)
  KO_seq_FractionXZone_i[is.na(KO_seq_FractionXZone_i)==TRUE] <- "Unaffiliated"
  #
  KO_seq_FractionXZone_i <- reshape2::melt(KO_seq_FractionXZone_i,id.vars = c("Label","nbGroup"))
  KO_seq_FractionXZone_i <- separate(KO_seq_FractionXZone_i,"variable",c("Fraction","Zone"),sep="X")
  KO_seq_FractionXZone_i <- KO_seq_FractionXZone_i %>% filter(Label != "Unaffiliated") %>% filter(!nbGroup %in% c("Unassigned","Lichen","Other multicellular Fungi", "Embryophyceae","Unaffiliated"))
  for (i in row.names(KO_seq_FractionXZone_i)) {
    Fractionx <- KO_seq_FractionXZone_i[i,"Fraction"]
    Zonex <- KO_seq_FractionXZone_i[i,"Zone"]
    nbGroupx <- KO_seq_FractionXZone_i[i,"nbGroup"]
    KO_seq_FractionXZone_i[i,"Sum"] <- sum(KO_seq_FractionXZone_i %>% filter(Fraction == Fractionx) %>% filter(Zone == Zonex) %>% filter(nbGroup == nbGroupx) %>% dplyr::select(value))
    KO_seq_FractionXZone_i[i,"Relative"] <- KO_seq_FractionXZone_i[i,"Sum"] *100 / sum(KO_seq_FractionXZone_i %>% filter(Fraction == Fractionx) %>% filter(Zone == Zonex) %>% dplyr::select(value))
    KO_seq_FractionXZone_i[i,"Proportion"] <- KO_seq_FractionXZone_i[i,"value"]*100/KO_seq_FractionXZone_i[i,"Sum"]}
  KO_seq_FractionXZone_i[,"nbGroup"] <- sapply(strsplit(as.character(KO_seq_FractionXZone_i$nbGroup),": "), `[`, 2)
  Metabo_seq_FractionXZone_k <- KO_seq_FractionXZone_i
  #plot
  color_ko_table_i <- color_ko_table %>% filter(Label %in% KO_seq_FractionXZone_i$Label)
  Total_ZoneXFraction_Sequence_fig <- ggplot(KO_seq_FractionXZone_i, mapping = aes(y= Proportion, x = Zone, fill = factor(Label, levels = color_ko_table_order), group = factor(Label, levels = color_ko_table_order)), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity",alpha=0.95,color="black") +
    scale_fill_manual(values = color_ko_table_i$palet) + facet_grid(Fraction~nbGroup,scales="free") + theme_unique_art() +
    labs(x="Zones",y="Relative abundance (RPKM %)",fill="Metabolisms/Cellular processes") + guides(fill=guide_legend(ncol=1)) +
    geom_label(aes(y = 106,label = paste0(round(Sum,1)," TPM","\n","(",round(Relative,1)," %)")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
  print(Total_ZoneXFraction_Sequence_fig)
  ggsave("Hist/Grp_ZoneXFraction_Function_abundance.svg", device = "svg", width = 18, height = 10)
  
  #plot
  color_ko_table_i <- color_ko_table %>% filter(Label %in% KO_seq_FractionXZone_i$Label)
  Total_ZoneXFraction_Sequence_fig <- ggplot(KO_seq_FractionXZone_i, mapping = aes(y= Proportion, x = nbGroup, fill = factor(Label, levels = color_ko_table_order), group = factor(Label, levels = color_ko_table_order)), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity",alpha=0.95,color="black") +
    scale_fill_manual(values = color_ko_table_i$palet) + theme_unique_art() + facet_grid(Fraction~Zone,scales="free") +
    labs(x="",y="Relative abundance (RPKM %)",fill="Metabolisms/Cellular processes") + theme(legend.position = "bottom") + guides(fill = guide_legend(title.position = "top",title.hjust = 0.5)) + theme(axis.title.x=element_blank()) +
    geom_label(aes(y = 106,label = paste0(round(Sum,1)," RPKM","\n","(",round(Relative,1)," %)")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
  print(Total_ZoneXFraction_Sequence_fig)
  ggsave("Hist/Grp_ZoneXFraction_Function_abundance_Four.svg", device = "svg", width = 14, height = 12)
  #      
# Without FrcationXZone aspect
  KO_seq_FractionXZone_i <- left_join(KO_seq_FractionXZone,ko2Met_Table_Final_i,by=c("Kegg.Onthology"="KO_id"))
  colnames(KO_seq_FractionXZone_i)[colnames(KO_seq_FractionXZone_i) == 'Pathway/Gene'] <- 'Label'
  KO_seq_FractionXZone_i$Total <- KO_seq_FractionXZone_i$SmallXMonimolimnion + KO_seq_FractionXZone_i$SmallXMixolimnion + KO_seq_FractionXZone_i$LargeXMonimolimnion + KO_seq_FractionXZone_i$LargeXMixolimnion
  KO_seq_FractionXZone_i <- KO_seq_FractionXZone_i %>% dplyr::select(-SmallXMixolimnion,-SmallXMonimolimnion,-LargeXMonimolimnion,-LargeXMixolimnion)
  #
  KO_seq_FractionXZone_i[,"Label"] <- sapply(strsplit(as.character(KO_seq_FractionXZone_i$Label),"_-_"), `[`, 1)
  KO_seq_FractionXZone_i[,"Label"] <- sapply(str_replace_all(as.character(KO_seq_FractionXZone_i$Label),"_"," "), `[`, 1)
  KO_seq_FractionXZone_i <- as.data.frame(KO_seq_FractionXZone_i)
  KO_seq_FractionXZone_i <- KO_seq_FractionXZone_i %>% dplyr::select(-"Kegg.Onthology") %>% group_by(nbGroup,Label) %>% summarise_all(sum)
  KO_seq_FractionXZone_i[is.na(KO_seq_FractionXZone_i)==TRUE] <- "Unaffiliated"
  #
  KO_seq_FractionXZone_i <- reshape2::melt(KO_seq_FractionXZone_i,id.vars = c("Label","nbGroup"))
  KO_seq_FractionXZone_i <- KO_seq_FractionXZone_i %>% filter(Label != "Unaffiliated") %>% filter(!nbGroup %in% c("Unassigned","Lichen","Other multicellular Fungi", "Embryophyceae","Unaffiliated"))
  #
  for (i in row.names(KO_seq_FractionXZone_i)) {
    nbGroupx <- KO_seq_FractionXZone_i[i,"nbGroup"]
    KO_seq_FractionXZone_i[i,"Sum"] <- sum(KO_seq_FractionXZone_i  %>% filter(nbGroup == nbGroupx) %>% dplyr::select(value))
    KO_seq_FractionXZone_i[i,"Relative"] <- KO_seq_FractionXZone_i[i,"Sum"] *100 / sum(KO_seq_FractionXZone_i %>% dplyr::select(value))
    KO_seq_FractionXZone_i[i,"Proportion"] <- KO_seq_FractionXZone_i[i,"value"]*100/KO_seq_FractionXZone_i[i,"Sum"]}
  KO_seq_FractionXZone_i[,"nbGroup"] <- sapply(strsplit(as.character(KO_seq_FractionXZone_i$nbGroup),": "), `[`, 2)
  #plot
  color_ko_table_i <- color_ko_table %>% filter(Label %in% KO_seq_FractionXZone_i$Label)
  Total_ZoneXFraction_Sequence_fig <- ggplot(KO_seq_FractionXZone_i, mapping = aes(y= Proportion, x = nbGroup, fill = factor(Label, levels = color_ko_table_order), group = factor(Label, levels = color_ko_table_order)), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity",alpha=0.95,color="black") +
    scale_fill_manual(values = color_ko_table_i$palet) + theme_unique_art() +
    labs(x="",y="Relative abundance (RPKM %)",fill="Metabolisms/Cellular processes") + theme(legend.position = "bottom") + guides(fill = guide_legend(title.position = "top",title.hjust = 0.5)) + theme(axis.title.x=element_blank()) +
    geom_label(aes(y = 106,label = paste0(round(Sum,1)," RPKM","\n","(",round(Relative,1)," %)")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
  print(Total_ZoneXFraction_Sequence_fig)
  ggsave("Hist/Grp_ZoneXFraction_Function_abundance_Total.svg", device = "svg", width = 14, height = 12)
    # KO Quaternary ------------------------------------------------------------
      KO_seq_FractionXZone <- dplyr::left_join(FractionXZone_seq,tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select("Row.names","nbGroup","Kegg.Onthology"),by="Row.names")
      KO_seq_FractionXZone <- separate_rows(KO_seq_FractionXZone,"Kegg.Onthology",sep = ",") # Duplicate row by KO 
      KO_seq_FractionXZone <- KO_seq_FractionXZone %>% dplyr::select(-"Row.names") %>% group_by(Kegg.Onthology,nbGroup) %>% summarise_all(sum)
      #
      #test
      #KO_seq_FractionXZone_i <- KO_seq_FractionXZone %>% filter(is.na(nbGroup)==FALSE) %>% filter(is.na(Kegg.Onthology)==FALSE) %>% filter(!nbGroup %in% c("Lichen","Other multicellular Fungi", "Embryophyceae","Unaffiliated"))
      #KO_seq_FractionXZone_i <- KO_seq_FractionXZone_i %>% filter(Kegg.Onthology %in% ko2Met_Table_Final_i$KO_id)
      #sum(KO_seq_FractionXZone_i[,"Somme"])
      #nrow(KO_seq_FractionXZone_i)
      # Quaternary plot Label ---------------------------------------------------
        KO_seq_FractionXZone_i <- left_join(KO_seq_FractionXZone,ko2Met_Table_Final_i,by=c("Kegg.Onthology"="KO_id"))
        colnames(KO_seq_FractionXZone_i)[colnames(KO_seq_FractionXZone_i) == 'Pathway/Gene'] <- 'Label'
        KO_seq_FractionXZone_i[,"Label"] <- sapply(strsplit(as.character(KO_seq_FractionXZone_i$Label),"_-_"), `[`, 1)
        KO_seq_FractionXZone_i[,"Label"] <- sapply(str_replace_all(as.character(KO_seq_FractionXZone_i$Label),"_"," "), `[`, 1)
        KO_seq_FractionXZone_i <- KO_seq_FractionXZone_i %>% filter(is.na(nbGroup)==FALSE) %>% filter(is.na(Kegg.Onthology)==FALSE) %>% filter(is.na(Label)==FALSE) %>% filter(!nbGroup %in% c("Lichen","Other multicellular Fungi", "Embryophyceae","Unaffiliated"))
        KO_seq_FractionXZone_i <- as.data.frame(KO_seq_FractionXZone_i) %>% dplyr::select(-"Kegg.Onthology") %>% group_by(nbGroup,Label) %>% summarise_all(sum)
        KO_seq_FractionXZone_i[,"Somme"] <- (KO_seq_FractionXZone_i[,"LargeXMonimolimnion"]+KO_seq_FractionXZone_i[,"LargeXMixolimnion"]+KO_seq_FractionXZone_i[,"SmallXMonimolimnion"]+KO_seq_FractionXZone_i[,"SmallXMixolimnion"])
        KO_seq_FractionXZone_i <- KO_seq_FractionXZone_i %>% filter(Somme > 10)
        ## Size
        KO_seq_FractionXZone_i[,"X"] <- (KO_seq_FractionXZone_i[,"SmallXMonimolimnion"]+KO_seq_FractionXZone_i[,"SmallXMixolimnion"])*100/(KO_seq_FractionXZone_i[,"LargeXMonimolimnion"]+KO_seq_FractionXZone_i[,"LargeXMixolimnion"]+KO_seq_FractionXZone_i[,"SmallXMonimolimnion"]+KO_seq_FractionXZone_i[,"SmallXMixolimnion"])
        ## Zone
        KO_seq_FractionXZone_i[,"Y"] <- (KO_seq_FractionXZone_i[,"SmallXMonimolimnion"]+KO_seq_FractionXZone_i[,"LargeXMonimolimnion"])*100/(KO_seq_FractionXZone_i[,"LargeXMonimolimnion"]+KO_seq_FractionXZone_i[,"LargeXMixolimnion"]+KO_seq_FractionXZone_i[,"SmallXMonimolimnion"]+KO_seq_FractionXZone_i[,"SmallXMixolimnion"])
        #Label nbgroup
        KO_seq_FractionXZone_i[,"nbGroup"] <- sapply(strsplit(as.character(KO_seq_FractionXZone_i$nbGroup),": "), `[`, 2)
        #
        # Plot #1
        color_ko_table_i <- color_ko_table %>% filter(Label %in% KO_seq_FractionXZone_i$Label)
        color_ko_table_i$fill <- c(rep("black",14),rep("#00000000",13))
        Quaternary_plot <- ggplot(KO_seq_FractionXZone_i, mapping = aes(y= Y, x = X)) + theme_unique_art() + 
          annotate("rect", xmin = 5, xmax = 95, ymin = 5, ymax = 95, color = "darkgrey", fill = "#00000000",linetype = "dashed") +
          annotate("rect", xmin = 12.5, xmax = 87.5, ymin = 12.5, ymax = 87.5, color = "grey", fill = "#00000000",linetype = "dashed") +
          annotate("rect", xmin = 20, xmax = 80, ymin = 20, ymax = 80, color = "lightgrey", fill = "#00000000",linetype = "dashed") +
          annotate("rect", xmin = 30, xmax = 70, ymin = 30, ymax = 70, color = "lightblue", fill = "#00000000",linetype = "dashed") +
          annotate("rect", xmin = 0, xmax = 100, ymin = 0, ymax = 100,
                   color = "black", fill = "#00000000") + 
          theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "right") + 
          geom_point(aes(size=Somme, color=factor(Label,levels = c(color_ko_table_i$Label)), fill=factor(Label,levels = c(color_ko_table_i$Label))), shape=21, stroke = 1) +
          scale_fill_manual(values = color_ko_table_i$palet) +
          scale_color_manual(values = color_ko_table_i$fill) +
          scale_size(range = c(3,10),breaks = c(500, 1000, 5000, 10000, 15000)) +
          labs(color="Metabolisms/Cellular processes",fill="Metabolisms/Cellular processes", size = "Total abundance (RPKM)") + 
          annotate("text", x = 5, y = 50, label = "Large", angle=90, fontface = "bold") +
          annotate("text", x = 95, y = 50, label = "Small", angle = -90, fontface = "bold") +
          annotate("text", x = 50, y = 5, label = "Mixolimnion", fontface = "bold") +
          annotate("text", x = 50, y = 95, label = "Monimolimnion", fontface = "bold") +
          geom_text_repel(data = subset(KO_seq_FractionXZone_i, Somme > 2000), size = 4,aes(label = "*"), 
                          box.padding = unit(0.45, "lines"),
                          max.overlaps = 1000) +
          guides(y.sec = "axis", size = guide_legend(title.position = "top",title.hjust = 0.5,ncol=1),
                 color = guide_legend(title.position = "top", override.aes = list(size=5),ncol=1),
                 fill = guide_legend(title.position = "top",ncol=1)) + 
          facet_wrap(nbGroup~.)
        print(Quaternary_plot)
        ggsave("Quaterneryplot/Quaterneryplot_abundance_Total_V.svg", device = "svg", width = 12, height = 12)
#
      # Quaternary plot KO ---------------------------------------------------
        KO_seq_FractionXZone_i <- KO_seq_FractionXZone %>% filter(is.na(nbGroup)==FALSE) %>% filter(is.na(Kegg.Onthology)==FALSE) %>% filter(!nbGroup %in% c("Lichen","Other multicellular Fungi", "Embryophyceae","Unaffiliated"))
        KO_seq_FractionXZone_i[,"Somme"] <- (KO_seq_FractionXZone_i[,"LargeXMonimolimnion"]+KO_seq_FractionXZone_i[,"LargeXMixolimnion"]+KO_seq_FractionXZone_i[,"SmallXMonimolimnion"]+KO_seq_FractionXZone_i[,"SmallXMixolimnion"])
        KO_seq_FractionXZone_i <- KO_seq_FractionXZone_i %>% filter(Kegg.Onthology %in% ko2Met_Table_Final_i$KO_id)
        KO_seq_FractionXZone_i <- KO_seq_FractionXZone_i %>% filter(Somme >= 100)
        ## Size
        KO_seq_FractionXZone_i[,"X"] <- (KO_seq_FractionXZone_i[,"SmallXMonimolimnion"]+KO_seq_FractionXZone_i[,"SmallXMixolimnion"])*100/(KO_seq_FractionXZone_i[,"LargeXMonimolimnion"]+KO_seq_FractionXZone_i[,"LargeXMixolimnion"]+KO_seq_FractionXZone_i[,"SmallXMonimolimnion"]+KO_seq_FractionXZone_i[,"SmallXMixolimnion"])
        ## Zone
        KO_seq_FractionXZone_i[,"Y"] <- (KO_seq_FractionXZone_i[,"SmallXMonimolimnion"]+KO_seq_FractionXZone_i[,"LargeXMonimolimnion"])*100/(KO_seq_FractionXZone_i[,"LargeXMonimolimnion"]+KO_seq_FractionXZone_i[,"LargeXMixolimnion"]+KO_seq_FractionXZone_i[,"SmallXMonimolimnion"]+KO_seq_FractionXZone_i[,"SmallXMixolimnion"])
        #Label nbgroup
        KO_seq_FractionXZone_i[,"nbGroup"] <- sapply(strsplit(as.character(KO_seq_FractionXZone_i$nbGroup),": "), `[`, 2)
        #
        KO_seq_FractionXZone_i[,"TXT"] <- NA
        for (i in row.names(KO_seq_FractionXZone_i)) { 
          if (KO_seq_FractionXZone_i[i,"Y"]>95 && KO_seq_FractionXZone_i[i,"Somme"]>100 ) { KO_seq_FractionXZone_i[i,"TXT"] <- KO_seq_FractionXZone_i[i,"Kegg.Onthology"]}
          if (KO_seq_FractionXZone_i[i,"X"]>95 && KO_seq_FractionXZone_i[i,"Somme"]>100 ) { KO_seq_FractionXZone_i[i,"TXT"] <- KO_seq_FractionXZone_i[i,"Kegg.Onthology"]}
          if (KO_seq_FractionXZone_i[i,"Y"]>87.5 && KO_seq_FractionXZone_i[i,"Somme"]>250 ) { KO_seq_FractionXZone_i[i,"TXT"] <- KO_seq_FractionXZone_i[i,"Kegg.Onthology"]}
          if (KO_seq_FractionXZone_i[i,"X"]>87.5 && KO_seq_FractionXZone_i[i,"Somme"]>250 ) { KO_seq_FractionXZone_i[i,"TXT"] <- KO_seq_FractionXZone_i[i,"Kegg.Onthology"]}
          if (KO_seq_FractionXZone_i[i,"Y"]>80 && KO_seq_FractionXZone_i[i,"Somme"]>500 ) { KO_seq_FractionXZone_i[i,"TXT"] <- KO_seq_FractionXZone_i[i,"Kegg.Onthology"]}
          if (KO_seq_FractionXZone_i[i,"X"]>80 && KO_seq_FractionXZone_i[i,"Somme"]>500 ) { KO_seq_FractionXZone_i[i,"TXT"] <- KO_seq_FractionXZone_i[i,"Kegg.Onthology"]}
          if (KO_seq_FractionXZone_i[i,"Y"]>70 && KO_seq_FractionXZone_i[i,"Somme"]>1000 ) { KO_seq_FractionXZone_i[i,"TXT"] <- KO_seq_FractionXZone_i[i,"Kegg.Onthology"]}
          if (KO_seq_FractionXZone_i[i,"X"]>70 && KO_seq_FractionXZone_i[i,"Somme"]>1000 ) { KO_seq_FractionXZone_i[i,"TXT"] <- KO_seq_FractionXZone_i[i,"Kegg.Onthology"]}
          
          if (KO_seq_FractionXZone_i[i,"Y"]<5 && KO_seq_FractionXZone_i[i,"Somme"]>100 ) { KO_seq_FractionXZone_i[i,"TXT"] <- KO_seq_FractionXZone_i[i,"Kegg.Onthology"]}
          if (KO_seq_FractionXZone_i[i,"X"]<5 && KO_seq_FractionXZone_i[i,"Somme"]>100 ) { KO_seq_FractionXZone_i[i,"TXT"] <- KO_seq_FractionXZone_i[i,"Kegg.Onthology"]}
          if (KO_seq_FractionXZone_i[i,"Y"]<12.5 && KO_seq_FractionXZone_i[i,"Somme"]>250 ) { KO_seq_FractionXZone_i[i,"TXT"] <- KO_seq_FractionXZone_i[i,"Kegg.Onthology"]}
          if (KO_seq_FractionXZone_i[i,"X"]<12.5 && KO_seq_FractionXZone_i[i,"Somme"]>250 ) { KO_seq_FractionXZone_i[i,"TXT"] <- KO_seq_FractionXZone_i[i,"Kegg.Onthology"]}
          if (KO_seq_FractionXZone_i[i,"Y"]<20 && KO_seq_FractionXZone_i[i,"Somme"]>500 ) { KO_seq_FractionXZone_i[i,"TXT"] <- KO_seq_FractionXZone_i[i,"Kegg.Onthology"]}
          if (KO_seq_FractionXZone_i[i,"X"]<20 && KO_seq_FractionXZone_i[i,"Somme"]>500 ) { KO_seq_FractionXZone_i[i,"TXT"] <- KO_seq_FractionXZone_i[i,"Kegg.Onthology"]}
          if (KO_seq_FractionXZone_i[i,"Y"]<30 && KO_seq_FractionXZone_i[i,"Somme"]>1000 ) { KO_seq_FractionXZone_i[i,"TXT"] <- KO_seq_FractionXZone_i[i,"Kegg.Onthology"]}
          if (KO_seq_FractionXZone_i[i,"X"]<30 && KO_seq_FractionXZone_i[i,"Somme"]>1000 ) { KO_seq_FractionXZone_i[i,"TXT"] <- KO_seq_FractionXZone_i[i,"Kegg.Onthology"]}
          
        }
        
        # Plot #2
        Quaternary_plot <- ggplot(KO_seq_FractionXZone_i, mapping = aes(y= Y, x = X)) + theme_unique_art() + 
          #annotate("rect", xmin = 5, xmax = 95, ymin = 5, ymax = 95, alpha = .2) +
          annotate("rect", xmin = 0, xmax = 100, ymin = 0, ymax = 100,
                   color = "black", fill = "#00000000") + 
          annotate("rect", xmin = 5, xmax = 95, ymin = 5, ymax = 95, color = "darkgrey", fill = "#00000000",linetype = "dashed") +
          annotate("rect", xmin = 12.5, xmax = 87.5, ymin = 12.5, ymax = 87.5, color = "grey", fill = "#00000000",linetype = "dashed") +
          annotate("rect", xmin = 20, xmax = 80, ymin = 20, ymax = 80, color = "lightgrey", fill = "#00000000",linetype = "dashed") +
          annotate("rect", xmin = 30, xmax = 70, ymin = 30, ymax = 70, color = "lightblue", fill = "#00000000",linetype = "dashed") +
          theme(legend.text=element_text(size=10),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "right") + 
          geom_point(aes(size=Somme), shape=21, stroke = 1, color = "#B71C1CFF") +
          scale_size(range = c(2,10),breaks = c(100, 250, 500, 1000, 2000)) +
          labs(color="Metabolisms/Cellular processes",fill="Metabolisms/Cellular processes", size = "Total abundance (RPKM)") + 
          geom_hline(yintercept=50, linetype="dashed", color = "darkgrey") +
          geom_vline(xintercept=50, linetype="dashed", color = "darkgrey") +
          annotate("text", x = -5, y = 50, label = "Large", angle=90, fontface = "bold") +
          annotate("text", x = 105, y = 50, label = "Small", angle = -90, fontface = "bold") +
          annotate("text", x = 50, y = -5, label = "Mixolimnion", fontface = "bold") +
          annotate("text", x = 50, y = 105, label = "Monimolimnion", fontface = "bold") +
          geom_label_repel(data = KO_seq_FractionXZone_i %>% filter(is.na(TXT) == FALSE),
                          size = 3.2,aes(label = TXT), 
                          box.padding = unit(0.45, "lines"),
                          max.overlaps = 1000) +
          guides(y.sec = "axis", size = guide_legend(title.position = "top",title.hjust = 0.5,ncol=1)) + 
          facet_wrap(nbGroup~.)
        print(Quaternary_plot)
        ggsave("Quaterneryplot/Quaterneryplot_abundance_Total_H.svg", device = "svg", width = 12, height = 12)
        write.table(KO_seq_FractionXZone_i,"Table/KO2Quaterneryplot.tsv",row.names = FALSE,col.names = TRUE,quote = FALSE, sep ="\t")
        # Table 2
        KO_seq_FractionXZone_j <- KO_seq_FractionXZone_i %>% filter(is.na(TXT)==FALSE) %>% arrange(nbGroup)
        write.table(KO_seq_FractionXZone_j,"Table/KO2Quaterneryplot_TXT.tsv",row.names = FALSE,col.names = TRUE,quote = FALSE, sep ="\t")
        
#
    # PCoA for Quaternary -----------------------------------------------------
      # ZoneXFraction -----------------------------------------------------------
         KO_seq_FractionXZone <- dplyr::left_join(FractionXZone_seq,tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select("Row.names","nbGroup","Kegg.Onthology","Last"),by="Row.names")
        #KO_seq_FractionXZone <- KO_seq_FractionXZone %>% filter(Last != "Dinobryon")
        #KO_seq_FractionXZone <- KO_seq_FractionXZone %>% filter(Last != "Ciliophora")
        KO_seq_FractionXZone <- separate_rows(KO_seq_FractionXZone,"Kegg.Onthology",sep = ",") # Duplicate row by KO 
        KO_seq_FractionXZone <- KO_seq_FractionXZone %>% dplyr::select(-"Row.names",-"Last") %>% group_by(Kegg.Onthology,nbGroup) %>% summarise_all(sum)
      # Start
        KO_seq_FractionXZone_i <- KO_seq_FractionXZone %>% filter(is.na(nbGroup)==FALSE) %>% filter(is.na(Kegg.Onthology)==FALSE) %>% filter(!nbGroup %in% c("Lichen","Other multicellular Fungi", "Embryophyceae","Unaffiliated"))
        KO_seq_FractionXZone_i <- reshape2::melt(KO_seq_FractionXZone_i,id.vars = c("nbGroup","Kegg.Onthology"))
        KO_seq_FractionXZone_i$CDT <- paste(KO_seq_FractionXZone_i$nbGroup,KO_seq_FractionXZone_i$variable,sep=" in ")
        dataCDT <- reshape2::dcast(data = KO_seq_FractionXZone_i,formula = Kegg.Onthology~CDT,fun.aggregate = sum,value.var = "value")
        row.names(dataCDT) <- dataCDT$Kegg.Onthology ; dataCDT <- dataCDT %>% select(-Kegg.Onthology)
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
        coord <- separate(coord,"variable",c("Fraction","Zone"),sep="X")
        # ggplot2
        palet_tree_grp <- c("#ad0635FF","#d84040FF","#4080c0FF","#60d0a0FF","#70a830FF","#98d048FF","#d7c11cFF","#40a078FF","#333b43FF")
        pcoa_plot_a <- ggplot(coord) + theme_unique_art() + 
          theme(legend.text=element_text(size=10),legend.position = "right") + 
          geom_point(aes(y= `PCoA 3`, x = `PCoA 1`,color = nbGroup, fill =Fraction, shape = Zone), size=3, stroke =2) + 
          #geom_circle(aes(x0 = `PCoA 1`, y0 = `PCoA 3`, fill= Fraction, r=0.005), colour = "black")  + 
          scale_shape_manual(values = c(21, 24)) +
          scale_color_manual(values = palet_tree_grp) + 
          scale_fill_manual(values = c("white","black")) + 
          labs(title="Abundance (RPKM)",color="Functional groups", shape = "Zone", x= paste0("Dim 1 [",Dim1Seq,"]"), y=paste0("Dim3 [",Dim3Seq,"]")) +
          theme(axis.title.x=element_blank()) +
          guides(fill = guide_legend(override.aes = list(shape = 21))) +
          scale_x_continuous(position = 'top')
        LEGEND_a <- get_legend(pcoa_plot_a)
        pcoa_plot_a <- pcoa_plot_a + theme(legend.position = "none")
        pcoa_plot_b <- ggplot(coord) + theme_unique_art() + 
          theme(legend.text=element_text(size=10),legend.position = "right") + 
          geom_point(aes(y= `PCoA 2`, x = `PCoA 1`,color = nbGroup, fill = Fraction, shape = Zone), size=3, stroke =2) + 
          theme(legend.position = "none") + 
          #geom_circle(aes(x0 = `PCoA 1`, y0 = `PCoA 2`, fill= Fraction, r=0.005), colour = "black") + 
          scale_shape_manual(values = c(21, 24)) +
          scale_color_manual(values = palet_tree_grp) + 
          scale_fill_manual(values = c("white","black")) + 
          labs(color="Functional groups", shape = "Zone", x= paste0("Dim1 [",Dim1Seq,"]"), y=paste0("Dim2 [",Dim2Seq,"]"))
      # Présence/absence
        dataCDT[dataCDT>0] <- 1
        pcoa <- ape::pcoa(vegdist(t(dataCDT), method = "bray"))
        coord <- pcoa$vectors[,1:3] ; colnames(coord) <- c("PCoA 1","PCoA 2","PCoA 3")
        Dim1rich <- paste0(round(pcoa$values[1,"Relative_eig"],3)*100,"%")
        Dim2rich <- paste0(round(pcoa$values[2,"Relative_eig"],3)*100,"%")
        Dim3rich <- paste0(round(pcoa$values[3,"Relative_eig"],3)*100,"%")
        coord <- as.data.frame(coord)
        coord$CDT <- row.names(coord)
        coord[,"nbGroup"] <- sapply(strsplit(as.character(coord$CDT)," in "), `[`, 1)
        coord[,"variable"] <- sapply(strsplit(as.character(coord$CDT)," in "), `[`, 2)
        coord <- separate(coord,"variable",c("Fraction","Zone"),sep="X")
        # ggplot2
        palet_tree_grp <- c("#ad0635FF","#d84040FF","#4080c0FF","#60d0a0FF","#70a830FF","#98d048FF","#d7c11cFF","#40a078FF","#333b43FF")
        pcoa_plot_c <- ggplot(coord) + theme_unique_art() + 
          theme(legend.text=element_text(size=10),legend.position = "right") + 
          geom_point(aes(y= `PCoA 3`, x = `PCoA 1`,color = nbGroup, fill=Fraction, shape = Zone), size=3, stroke = 2) + 
          #geom_circle(aes(x0 = `PCoA 1`, y0 = `PCoA 3`, fill= Fraction, r=0.005), colour = "black") +
          scale_shape_manual(values = c(21, 24)) +
          scale_color_manual(values = palet_tree_grp) + 
          scale_fill_manual(values = c("white","black")) + 
          labs(title="Presence/Absence",color="Functional groups", shape = "Zone", x= paste0("Dim 1 [",Dim1rich,"]"), y=paste0("Dim3 [",Dim3rich,"]")) +
          theme(axis.title.x=element_blank()) +
          scale_x_continuous(position = 'top') + scale_y_continuous(position = 'right') + 
          theme(legend.position = "none")
        pcoa_plot_d <- ggplot(coord) + theme_unique_art() + 
          theme(legend.text=element_text(size=10),legend.position = "right") + 
          geom_point(aes(y= `PCoA 2`, x = `PCoA 1`,color = nbGroup, fill=Fraction, shape = Zone), size=3, stroke = 2) + 
          theme(legend.position = "none") + 
          #geom_circle(aes(x0 = `PCoA 1`, y0 = `PCoA 2`, fill= Fraction, r=0.005), colour = "black") + 
          scale_shape_manual(values = c(21, 24)) +
          scale_color_manual(values = palet_tree_grp) + 
          scale_fill_manual(values = c("white","black")) + 
          labs(color="Functional groups", shape = "Zone", x= paste0("Dim1 [",Dim1rich,"]"), y=paste0("Dim2 [",Dim2rich,"]")) + 
          scale_y_continuous(position = 'right')
        
      # Coplot
        svglite("Quaterneryplot/GRPandCDT.svg",width = 12.00,height = 9.00)
        b_plot <- plot_grid(pcoa_plot_a,pcoa_plot_c,NULL,NULL,NULL,LEGEND_a, pcoa_plot_b, pcoa_plot_d,NULL, ncol = 3, nrow = 3, rel_widths = c(3,3,2),rel_heights = c(3,0,3))
        print(b_plot)
        dev.off()
      # ZoneXFractionXCycle -----------------------------------------------------------
        KO_seq_FractionXZoneXCycle <- dplyr::left_join(CycleXFractionXZone_seq,tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select("Row.names","nbGroup","Kegg.Onthology","Last"),by="Row.names")
        #KO_seq_FractionXZoneXCycle <- KO_seq_FractionXZoneXCycle %>% filter(Last != "Dinobryon")
        #KO_seq_FractionXZoneXCycle <- KO_seq_FractionXZoneXCycle %>% filter(Last != "Ciliophora")
        KO_seq_FractionXZoneXCycle <- separate_rows(KO_seq_FractionXZoneXCycle,"Kegg.Onthology",sep = ",") # Duplicate row by KO 
        KO_seq_FractionXZoneXCycle <- KO_seq_FractionXZoneXCycle %>% dplyr::select(-"Row.names",-"Last") %>% group_by(Kegg.Onthology,nbGroup) %>% summarise_all(sum)
        # Start
        KO_seq_FractionXZoneXCycle_i <- KO_seq_FractionXZoneXCycle %>% filter(is.na(nbGroup)==FALSE) %>% filter(is.na(Kegg.Onthology)==FALSE) %>% filter(!nbGroup %in% c("Lichen","Other multicellular Fungi", "Embryophyceae","Unaffiliated"))
        KO_seq_FractionXZoneXCycle_i <- reshape2::melt(KO_seq_FractionXZoneXCycle_i,id.vars = c("nbGroup","Kegg.Onthology"))
        KO_seq_FractionXZoneXCycle_i$CDT <- paste(KO_seq_FractionXZoneXCycle_i$nbGroup,KO_seq_FractionXZoneXCycle_i$variable,sep=" in ")
        dataCDT <- reshape2::dcast(data = KO_seq_FractionXZoneXCycle_i,formula = Kegg.Onthology~CDT,fun.aggregate = sum,value.var = "value")
        row.names(dataCDT) <- dataCDT$Kegg.Onthology ; dataCDT <- dataCDT %>% select(-Kegg.Onthology)
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
        coord <- separate(coord,"variable",c("Fraction","Zone","Cycle"),sep="X")
        # ggplot2
        palet_tree_grp <- c("#ad0635FF","#d84040FF","#4080c0FF","#60d0a0FF","#70a830FF","#98d048FF","#d7c11cFF","#40a078FF","#333b43FF")
        pcoa_plot_a <- ggplot(coord) + theme_unique_art() + 
          theme(legend.text=element_text(size=10),legend.position = "right") + 
          geom_point(aes(y= `PCoA 3`, x = `PCoA 1`,color = nbGroup, fill = Cycle, shape = Zone, size = Fraction), stroke =2) + 
          #geom_circle(aes(x0 = `PCoA 1`, y0 = `PCoA 3`, fill= Fraction, r=0.005), colour = "black")  + 
          scale_shape_manual(values = c(21, 24)) +
          scale_size_manual(values = c(4, 2)) +
          scale_color_manual(values = palet_tree_grp) + 
          scale_fill_manual(values = c("white","black")) + 
          labs(title="Abundance (RPKM)",color="Functional groups", shape = "Zone (shape)", size = "Fraction (size)", fill = "Cycle (fill)", x= paste0("Dim 1 [",Dim1Seq,"]"), y=paste0("Dim3 [",Dim3Seq,"]")) +
          theme(axis.title.x=element_blank()) +
          guides(fill = guide_legend(override.aes = list(shape = 21, size = 3)), shape = guide_legend(override.aes = list(size = 3))) +
          scale_x_continuous(position = 'top') +
          scale_y_continuous(position = 'right')
        LEGEND_a <- get_legend(pcoa_plot_a)
        pcoa_plot_a <- pcoa_plot_a + theme(legend.position = "none")
        pcoa_plot_b <- ggplot(coord) + theme_unique_art() + 
          theme(legend.text=element_text(size=10),legend.position = "right") + 
          geom_point(aes(y= `PCoA 2`, x = `PCoA 1`,color = nbGroup, fill = Cycle, shape = Zone, size = Fraction), stroke =2) + 
          theme(legend.position = "none") + 
          #geom_circle(aes(x0 = `PCoA 1`, y0 = `PCoA 2`, fill= Fraction, r=0.005), colour = "black") + 
          scale_shape_manual(values = c(21, 24)) +
          scale_size_manual(values = c(4, 2)) +
          scale_color_manual(values = palet_tree_grp) + 
          scale_fill_manual(values = c("white","black")) + 
          labs(color="Functional groups", shape = "Zone", x= paste0("Dim1 [",Dim1Seq,"]"), y=paste0("Dim2 [",Dim2Seq,"]")) +
          scale_y_continuous(position = 'right')
        # Présence/absence
        dataCDT[dataCDT>0] <- 1
        pcoa <- ape::pcoa(vegdist(t(dataCDT), method = "bray"))
        coord <- pcoa$vectors[,1:3] ; colnames(coord) <- c("PCoA 1","PCoA 2","PCoA 3")
        Dim1rich <- paste0(round(pcoa$values[1,"Relative_eig"],3)*100,"%")
        Dim2rich <- paste0(round(pcoa$values[2,"Relative_eig"],3)*100,"%")
        Dim3rich <- paste0(round(pcoa$values[3,"Relative_eig"],3)*100,"%")
        coord <- as.data.frame(coord)
        coord$CDT <- row.names(coord)
        coord[,"nbGroup"] <- sapply(strsplit(as.character(coord$CDT)," in "), `[`, 1)
        coord[,"variable"] <- sapply(strsplit(as.character(coord$CDT)," in "), `[`, 2)
        coord <- separate(coord,"variable",c("Fraction","Zone","Cycle"),sep="X")
        # ggplot2
        palet_tree_grp <- c("#ad0635FF","#d84040FF","#4080c0FF","#60d0a0FF","#70a830FF","#98d048FF","#d7c11cFF","#40a078FF","#333b43FF")
        pcoa_plot_c <- ggplot(coord) + theme_unique_art() + 
          theme(legend.text=element_text(size=10),legend.position = "right") + 
          geom_point(aes(y= `PCoA 3`, x = `PCoA 1`,color = nbGroup, fill = Cycle, shape = Zone, size = Fraction), stroke = 2) + 
          #geom_circle(aes(x0 = `PCoA 1`, y0 = `PCoA 3`, fill= Fraction, r=0.005), colour = "black") +
          scale_shape_manual(values = c(21, 24)) +
          scale_size_manual(values = c(4, 2)) +
          scale_color_manual(values = palet_tree_grp) + 
          scale_fill_manual(values = c("white","black")) + 
          labs(title="Presence/Absence",color="Functional groups", shape = "Zone", x= paste0("Dim 1 [",Dim1rich,"]"), y=paste0("Dim3 [",Dim3rich,"]")) +
          theme(axis.title.x=element_blank()) +
          scale_x_continuous(position = 'top') + scale_y_continuous(position = 'left') + 
          theme(legend.position = "none")
        pcoa_plot_d <- ggplot(coord) + theme_unique_art() + 
          theme(legend.text=element_text(size=10),legend.position = "right") + 
          geom_point(aes(y= `PCoA 2`, x = `PCoA 1`,color = nbGroup, fill = Cycle, shape = Zone, size = Fraction), stroke = 2) + 
          theme(legend.position = "none") + 
          #geom_circle(aes(x0 = `PCoA 1`, y0 = `PCoA 2`, fill= Fraction, r=0.005), colour = "black") + 
          scale_shape_manual(values = c(21, 24)) +
          scale_size_manual(values = c(4, 2)) +
          scale_color_manual(values = palet_tree_grp) + 
          scale_fill_manual(values = c("white","black")) + 
          labs(color="Functional groups", shape = "Zone", x= paste0("Dim1 [",Dim1rich,"]"), y=paste0("Dim2 [",Dim2rich,"]")) + 
          scale_y_continuous(position = 'left')
        
        # Coplot
        svglite("Quaterneryplot/GRPandCDTS.svg",width = 12.00,height = 9.00)
        b_plot <- plot_grid(pcoa_plot_c, pcoa_plot_a,NULL,NULL,NULL,LEGEND_a, pcoa_plot_d, pcoa_plot_b, NULL, ncol = 3, nrow = 3, rel_widths = c(3,3,2),rel_heights = c(3,0,3))
        print(b_plot)
        dev.off()
        
    # Taxo ------------------------------------------------------------
  Taxo_seq_FractionXZone <- dplyr::left_join(FractionXZone_seq,tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select("Row.names","nbGroup","Supergroup","Division"),by="Row.names")
  Taxo_seq_FractionXZone_i <- Taxo_seq_FractionXZone %>% dplyr::select(-"Row.names") %>% group_by(Supergroup,Division,nbGroup) %>% summarise_all(sum)
  #
  Taxo_seq_FractionXZone_i <- reshape2::melt(Taxo_seq_FractionXZone_i,id.vars = c("nbGroup","Supergroup","Division"))
  Taxo_seq_FractionXZone_i <- Taxo_seq_FractionXZone_i %>% mutate(Division = ifelse(is.na(Division) == TRUE, paste0('Unaffiliated_',Supergroup),Division))
  Taxo_seq_FractionXZone_i[Taxo_seq_FractionXZone_i=="Unaffiliated_NA"] <- NA
  Taxo_seq_FractionXZone_i <- separate(Taxo_seq_FractionXZone_i,"variable",c("Fraction","Zone"),sep="X")
  Taxo_seq_FractionXZone_i[,"Division"][is.na(Taxo_seq_FractionXZone_i[,"Division"]) == TRUE] <- "Unaffiliated"
  Taxo_seq_FractionXZone_i[,"nbGroup"][is.na(Taxo_seq_FractionXZone_i[,"nbGroup"]) == TRUE] <- "Unassigned"
  
  Taxo_seq_FractionXZone_i <- Taxo_seq_FractionXZone_i %>% filter(Division != "Unaffiliated") %>% filter(!nbGroup %in% c("Unassigned","Lichen","Other multicellular Fungi", "Embryophyceae","Unaffiliated"))
  for (i in row.names(Taxo_seq_FractionXZone_i)) {
    Fractionx <- Taxo_seq_FractionXZone_i[i,"Fraction"]
    Zonex <- Taxo_seq_FractionXZone_i[i,"Zone"]
    nbGroupx <- Taxo_seq_FractionXZone_i[i,"nbGroup"]
    Taxo_seq_FractionXZone_i[i,"Sum"] <- sum(Taxo_seq_FractionXZone_i %>% filter(Fraction == Fractionx) %>% filter(Zone == Zonex) %>% filter(nbGroup == nbGroupx) %>% dplyr::select(value))
    Taxo_seq_FractionXZone_i[i,"Relative"] <- Taxo_seq_FractionXZone_i[i,"Sum"] *100 / sum(Taxo_seq_FractionXZone_i %>% filter(Fraction == Fractionx) %>% filter(Zone == Zonex) %>% dplyr::select(value))
    Taxo_seq_FractionXZone_i[i,"Proportion"] <- Taxo_seq_FractionXZone_i[i,"value"]*100/Taxo_seq_FractionXZone_i[i,"Sum"]}
  Taxo_seq_FractionXZone_i[,"nbGroup"] <- sapply(strsplit(as.character(Taxo_seq_FractionXZone_i$nbGroup),": "), `[`, 2)
  color_tax_table_order <- unique(Taxo_seq_FractionXZone_i$Division)
  
  #plot
  Total_ZoneXFraction_Sequence_fig <- ggplot(Taxo_seq_FractionXZone_i, mapping = aes(y= Proportion, x = Zone, fill = factor(Division, levels = color_tax_table_order), group = factor(Division, levels = color_tax_table_order)), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity",alpha=0.95,color="black") +
    scale_fill_manual(values = paletspe_i, na.value = NA) + facet_grid(Fraction~nbGroup,scales="free") + theme_unique_art() +
    labs(x="Zones",y="Unigenes %",fill="Division") + guides(fill=guide_legend(ncol=1)) +
    geom_label(aes(y = 106,label = paste0(round(Sum,1)," TPM","\n","(",round(Relative,1)," %)")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B99")
  print(Total_ZoneXFraction_Sequence_fig)
  ggsave("Hist/Grp_ZoneXFraction_Taxonomy_abundance.svg", device = "svg", width = 18, height = 10)
  #    


# Conditions --------------------------------------------------------------
  KO_seq_FractionXZoneXCycle <- dplyr::left_join(CycleXFractionXZone_seq,tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select("Row.names","nbGroup","Kegg.Onthology"),by="Row.names")
  KO_seq_FractionXZoneXCycle <- separate_rows(KO_seq_FractionXZoneXCycle,"Kegg.Onthology",sep = ",") # Duplicate row by KO 
  #
  # Condition Know fct°---------------------------------------------------------------
  KO_seq_FractionXZoneXCycle_i <- left_join(KO_seq_FractionXZoneXCycle,ko2Met_Table_Final_i,by=c("Kegg.Onthology"="KO_id"))
  colnames(KO_seq_FractionXZoneXCycle_i)[colnames(KO_seq_FractionXZoneXCycle_i) == 'Pathway/Gene'] <- 'Label'
  KO_seq_FractionXZoneXCycle_i[,"Label"] <- sapply(strsplit(as.character(KO_seq_FractionXZoneXCycle_i$Label),"_-_"), `[`, 1)
  KO_seq_FractionXZoneXCycle_i[,"Label"] <- sapply(str_replace_all(as.character(KO_seq_FractionXZoneXCycle_i$Label),"_"," "), `[`, 1)
  KO_seq_FractionXZoneXCycle_i <- KO_seq_FractionXZoneXCycle_i %>% filter(is.na(nbGroup)==FALSE) %>% filter(is.na(Kegg.Onthology)==FALSE) %>% filter(is.na(Label)==FALSE)
  #
  KO_seq_FractionXZoneXCycle_i <- KO_seq_FractionXZoneXCycle_i %>% dplyr::select(-"Row.names",-"Kegg.Onthology",-"Label") %>% group_by(nbGroup) %>% summarise_all(sum)
  KO_seq_FractionXZoneXCycle_i <- KO_seq_FractionXZoneXCycle_i %>% filter(is.na(nbGroup)==FALSE) %>% filter(!nbGroup %in% c("Lichen","Other multicellular Fungi", "Embryophyceae","Unaffiliated"))
  #
  KO_seq_FractionXZoneXCycle_i <- reshape2::melt(KO_seq_FractionXZoneXCycle_i,id.vars = c("nbGroup"))
  KO_seq_FractionXZoneXCycle_i <- separate(KO_seq_FractionXZoneXCycle_i,"variable",c("Fraction","Zone","Cycle"),sep="X")
  for (i in row.names(KO_seq_FractionXZoneXCycle_i)) {
    Fractionx <- KO_seq_FractionXZoneXCycle_i[i,"Fraction"]
    Zonex <- KO_seq_FractionXZoneXCycle_i[i,"Zone"]
    Cyclex <- KO_seq_FractionXZoneXCycle_i[i,"Cycle"]
    nbGroupx <- KO_seq_FractionXZoneXCycle_i[i,"nbGroup"]
    KO_seq_FractionXZoneXCycle_i[i,"Sum"] <- sum(KO_seq_FractionXZoneXCycle_i %>% filter(Fraction == Fractionx) %>% filter(Zone == Zonex) %>% filter(Cycle == Cyclex) %>% dplyr::select(value))
    KO_seq_FractionXZoneXCycle_i[i,"Relative"] <- KO_seq_FractionXZoneXCycle_i[i,"value"] *100 / sum(KO_seq_FractionXZoneXCycle_i %>% filter(Fraction == Fractionx) %>% filter(Zone == Zonex) %>% filter(Cycle == Cyclex) %>% dplyr::select(value))
    KO_seq_FractionXZoneXCycle_i[i,"Proportion"] <- KO_seq_FractionXZoneXCycle_i[i,"value"]*100/KO_seq_FractionXZoneXCycle_i[i,"Sum"]}
  #KO_seq_FractionXZoneXCycle_i[,"nbGroup"] <- sapply(strsplit(as.character(KO_seq_FractionXZoneXCycle_i$nbGroup),": "), `[`, 2)
  #Plot
  palet_tree_grp <- c("#ad0635FF","#d84040FF","#4080c0FF","#60d0a0FF","#70a830FF","#98d048FF","#d7c11cFF","#40a078FF","#333b43FF")
  KO_seq_FractionXZoneXCycle_i_fig <- ggplot(KO_seq_FractionXZoneXCycle_i, mapping = aes(y= value, x = Fraction, fill = nbGroup, group = nbGroup), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity",alpha=0.95,color="black") +
    facet_grid(Zone~Cycle,scales="free") + theme_unique_art() + 
    scale_fill_manual(values = palet_tree_grp) + 
    labs(x="Fraction",y="Abundance (RPKM)",fill="Functional groups") + guides(fill=guide_legend(ncol=1))
  print(KO_seq_FractionXZoneXCycle_i_fig)
  ggsave("Hist/Grp_conditions_KnowFunction_abundance.svg", device = "svg", width = 8, height = 8)
#
  # Condition All fct°---------------------------------------------------------------
  KO_seq_FractionXZoneXCycle_i <- KO_seq_FractionXZoneXCycle %>% dplyr::select(-"Row.names",-"Kegg.Onthology") %>% group_by(nbGroup) %>% summarise_all(sum)
  KO_seq_FractionXZoneXCycle_i <- KO_seq_FractionXZoneXCycle_i %>% filter(is.na(nbGroup)==FALSE) %>% filter(!nbGroup %in% c("Lichen","Other multicellular Fungi", "Embryophyceae","Unaffiliated"))
  #
  KO_seq_FractionXZoneXCycle_i <- reshape2::melt(KO_seq_FractionXZoneXCycle_i,id.vars = c("nbGroup"))
  KO_seq_FractionXZoneXCycle_i <- separate(KO_seq_FractionXZoneXCycle_i,"variable",c("Fraction","Zone","Cycle"),sep="X")
  for (i in row.names(KO_seq_FractionXZoneXCycle_i)) {
    Fractionx <- KO_seq_FractionXZoneXCycle_i[i,"Fraction"]
    Zonex <- KO_seq_FractionXZoneXCycle_i[i,"Zone"]
    Cyclex <- KO_seq_FractionXZoneXCycle_i[i,"Cycle"]
    nbGroupx <- KO_seq_FractionXZoneXCycle_i[i,"nbGroup"]
    KO_seq_FractionXZoneXCycle_i[i,"Sum"] <- sum(KO_seq_FractionXZoneXCycle_i %>% filter(Fraction == Fractionx) %>% filter(Zone == Zonex) %>% filter(Cycle == Cyclex) %>% dplyr::select(value))
    KO_seq_FractionXZoneXCycle_i[i,"Relative"] <- KO_seq_FractionXZoneXCycle_i[i,"value"] *100 / sum(KO_seq_FractionXZoneXCycle_i %>% filter(Fraction == Fractionx) %>% filter(Zone == Zonex) %>% filter(Cycle == Cyclex) %>% dplyr::select(value))
    KO_seq_FractionXZoneXCycle_i[i,"Proportion"] <- KO_seq_FractionXZoneXCycle_i[i,"value"]*100/KO_seq_FractionXZoneXCycle_i[i,"Sum"]}
  #KO_seq_FractionXZoneXCycle_i[,"nbGroup"] <- sapply(strsplit(as.character(KO_seq_FractionXZoneXCycle_i$nbGroup),": "), `[`, 2)
  #Plot
  palet_tree_grp <- c("#ad0635FF","#d84040FF","#4080c0FF","#60d0a0FF","#70a830FF","#98d048FF","#d7c11cFF","#40a078FF","#333b43FF")
  KO_seq_FractionXZoneXCycle_i_fig <- ggplot(KO_seq_FractionXZoneXCycle_i, mapping = aes(y= value, x = Fraction, fill = nbGroup, group = nbGroup), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity",alpha=0.95,color="black") +
    facet_grid(Zone~Cycle,scales="free") + theme_unique_art() + 
    scale_fill_manual(values = palet_tree_grp) + 
    labs(x="Fraction",y="Abundance (RPKM)",fill="Functional groups") + guides(fill=guide_legend(ncol=1))
  print(KO_seq_FractionXZoneXCycle_i_fig)
  ggsave("Hist/Grp_conditions_AllFunction_abundance.svg", device = "svg", width = 8, height = 8)
  #
  # Taxonomy All fct°---------------------------------------------------------------
  KO_seq_FractionXZoneXCycle <- dplyr::left_join(CycleXFractionXZone_seq,tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select("Row.names","nbGroup","Kegg.Onthology","Division"),by="Row.names")
  KO_seq_FractionXZoneXCycle <- separate_rows(KO_seq_FractionXZoneXCycle,"Kegg.Onthology",sep = ",") # Duplicate row by KO 
  #
  KO_seq_FractionXZoneXCycle_i <- KO_seq_FractionXZoneXCycle %>% dplyr::select(-"Row.names",-"Kegg.Onthology") %>% group_by(nbGroup,Division) %>% summarise_all(sum)
  KO_seq_FractionXZoneXCycle_i <- KO_seq_FractionXZoneXCycle_i %>% filter(is.na(nbGroup)==FALSE) %>% filter(!nbGroup %in% c("Lichen","Other multicellular Fungi", "Embryophyceae","Unaffiliated"))
  KO_seq_FractionXZoneXCycle_i[,"Division"][is.na(KO_seq_FractionXZoneXCycle_i[,"Division"])==TRUE] <- "Unaffiliated"
  #
  KO_seq_FractionXZoneXCycle_i <- reshape2::melt(KO_seq_FractionXZoneXCycle_i,id.vars = c("nbGroup","Division"))
  KO_seq_FractionXZoneXCycle_i <- separate(KO_seq_FractionXZoneXCycle_i,"variable",c("Fraction","Zone","Cycle"),sep="X")
  for (i in row.names(KO_seq_FractionXZoneXCycle_i)) {
    Fractionx <- KO_seq_FractionXZoneXCycle_i[i,"Fraction"]
    Zonex <- KO_seq_FractionXZoneXCycle_i[i,"Zone"]
    Cyclex <- KO_seq_FractionXZoneXCycle_i[i,"Cycle"]
    nbGroupx <- KO_seq_FractionXZoneXCycle_i[i,"nbGroup"]
    Divisionx <- KO_seq_FractionXZoneXCycle_i[i,"Divisionx"]
    KO_seq_FractionXZoneXCycle_i[i,"Sum"] <- sum(KO_seq_FractionXZoneXCycle_i %>% filter(Fraction == Fractionx) %>% filter(Zone == Zonex) %>% filter(Cycle == Cyclex) %>% filter(nbGroup == nbGroupx) %>% dplyr::select(value))/nrow(KO_seq_FractionXZoneXCycle_i %>% filter(Fraction == Fractionx) %>% filter(Zone == Zonex) %>% filter(Cycle == Cyclex) %>% filter(nbGroup == nbGroupx))
    KO_seq_FractionXZoneXCycle_i[i,"Relative"] <- KO_seq_FractionXZoneXCycle_i[i,"value"] *100 / sum(KO_seq_FractionXZoneXCycle_i %>% filter(Fraction == Fractionx) %>% filter(Zone == Zonex) %>% filter(Cycle == Cyclex) %>% dplyr::select(value))
    KO_seq_FractionXZoneXCycle_i[i,"Proportion"] <- KO_seq_FractionXZoneXCycle_i[i,"value"]*100/KO_seq_FractionXZoneXCycle_i[i,"Sum"]}
  KO_seq_FractionXZoneXCycle_i[,"nbGroup"] <- sapply(strsplit(as.character(KO_seq_FractionXZoneXCycle_i$nbGroup),": "), `[`, 2)
  #Plot
  color_tax_table_order <- unique(Taxo_seq_FractionXZone_i$Division)
  paletspe_i <- rep(c("#e53a35","#f44336","#ef5350",
                      "#303f9f",#"#303f9fFF",
                      "#3065c0","#3f88e5",
                      "#2e7d32",
                      "#fbc02c","#fdd835",
                      "#616161","#757575",
                      "#1876d2",
                      "#3897a7","#41acc1"),
                    4)
  KO_seq_FractionXZoneXCycle_i_fig <- ggplot() + 
    geom_bar(data = KO_seq_FractionXZoneXCycle_i,
             mapping = aes(y= value, x = Fraction, fill = factor(Division, levels = color_tax_table_order), group = nbGroup),
             stat="identity",alpha=0.95) + 
    geom_bar(data = KO_seq_FractionXZoneXCycle_i %>% select(-Division) %>% group_by(nbGroup, Zone, Fraction, Cycle) %>% summarise_all(sum), 
             mapping = aes(y= Sum, x = Fraction, fill = NA, color = nbGroup),
             stat="identity") + 
    geom_text(data = KO_seq_FractionXZoneXCycle_i %>% select(-Division) %>% group_by(nbGroup, Zone, Fraction, Cycle) %>% summarise_all(sum) %>% arrange(desc(nbGroup)),
              mapping = aes(y= Sum, x = Fraction, label = nbGroup), 
              size = 2.5, position = position_stack(vjust = 0.5),
              fontface = "bold") +
    facet_grid(Zone~Cycle,scales="free") + theme_unique_art() + 
    scale_fill_manual(values = paletspe_i, na.value = NA) + 
    scale_color_manual(values = rep("black",7)) + 
    labs(x="Fraction",y="Abundance (RPKM)",fill="Functional groups") + guides(fill=guide_legend(ncol=1)) +
    guides(color="none")
  print(KO_seq_FractionXZoneXCycle_i_fig)
  ggsave("Hist/GrpXDivisions_conditions_AllFunction_abundance.svg", device = "svg", width = 8, height = 10)
  
  
  #
  # Taxonomy All fct° without Cycle ---------------------------------------------------------------
  KO_seq_FractionXZone <- dplyr::left_join(FractionXZone_seq,tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select("Row.names","nbGroup","Kegg.Onthology","Division", "Supergroup"),by="Row.names")
  #KO_seq_FractionXZone <- separate_rows(KO_seq_FractionXZone,"Kegg.Onthology",sep = ",") # Duplicate row by KO 
  KO_seq_FractionXZone <- KO_seq_FractionXZone %>% filter(is.na(Kegg.Onthology)==FALSE)
  #
  KO_seq_FractionXZone <- KO_seq_FractionXZone %>% mutate(Division = ifelse(is.na(Division) == TRUE, paste0('Unaffiliated ',Supergroup),Division))
  KO_seq_FractionXZone[KO_seq_FractionXZone=="Unaffiliated_NA"] <- NA
  color_tax_table_order <- unique(KO_seq_FractionXZone %>% arrange(Supergroup,Division) %>% select(Division))
  #
  KO_seq_FractionXZone_i <- KO_seq_FractionXZone %>% dplyr::select(-"Row.names",-"Kegg.Onthology",-"Supergroup") %>% group_by(nbGroup,Division) %>% summarise_all(sum)
  KO_seq_FractionXZone_i[,"Division"][is.na(KO_seq_FractionXZone_i[,"Division"])==TRUE] <- "Unaffiliated"
  KO_seq_FractionXZone_i <- KO_seq_FractionXZone_i %>% filter(is.na(nbGroup)==FALSE) %>% filter(!nbGroup %in% c("Lichen","Other multicellular Fungi", "Embryophyceae","Unaffiliated"))  #
  KO_seq_FractionXZone_i <- reshape2::melt(KO_seq_FractionXZone_i,id.vars = c("nbGroup","Division"))
  KO_seq_FractionXZone_i <- separate(KO_seq_FractionXZone_i,"variable",c("Fraction","Zone"),sep="X")
  for (i in row.names(KO_seq_FractionXZone_i)) {
    Fractionx <- KO_seq_FractionXZone_i[i,"Fraction"]
    Zonex <- KO_seq_FractionXZone_i[i,"Zone"]
    nbGroupx <- KO_seq_FractionXZone_i[i,"nbGroup"]
    Divisionx <- KO_seq_FractionXZone_i[i,"Divisionx"]
    KO_seq_FractionXZone_i[i,"Sum"] <- sum(KO_seq_FractionXZone_i %>% filter(Fraction == Fractionx) %>% filter(Zone == Zonex) %>% filter(nbGroup == nbGroupx) %>% dplyr::select(value))/nrow(KO_seq_FractionXZone_i %>% filter(Fraction == Fractionx) %>% filter(Zone == Zonex) %>% filter(nbGroup == nbGroupx))
    KO_seq_FractionXZone_i[i,"Relative"] <- KO_seq_FractionXZone_i[i,"value"] *100 / sum(KO_seq_FractionXZone_i %>% filter(Fraction == Fractionx) %>% filter(Zone == Zonex) %>% dplyr::select(value))
    KO_seq_FractionXZone_i[i,"Proportion"] <- KO_seq_FractionXZone_i[i,"value"]*100/KO_seq_FractionXZone_i[i,"Sum"]}
  KO_seq_FractionXZone_i[,"nbGroup"] <- sapply(strsplit(as.character(KO_seq_FractionXZone_i$nbGroup),": "), `[`, 2)
  #Plot
  color_tax_table_order <- color_tax_table_order %>% filter(Division %in% KO_seq_FractionXZone_i$Division)
  color_tax_table_order <- c(color_tax_table_order$Division)
  paletspe_i <- rep(c("#b92b2c","#d43b3b","#f24d43",
                      "#8e2059","#c42863",
                      "#3e8dd1","#4bacf4",
                      "#2e7d32",
                      "#f5a520","#f8d835",#"#fdd835",
                      "#616161","#757575",
                      "#1876d2",
                      "#3897a7","#41acc1"),
                    4)
  KO_seq_FractionXZone_i_fig <- ggplot() + 
    geom_bar(data = KO_seq_FractionXZone_i,
             mapping = aes(y= value, x = Fraction, fill = factor(Division, levels = color_tax_table_order), group = nbGroup),
             stat="identity",alpha=1) + 
    geom_bar(data = KO_seq_FractionXZone_i %>% select(-Division) %>% group_by(nbGroup, Zone, Fraction) %>% summarise_all(sum), 
             mapping = aes(y= Sum, x = Fraction, color = nbGroup), fill = NA,
             stat="identity") + 
    geom_text(data = KO_seq_FractionXZone_i %>% select(-Division) %>% group_by(nbGroup, Zone, Fraction) %>% summarise_all(sum) %>% arrange(desc(nbGroup)),
              mapping = aes(y= Sum, x = Fraction, label = nbGroup), 
              size = 2.5, position = position_stack(vjust = 0.5),
              fontface = "bold") +
    facet_grid(.~Zone,scales="free") + theme_unique_art() + 
    scale_fill_manual(values = paletspe_i, na.value = NA) + 
    scale_color_manual(values = rep("black",7)) + 
    labs(x="Fraction",y="Abundance (RPKM)",fill="Divisions") + guides(fill=guide_legend(ncol=1)) +
    guides(color="none")
  print(KO_seq_FractionXZone_i_fig)
  ggsave("Hist/GrpXDivisions_conditions_AllFunction*_abundance.svg", device = "svg", width = 8, height = 6)
  write.table(KO_seq_FractionXZone_i,"Table/FunctionalGroupContribXCDTS.tsv",row.names = FALSE,col.names = TRUE,quote = FALSE, sep ="\t")
  
  #
# Table unigene-Taxo-GRP --------------------------------------------------
  x <- tableVinput_abundance_normalize_annots_tax_Funct 
  x$Somme <- rowSums(x[,all_of(MetaT_Summary$Echantillon)])
  x$LargeXMixolimnion <- rowSums(x[,all_of(LargeXMixolimnion)])
  x$SmallXMixolimnion <- rowSums(x[,all_of(SmallXMixolimnion)])
  x$LargeXMonimolimnion <- rowSums(x[,all_of(LargeXMonimolimnion)])
  x$SmallXMonimolimnion <- rowSums(x[,all_of(SmallXMonimolimnion)])
  x <- x %>% dplyr::select(-all_of(MetaT_Summary$Echantillon),-ko,-ko_BestHitNotSignif,-PR2_taxonomy)
  x <- x %>% filter(is.na(Kegg.Onthology) == FALSE) %>% filter(is.na(nbGroup) == FALSE)
  x <- separate_rows(x,"Kegg.Onthology",sep = ",") # Duplicate row by KO 
  x <- x %>% filter(Kegg.Onthology %in% all_of(KO_seq_FractionXZone_i$TXT))
  x <- x %>% filter(Somme > 10) %>% arrange(Somme)
  write.table(x,"Table/GRPXFUNXTAX.tsv",row.names = FALSE,col.names = TRUE,quote = FALSE, sep ="\t")
#
# Figure Metabolism -------------------------------------------------------
  # Relative ----------------------------------------------------------------


    # Color
  paletspe <- c("#b71c1d","#c62828","#d32f2e","#e53935","#f24336","#ef5350",
                "#880e4f","#ad1557","#c2185a",
                "#4a148c","#6a1c9a",
                "#33691d","#558b2f","#689f38","#7cb342","#8bc34a","#9ccc65",
                "#164c40","#21695c","#28796c",
                "#2447a1","#3065c0","#3775d2","#3f88e5","#4496f3","#48a5f5",
                "#f6c02d","#f8d835","#f9eb3b","#faee58",#"#fbf176",
                "#4527a0",
                "#1b237e","#283593","#303e9f","#3949ab",
                #"#e64a18","#f2511e","#f25823","#f37043",
                "#2f838f","#3897a7","#41acc1","#48bcd4","#4cc6da",
                #"#424242","#616161","#757575","#9e9e9e","#bdbdbd",
                "#bdbdbd")
  # Taxa
  uniqueTax <- unique(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select("Supergroup","Division") %>% arrange(Supergroup))
  uniqueTax <- uniqueTax %>% mutate(Division = ifelse(is.na(Division) == TRUE, paste0('Unaffiliated_',Supergroup),Division))
  uniqueTax[uniqueTax=="Unaffiliated_NA"] <- "Unaffiliated"
  uniqueTax <- uniqueTax %>% mutate(Division = ifelse(grepl(Division,pattern = "_X") == TRUE, paste0('Unaffiliated_',Supergroup),Division)) %>% distinct(Division,Supergroup) %>% arrange(Supergroup,Division)
  
  uniqueTax$color <- unique(paletspe)
  uniqueTax <- rbind(uniqueTax,data.frame(Supergroup="ZZZ",Division="Other<0.005%",color="#f36f02"))
  # launch loop
  Metabo_seq_FractionXZone <- dplyr::left_join(FractionXZone_seq,tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select("Row.names","nbGroup","Kegg.Onthology","Supergroup","Division"),by="Row.names")
  Metabo_seq_FractionXZone <- separate_rows(Metabo_seq_FractionXZone,"Kegg.Onthology",sep = ",") # Duplicate row by KO 
  Metabo_seq_FractionXZone <- Metabo_seq_FractionXZone %>% dplyr::select(-"Row.names") %>% group_by(Kegg.Onthology,nbGroup,Supergroup,Division) %>% summarise_all(sum)
  Metabo_seq_FractionXZone <- left_join(Metabo_seq_FractionXZone,ko2Met_Table_Final_i,by=c("Kegg.Onthology"="KO_id"))
  colnames(Metabo_seq_FractionXZone)[colnames(Metabo_seq_FractionXZone) == 'Pathway/Gene'] <- 'Label'
  Metabo_seq_FractionXZone[,"Label"] <- sapply(strsplit(as.character(Metabo_seq_FractionXZone$Label),"_-_"), `[`, 1)
  Metabo_seq_FractionXZone[,"Label"] <- sapply(str_replace_all(as.character(Metabo_seq_FractionXZone$Label),"_"," "), `[`, 1)
  Metabo_seq_FractionXZone <- as.data.frame(Metabo_seq_FractionXZone)
  for (i in unique(Metabo_seq_FractionXZone$Label)) { if (is.na(i) == FALSE) { print (i)
    j <- str_replace_all(i," ","_")
    # Select i
    Metabo_seq_FractionXZone_i <- Metabo_seq_FractionXZone %>% filter(Label==i)
    Metabo_seq_FractionXZone_i <- Metabo_seq_FractionXZone_i %>% dplyr::select(-"Kegg.Onthology") %>% group_by(nbGroup,Label,Supergroup,Division) %>% summarise_all(sum)
    Metabo_seq_FractionXZone_i[,"nbGroup"][is.na(Metabo_seq_FractionXZone_i[,"nbGroup"])==TRUE] <- "Unassigned"
    # Other
    Metabo_seq_FractionXZone_i$Somme <- rowSums(Metabo_seq_FractionXZone_i[,5:8])
    Metabo_seq_FractionXZone_i[,"Division"][Metabo_seq_FractionXZone_i[,"Somme"]<0.00005*colSums(Metabo_seq_FractionXZone_i[,"Somme"])] <- "Other<0.005%"
    Metabo_seq_FractionXZone_i <- Metabo_seq_FractionXZone_i %>% dplyr::select(-Somme)
    Metabo_seq_FractionXZone_i[,"Supergroup"][Metabo_seq_FractionXZone_i[,"Division"]=="Other<0.005%"] <- "Other<0.005%"
    Metabo_seq_FractionXZone_i <- Metabo_seq_FractionXZone_i %>% group_by(Supergroup,Division,nbGroup,Label) %>% summarise_all(sum)
    #
    Metabo_seq_FractionXZone_i <- reshape2::melt(Metabo_seq_FractionXZone_i,id.vars = c("Label","nbGroup","Supergroup","Division"))
    Metabo_seq_FractionXZone_i <- separate(Metabo_seq_FractionXZone_i,"variable",c("Fraction","Zone"),sep="X")
    Metabo_seq_FractionXZone_i <- Metabo_seq_FractionXZone_i %>% mutate(Division = ifelse(is.na(Division) == TRUE, paste0('Unaffiliated_',Supergroup),Division))
    Metabo_seq_FractionXZone_i[Metabo_seq_FractionXZone_i=="Unaffiliated_NA"] <- "Unaffiliated"
    Metabo_seq_FractionXZone_i[,"Supergroup"][Metabo_seq_FractionXZone_i[,"Supergroup"]=="Unaffiliated_Unaffiliated"] <- "Unaffiliated"
    Metabo_seq_FractionXZone_i <- Metabo_seq_FractionXZone_i %>% filter(Division != "Unaffiliated")
    Metabo_seq_FractionXZone_i[,"nbGroup"] <- sapply(strsplit(as.character(Metabo_seq_FractionXZone_i$nbGroup),": "), `[`, 2)
    Metabo_seq_FractionXZone_i[,"nbGroup"][is.na(Metabo_seq_FractionXZone_i[,"nbGroup"])==TRUE] <- "Unassigned"
    Metabo_seq_FractionXZone_i[,"Supergroup"][Metabo_seq_FractionXZone_i[,"Division"]=="Other<0.005%"] <- "ZZZ"
    # Prop
    for (i in row.names(Metabo_seq_FractionXZone_i)) {
      Fractionx <- Metabo_seq_FractionXZone_i[i,"Fraction"]
      Zonex <- Metabo_seq_FractionXZone_i[i,"Zone"]
      nbGroupx <- Metabo_seq_FractionXZone_i[i,"nbGroup"]
      Metabo_seq_FractionXZone_i[i,"Sum"] <- sum(Metabo_seq_FractionXZone_i %>% filter(Fraction == Fractionx) %>% filter(Zone == Zonex) %>% filter(nbGroup == nbGroupx) %>% dplyr::select(value))
      Metabo_seq_FractionXZone_i[i,"Relative"] <- Metabo_seq_FractionXZone_i[i,"Sum"] *100 / sum(Metabo_seq_FractionXZone_i %>% filter(Fraction == Fractionx) %>% filter(Zone == Zonex) %>% dplyr::select(value))
      Metabo_seq_FractionXZone_i[i,"Proportion"] <- Metabo_seq_FractionXZone_i[i,"value"]*100/Metabo_seq_FractionXZone_i[i,"Sum"]
    }
    Metabo_seq_FractionXZone_i[,"Proportion"][is.na(Metabo_seq_FractionXZone_i[,"Proportion"])==TRUE] <- 0
    Metabo_seq_FractionXZone_i <- Metabo_seq_FractionXZone_i %>% filter(round(Relative,1) != 0)
    # Color
    Metabo_seq_FractionXZone_i <- merge(x=Metabo_seq_FractionXZone_i %>% arrange(Supergroup),y=uniqueTax %>% dplyr::select(-Supergroup),by="Division")
    Metabo_seq_FractionXZone_i <- Metabo_seq_FractionXZone_i %>% arrange(Supergroup)
    color_tax_table_order <- unique(Metabo_seq_FractionXZone_i$Division)
    #plot
    Total_ZoneXFraction_Sequence_fig <- ggplot(Metabo_seq_FractionXZone_i, mapping = aes(y= Proportion, x = Zone, fill = factor(Division, levels = color_tax_table_order), group = factor(Division, levels = color_tax_table_order)), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity",alpha=0.95,color="black") +
      scale_fill_manual(values = unique(Metabo_seq_FractionXZone_i$color)) + facet_grid(Fraction~nbGroup,scales="free") + theme_unique_art() +
      labs(x="Zones",y="Unigenes %",fill="Division") + guides(fill=guide_legend(ncol=1)) +
      geom_label(aes(y = 106,label = paste0(round(Sum,1)," TPM","\n","(",round(Relative,1)," %)")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B")
    print(Total_ZoneXFraction_Sequence_fig)
    ggsave(paste0("Hist_Meta/Grp_ZoneXFraction_",j,"_abundance.svg"), device = "svg", width = 18, height = 10)
    #
    }}
#
  # Proportion --------------------------------------------------------------
  # Color
  paletspe <- c("#b71c1d","#c62828","#d32f2e","#e53935","#f24336","#ef5350",
                "#880e4f","#ad1557","#c2185a",
                "#4a148c","#6a1c9a",
                "#33691d","#558b2f","#689f38","#7cb342","#8bc34a","#9ccc65",
                "#164c40","#21695c","#28796c",
                "#2447a1","#3065c0","#3775d2","#3f88e5","#4496f3","#48a5f5",
                "#f6c02d","#f8d835","#f9eb3b","#faee58",#"#fbf176",
                "#4527a0",
                "#1b237e","#283593","#303e9f","#3949ab",
                #"#e64a18","#f2511e","#f25823","#f37043",
                "#2f838f","#3897a7","#41acc1","#48bcd4","#4cc6da",
                #"#424242","#616161","#757575","#9e9e9e","#bdbdbd",
                "#bdbdbd")
  # Taxa
  uniqueTax <- unique(tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select("Supergroup","Division") %>% arrange(Supergroup))
  uniqueTax <- uniqueTax %>% mutate(Division = ifelse(is.na(Division) == TRUE, paste0('Unaffiliated_',Supergroup),Division))
  uniqueTax[uniqueTax=="Unaffiliated_NA"] <- "Unaffiliated"
  uniqueTax <- uniqueTax %>% mutate(Division = ifelse(grepl(Division,pattern = "_X") == TRUE, paste0('Unaffiliated_',Supergroup),Division)) %>% distinct(Division,Supergroup) %>% arrange(Supergroup,Division)
  
  uniqueTax$color <- unique(paletspe)
  uniqueTax <- rbind(uniqueTax,data.frame(Supergroup="ZZZ",Division="Other<0.005%",color="#f36f02"))
  # launch loop
  Metabo_seq_FractionXZone <- dplyr::left_join(FractionXZone_seq,tableVinput_abundance_normalize_annots_tax_Funct %>% dplyr::select("Row.names","nbGroup","Kegg.Onthology","Supergroup","Division"),by="Row.names")
  Metabo_seq_FractionXZone <- separate_rows(Metabo_seq_FractionXZone,"Kegg.Onthology",sep = ",") # Duplicate row by KO 
  Metabo_seq_FractionXZone <- Metabo_seq_FractionXZone %>% dplyr::select(-"Row.names") %>% group_by(Kegg.Onthology,nbGroup,Supergroup,Division) %>% summarise_all(sum)
  Metabo_seq_FractionXZone <- left_join(Metabo_seq_FractionXZone,ko2Met_Table_Final_i,by=c("Kegg.Onthology"="KO_id"))
  colnames(Metabo_seq_FractionXZone)[colnames(Metabo_seq_FractionXZone) == 'Pathway/Gene'] <- 'Label'
  Metabo_seq_FractionXZone[,"Label"] <- sapply(strsplit(as.character(Metabo_seq_FractionXZone$Label),"_-_"), `[`, 1)
  Metabo_seq_FractionXZone[,"Label"] <- sapply(str_replace_all(as.character(Metabo_seq_FractionXZone$Label),"_"," "), `[`, 1)
  Metabo_seq_FractionXZone <- as.data.frame(Metabo_seq_FractionXZone)
  for (meta in unique(Metabo_seq_FractionXZone$Label)) { if (is.na(meta) == FALSE) { print (meta)
    meta_ <- str_replace_all(meta," ","_")
    # Select meta
    Metabo_seq_FractionXZone_i <- Metabo_seq_FractionXZone %>% filter(Label==meta)
    Metabo_seq_FractionXZone_i <- Metabo_seq_FractionXZone_i %>% dplyr::select(-"Kegg.Onthology") %>% group_by(nbGroup,Label,Supergroup,Division) %>% summarise_all(sum)
    Metabo_seq_FractionXZone_i[,"nbGroup"][is.na(Metabo_seq_FractionXZone_i[,"nbGroup"])==TRUE] <- "Unassigned"
    Metabo_seq_FractionXZone_i <- as.data.frame(Metabo_seq_FractionXZone_i) %>% filter(nbGroup != "Unassigned")
    # Other
    if (nrow(Metabo_seq_FractionXZone_i) > 0) {
    Metabo_seq_FractionXZone_i$Somme <- rowSums(Metabo_seq_FractionXZone_i[,5:8])
    Metabo_seq_FractionXZone_i[,"Division"][Metabo_seq_FractionXZone_i[,"Somme"]<0.00005*sum(Metabo_seq_FractionXZone_i[,"Somme"])] <- "Other<0.005%"
    Metabo_seq_FractionXZone_i <- Metabo_seq_FractionXZone_i %>% dplyr::select(-Somme)
    Metabo_seq_FractionXZone_i[,"Supergroup"][Metabo_seq_FractionXZone_i[,"Division"]=="Other<0.005%"] <- "Other<0.005%"
    Metabo_seq_FractionXZone_i <- Metabo_seq_FractionXZone_i %>% group_by(Supergroup,Division,nbGroup,Label) %>% summarise_all(sum)
    #
    Metabo_seq_FractionXZone_i <- reshape2::melt(Metabo_seq_FractionXZone_i,id.vars = c("Label","nbGroup","Supergroup","Division"))
    Metabo_seq_FractionXZone_i <- separate(Metabo_seq_FractionXZone_i,"variable",c("Fraction","Zone"),sep="X")
    Metabo_seq_FractionXZone_i <- Metabo_seq_FractionXZone_i %>% mutate(Division = ifelse(is.na(Division) == TRUE, paste0('Unaffiliated_',Supergroup),Division))
    Metabo_seq_FractionXZone_i[Metabo_seq_FractionXZone_i=="Unaffiliated_NA"] <- "Unaffiliated"
    Metabo_seq_FractionXZone_i[,"Supergroup"][Metabo_seq_FractionXZone_i[,"Supergroup"]=="Unaffiliated_Unaffiliated"] <- "Unaffiliated"
    Metabo_seq_FractionXZone_i <- Metabo_seq_FractionXZone_i %>% filter(Division != "Unaffiliated")
    Metabo_seq_FractionXZone_i[,"nbGroup"] <- sapply(strsplit(as.character(Metabo_seq_FractionXZone_i$nbGroup),": "), `[`, 2)
    Metabo_seq_FractionXZone_i[,"nbGroup"][is.na(Metabo_seq_FractionXZone_i[,"nbGroup"])==TRUE] <- "Unassigned"
    Metabo_seq_FractionXZone_i[,"Supergroup"][Metabo_seq_FractionXZone_i[,"Division"]=="Other<0.005%"] <- "ZZZ"
    # Prop
    for (i in row.names(Metabo_seq_FractionXZone_i)) { if ( !(Metabo_seq_FractionXZone_i[i,"nbGroup"] %in% c("Unassigned","Lichen","Other multicellular Fungi", "Embryophyceae"))) {
      Fractionx <- Metabo_seq_FractionXZone_i[i,"Fraction"]
      Zonex <- Metabo_seq_FractionXZone_i[i,"Zone"]
      nbGroupx <- Metabo_seq_FractionXZone_i[i,"nbGroup"]
      Metabo_seq_FractionXZone_i[i,"Sum"] <- sum(Metabo_seq_FractionXZone_i %>% filter(Fraction == Fractionx) %>% filter(Zone == Zonex) %>% filter(nbGroup == nbGroupx) %>% dplyr::select(value))
      Metabo_seq_FractionXZone_i[i,"Relative"] <- Metabo_seq_FractionXZone_i[i,"Sum"] *100 / Metabo_seq_FractionXZone_k %>% filter(Label == meta) %>% filter(Fraction == Fractionx) %>% filter(Zone == Zonex) %>% filter(nbGroup == nbGroupx) %>% dplyr::select(Sum)
      Metabo_seq_FractionXZone_i[i,"Proportion"] <- Metabo_seq_FractionXZone_i[i,"value"]*100/Metabo_seq_FractionXZone_i[i,"Sum"]
    }}
    Metabo_seq_FractionXZone_i[,"Proportion"][is.na(Metabo_seq_FractionXZone_i[,"Proportion"])==TRUE] <- 0
    Metabo_seq_FractionXZone_i <- Metabo_seq_FractionXZone_i %>% filter(round(Relative,1) != 0)
    # Color
    Metabo_seq_FractionXZone_i <- merge(x=Metabo_seq_FractionXZone_i %>% arrange(Supergroup),y=uniqueTax %>% dplyr::select(-Supergroup),by="Division")
    Metabo_seq_FractionXZone_i <- Metabo_seq_FractionXZone_i %>% arrange(Supergroup)
    color_tax_table_order <- unique(Metabo_seq_FractionXZone_i$Division)
    #plot
    Total_ZoneXFraction_Sequence_fig <- ggplot(Metabo_seq_FractionXZone_i, mapping = aes(y= Proportion, x = Zone, fill = factor(Division, levels = color_tax_table_order), group = factor(Division, levels = color_tax_table_order)), Rowv = NA, col = colMain, scale = "column") + geom_bar(stat="identity",alpha=0.95,color="black") +
      scale_fill_manual(values = unique(Metabo_seq_FractionXZone_i$color)) + facet_grid(Fraction~nbGroup,scales="free") + theme_unique_art() +
      labs(title= meta,x="Zones",y="Unigenes %",fill="Division") + guides(fill=guide_legend(ncol=1)) +
      geom_label(aes(y = 106,label = paste0(round(Sum,1)," TPM","\n","(",round(Relative,1)," %)")),color = "white",size = 3,show.legend = FALSE, fill = "#3B3B3B")
    print(Total_ZoneXFraction_Sequence_fig)
    ggsave(paste0("Hist_Meta/PRP_ZoneXFraction_",meta_,"_abundance.svg"), device = "svg", width = 18, height = 10)
    #
  }}}
  