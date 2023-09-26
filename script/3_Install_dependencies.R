  # !/usr/bin/env Rscript
#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 02/12/2022
#
# Install dependencies in conda environment 
# Set directory, detect R or Rstudio -----------------------------------------------------------
se <- Sys.getenv()
if (length(se[grepl("rstudio",se,ignore.case=TRUE)]) > 0 ) { inputmode <- TRUE }
if (length(se[grepl("rstudio",se,ignore.case=TRUE)]) == 0 ) { inputmode <- FALSE }
#
if (inputmode == TRUE) {
  current <- dirname(rstudioapi::getSourceEditorContext()$path)
  setwd(current)}
#
# Set directory and import package -----------------------------------------------------------
pkgs <- c("GUniFrac","ggplot2","tidyr","dplyr","cowplot","ggrepel","ggsci","scales","varhandle","treemap","VennDiagram","FactoMineR","RColorBrewer","factoextra","reshape2","ggpubr","hrbrthemes","svglite","BiocManager", "stringr","paletteer","gtools","vegan","SciViews","cluster","NbClust","tibble","treemapify","psych","gplots","ggExtra","devtools","data.table","rgl","car","scatterplot3d","R.utils","ggh4x")
#
install.packages(pkgs, 
                 repos = paste0("file:",getwd(),"/REnv_Monjot_2023A_packages"),
                 type = "source")
library("devtools")
library("BiocManager")
## Elementalist
devtools::install_github("teunbrand/elementalist",upgrade = "never")
## DADA2
devtools::install_github("benjjneb/dada2", ref="v1.18",upgrade = "never")
## DESeq2
BiocManager::install("DESeq2", update = FALSE)
## ComplexHeatmap
BiocManager::install("ComplexHeatmap", update = FALSE)
## SARTools
devtools::install_github("PF2-pasteur-fr/SARTools", build_opts="--no-resave-data",upgrade = "never")
#
dependencies <- c(pkgs,"elementalist","dada2","DESeq2","ComplexHeatmap")
check <- as.logical(dependencies %in% rownames(installed.packages()))
ok <- as.logical(rep(TRUE, each = length(dependencies)))
w <- 0
for (checkpackage in check) { 
  if (checkpackage == TRUE) { w <- w+1 }}
if ( w == length(dependencies)) { message("Dependencies installed !") }
if ( w != length(dependencies)) { message("Dependencies not installed !") }
#
