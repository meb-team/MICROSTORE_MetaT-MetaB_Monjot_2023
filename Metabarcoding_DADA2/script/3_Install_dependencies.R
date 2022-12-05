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
pkgs <- c("GUniFrac","ggplot2","dplyr","cowplot","ggrepel","ggsci","scales","varhandle","treemap","VennDiagram","FactoMineR","factoextra","reshape2","ggpubr","hrbrthemes","svglite", "stringr","paletteer","gtools","vegan","SciViews","cluster","NbClust","tibble","treemapify","psych","gplots","ggExtra","devtools","data.table","rgl","car","scatterplot3d")
#
getOption("Ncpus", 1L)
halfcore <- parallel::detectCores()/2
options(Ncpus = halfcore)
install.packages(pkgs, 
                 repos = paste0("file:",getwd(),"/REnv4_Rpackages"),
                 type = "source")
library("devtools")
## Elementalist
devtools::install_github("teunbrand/elementalist")
## DADA2
devtools::install_github("benjjneb/dada2", ref="v1.18")
#
dependencies <- c(pkgs,"elementalist","dada2")
check <- as.logical(dependencies %in% rownames(installed.packages()))
ok <- as.logical(rep(TRUE, each = length(dependencies)))
w <- 0
for (checkpackage in check) { 
  if (checkpackage == TRUE) { w <- w+1 }}
if ( w == length(dependencies)) { message("Dependencies installed !") }
if ( w != length(dependencies)) { message("Dependencies not installed !") }
#