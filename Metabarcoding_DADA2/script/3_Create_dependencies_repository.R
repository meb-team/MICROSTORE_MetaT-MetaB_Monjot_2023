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
setwd("..")
#
getOption("Ncpus", 1L)
halfcore <- parallel::detectCores()/2
options(Ncpus = halfcore)
#
install.packages("gert", repos="http://cran.rstudio.com/")
install.packages("FactoMineR", repos="http://cran.rstudio.com/")
install.packages("factoextra", repos="http://cran.rstudio.com/")
install.packages("reshape2", repos="http://cran.rstudio.com/")
install.packages("ggpubr", repos="http://cran.rstudio.com/")
install.packages("hrbrthemes", repos="http://cran.rstudio.com/")
install.packages("svglite", repos="http://cran.rstudio.com/")
install.packages("stringr", repos="http://cran.rstudio.com/")
install.packages("GUniFrac", repos="http://cran.rstudio.com/")
install.packages("ggplot2", repos="http://cran.rstudio.com/")
install.packages("dplyr", repos="http://cran.rstudio.com/")
install.packages("cowplot", repos="http://cran.rstudio.com/")
install.packages("ggrepel", repos="http://cran.rstudio.com/")
install.packages("ggsci", repos="http://cran.rstudio.com/")
install.packages("scales", repos="http://cran.rstudio.com/")
install.packages("varhandle", repos="http://cran.rstudio.com/")
install.packages("treemap", repos="http://cran.rstudio.com/")
install.packages("VennDiagram", repos="http://cran.rstudio.com/")
install.packages("paletteer", repos="http://cran.rstudio.com/")
install.packages("devtools", repos="http://cran.rstudio.com/")
install.packages("gtools", repos="http://cran.rstudio.com/")
install.packages("SciViews", repos="http://cran.rstudio.com/")
install.packages("cluster", repos="http://cran.rstudio.com/")
install.packages("NbClust", repos="http://cran.rstudio.com/")
install.packages("tibble", repos="http://cran.rstudio.com/")
install.packages("treemapify", repos="http://cran.rstudio.com/")
install.packages("psych", repos="http://cran.rstudio.com/")
install.packages("gplots", repos="http://cran.rstudio.com/")
install.packages("ggExtra", repos="http://cran.rstudio.com/")
install.packages("rgl", repos="http://cran.rstudio.com/")
install.packages("data.table", repos="http://cran.rstudio.com/")
install.packages("car", repos="http://cran.rstudio.com/")
install.packages("scatterplot3d", repos="http://cran.rstudio.com/")



devtools::install_github("teunbrand/elementalist")
## DADA2
devtools::install_github("benjjneb/dada2", ref="v1.18")
#install.packages("BiocManager", repos="http://cran.rstudio.com/")
#BiocManager::install(version = '3.12', ask = FALSE)
#BiocManager::install("dada2", version = "3.12", ask = FALSE)

dependencies <- c("GUniFrac","ggplot2","dplyr","cowplot","ggrepel","ggsci","scales","varhandle","treemap","VennDiagram","FactoMineR","factoextra","reshape2","ggpubr","hrbrthemes","svglite", "stringr","dada2","paletteer","elementalist","gtools","vegan","SciViews","cluster","NbClust","tibble","treemapify","psych","gplots","ggExtra","devtools","data.table","rgl","car","scatterplot3d")

check <- as.logical(dependencies %in% rownames(installed.packages()))
ok <- as.logical(rep(TRUE, each = length(dependencies)))
w <- 0
for (checkpackage in check) { 
  if (checkpackage == TRUE) { w <- w+1 }}
if ( w == length(dependencies)) { message("Dependencies installed !") }
if ( w != length(dependencies)) { message("Dependencies not installed !") }


install.packages("miniCRAN", repos="http://cran.rstudio.com/", quiet=TRUE)
library("miniCRAN")
pkgList <- rownames(installed.packages())
dir.create(pth <- file.path("REnv4_Rpackages"))
makeRepo(pkgList, path = pth,repos="http://cran.rstudio.com/", type = "source")

#
