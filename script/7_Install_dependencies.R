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
# Install dependencies in conda environment 
# Set directory and import package -----------------------------------------------------------
system("conda install -y -c bioconda r-factominer ; conda install -y -c conda-forge r-factoextra ; conda install -y -c r r-reshape2 ; conda install -y -c conda-forge r-ggpubr ; conda install -y -c conda-forge r-hrbrthemes ; conda install -y -c conda-forge r-svglite")
system("conda install -y -c conda-forge imagemagick")
#system("conda install -y -c bioconda bioconductor-phyloseq")

setwd("..")
install.packages("parallel", repos="http://cran.rstudio.com/")
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
## PANAM dependencies
#install.packages("vegan", repos="http://cran.rstudio.com/")
#install.packages("picante", repos="http://cran.rstudio.com/")
#install.packages("cluster", repos="http://cran.rstudio.com/")
#install.packages("MASS", repos="http://cran.rstudio.com/")
#install.packages("BiocManager", repos="http://cran.rstudio.com/")
#BiocManager::install("phyloseq")

#dependencies <- c("parallel","GUniFrac","ggplot2","dplyr","cowplot","ggrepel","ggsci","scales","varhandle","treemap","FactoMineR","factoextra","reshape2","ggpubr","hrbrthemes","svglite","vegan","picante","cluster","MASS","phyloseq")
dependencies <- c("parallel","GUniFrac","ggplot2","dplyr","cowplot","ggrepel","ggsci","scales","varhandle","treemap","VennDiagram","FactoMineR","factoextra","reshape2","ggpubr","hrbrthemes","svglite", "stringr")
check <- as.logical(dependencies %in% rownames(installed.packages()))
ok <- as.logical(rep(TRUE, each = length(dependencies)))
w <- 0
for (checkpackage in check) { 
  if (checkpackage == TRUE) { w <- w+1 }}
if ( w == length(dependencies)) { message("Dependencies installed !") }
if ( w != length(dependencies)) { message("Dependencies not installed !") }


