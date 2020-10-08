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
# Import packages

# Set directory and import package -----------------------------------------------------------
system("conda install -y -c bioconda r-factominer ; conda install -y -c conda-forge r-factoextra ; conda install -y -c r r-reshape2 ; conda install -y -c conda-forge r-ggpubr ; conda install -y -c conda-forge r-hrbrthemes ; conda install -y -c conda-forge r-svglite")

setwd("..")
install.packages("GUniFrac", repos="http://cran.rstudio.com/")
install.packages("ggplot2", repos="http://cran.rstudio.com/")
install.packages("dplyr", repos="http://cran.rstudio.com/")
install.packages("cowplot", repos="http://cran.rstudio.com/")
install.packages("ggrepel", repos="http://cran.rstudio.com/")
install.packages("ggsci", repos="http://cran.rstudio.com/")
install.packages("scales", repos="http://cran.rstudio.com/")
install.packages("varhandle", repos="http://cran.rstudio.com/")
