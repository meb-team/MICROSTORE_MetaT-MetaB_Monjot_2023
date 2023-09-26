#!/bin/bash
#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 15/04/2021
#
#
#
cd script/
PATHCONDA=$(conda info | grep -i 'base environment' | awk -F" " '{print $4}')
source $PATHCONDA'/etc/profile.d/conda.sh'
conda env create -f environment_REnv_Monjot_2023A.yml ; conda activate REnv_Monjot_2023A
wget https://github.com/marbl/Krona/releases/download/v2.8/KronaTools-2.8.tar
tar -xvf KronaTools-2.8.tar
perl KronaTools-2.8/install.pl --prefix ./KronaTools-2.8
rm KronaTools-2.8.tar

Rscript 3_Install_dependencies.R
