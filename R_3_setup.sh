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
BEFORE=$SECONDS
cd script/
PATHCONDA=$(conda info | grep -i 'base environment' | awk -F" " '{print $4}')
source $PATHCONDA'/etc/profile.d/conda.sh'
conda activate REnv_Monjot_2023A
Rscript 3_Install_dependencies.R
#
ELAPSED=$((($SECONDS-$BEFORE)/60))
echo "R setup finished and takes "$ELAPSED" minutes !"
