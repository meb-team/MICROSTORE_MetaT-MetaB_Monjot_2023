#!/bin/bash
#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 24/05/2022
#
BEFORE=$SECONDS
# Preprocess
echo "Pre-processing start"
cd script/
if [ $(ls ../ | grep "dataDADA2" | wc -l) -eq 0 ]
then
bash 1_Pre-process.sh
fi

# End
ELAPSED=$((($SECONDS-$BEFORE)/60))
echo "Pre-processing stage is completed and takes "$ELAPSED" minutes"
