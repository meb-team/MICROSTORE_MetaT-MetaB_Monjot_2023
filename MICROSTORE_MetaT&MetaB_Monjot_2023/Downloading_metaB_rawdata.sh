#!/bin/bash
#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 22/09/2023
#
BEFORE=$SECONDS
# Downloading
echo "Downloading start"
cd script/
if [ $(ls ../../ | grep "^reads$" | wc -l) -eq 0 ]
then
bash 0_download_metaB_data.sh
fi

# End
ELAPSED=$((($SECONDS-$BEFORE)/60))
echo "Downloading stage is completed and takes "$ELAPSED" minutes"
