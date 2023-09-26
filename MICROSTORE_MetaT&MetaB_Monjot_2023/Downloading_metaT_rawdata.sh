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
cd rawdata/
wget "https://zenodo.org/record/8376851/files/main_table.mapping.unique.raw.noHuman.noConta.noMetazoa.annot.tsv"
wget "https://zenodo.org/record/8376851/files/table_taxonomy.perUnigene.allUnigenes.tsv"
# End
ELAPSED=$((($SECONDS-$BEFORE)/60))
echo "Downloading stage is completed and takes "$ELAPSED" minutes"
