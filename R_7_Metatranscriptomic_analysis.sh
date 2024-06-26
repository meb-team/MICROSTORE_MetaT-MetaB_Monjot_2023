#!/bin/bash
#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 11/05/2022
#
BEFORE=$SECONDS
# R CStack
ulimit -s unlimited
# Activate conda environment
PATHCONDA=$(conda info | grep -i 'base environment' | awk -F" " '{print $4}')
source $PATHCONDA'/etc/profile.d/conda.sh'
conda activate REnv_Monjot_2023A
if [ $(echo $1 | grep "." |wc -l ) == 1 ]
then
inifile=$1
else
echo 'Enter initialization file path (.ini) : '
read inifile
fi

## Set variable
### Composition variables
INPUT=$(cat $inifile | grep "INPUT" | awk -F"\t" '{ print $2 }')
echo "Unigene table: "$INPUT
OUTPUT=$(cat $inifile | grep "OUTPUT" | awk -F"\t" '{ print $2 }')
echo "Output: "$OUTPUT
TAX=$(cat $inifile | grep "TAX" | awk -F"\t" '{ print $2 }')
echo "Unigene Taxonomy: "$TAX
DATABASE=$(cat $inifile | grep "DATABASE" | awk -F"\t" '{ print $2 }')
echo "Database path used in the trait study: "$DATABASE

# set script directory
cd script/

## Launch Duplicat analysis script
echo "Start R processing: Duplicat analysis 0/2"
nohup Rscript 9_Duplicat_Metatrans_analize.R $INPUT $OUTPUT  >> Ranalyse.log
echo "Duplicat analysis completed 1/2"

## Launch KO metatranscriptomics analysis
echo "Start R processing: KO analysis 1/2"
nohup Rscript 10_KO_Metatrans_DESeq2.R $INPUT $OUTPUT $TAX $DATABASE >> Ranalyse.log
echo "KO analysis completed 2/2"

ELAPSED=$((($SECONDS-$BEFORE)/60))
echo "R process finished and takes "$ELAPSED" minutes !"
