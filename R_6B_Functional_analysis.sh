#!/bin/bash
#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 19/04/2021
#
BEFORE=$SECONDS
PATHCONDA=$(conda info | grep -i 'base environment' | awk -F" " '{print $4}')
source $PATHCONDA'/etc/profile.d/conda.sh'
conda activate REnv_Monjot_2023A
## INPUT
if [ $(echo $1 | grep "." |wc -l ) == 1 ]
then
INPUT=$1
else
echo 'Enter Trait table path : '
read INPUT
fi
## OUTPUT
if [ $(echo $2 | grep "." |wc -l ) == 1 ]
then
OUTPUT=$2
else
echo 'Enter Result path : '
read OUTPUT
fi
## REGION
if [ $(echo $3 | grep "." |wc -l ) == 1 ]
then
REGION=$3
else
echo 'Enter Region (V4 or V9) : '
read REGION
fi

# set script directory
cd script/

## Launch Composition script
echo "Start Functionnal analysis"
nohup Rscript 8B_Analyse_Function_Table.R $INPUT $OUTPUT $REGION >> Ranalyse.log
echo "Functionnal analysis completed"

ELAPSED=$((($SECONDS-$BEFORE)/60))
echo "R process finished and takes "$ELAPSED" minutes !"
