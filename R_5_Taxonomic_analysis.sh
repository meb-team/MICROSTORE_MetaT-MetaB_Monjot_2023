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

if [ $(echo $1 | grep "." |wc -l ) == 1 ]
then
inifile=$1
else
echo 'Enter initialization file path (.ini) : '
read inifile
fi

## Set variable
### Composition variables
RESULT=$(cat $inifile | grep "RESULT" | awk -F"\t" '{ print $2 }')
echo "DADA2 result directory: "$RESULT
OUTCOMP=$(cat $inifile | grep "OUTCOMP" | awk -F"\t" '{ print $2 }')
echo "Output: "$OUTCOMP
REGION=$(cat $inifile | grep "REGION" | awk -F"\t" '{ print $2 }')
echo "Region: "$REGION
FILTER=$(cat $inifile | grep "FILTER" | awk -F"\t" '{ print $2 }')
echo "Bokulich filter: "$FILTER
RAREFY=$(cat $inifile | grep "RAREFY" | awk -F"\t" '{ print $2 }')
echo "Rarefy global data: "$RAREFY
COMPOSITION=$(cat $inifile | grep "COMPOSITION" | awk -F"\t" '{ print $2 }')
echo "Process composition: "$COMPOSITION
RARECURVE=$(cat $inifile | grep "RARECURVE" | awk -F"\t" '{ print $2 }')
echo "Calcul rarecurve: "$RARECURVE
IDENTIFIER=$(cat $inifile | grep "IDENTIFIER" | awk -F"\t" '{ print $2 }')
echo "ID: "$IDENTIFIER

# set script directory
cd script/

## Launch Composition script
if [ $COMPOSITION == "yes" ]
then
    echo "Start Composition analysis"
    nohup Rscript 5_Analyse_Composition_ASV_DADA2.R $RESULT $OUTCOMP $REGION $FILTER $RAREFY >> Ranalyse.log
    echo "Composition analysis completed"
fi

## Launch taxonomic table generation
echo "Start taxonomic table generation"
nohup Rscript 6_Preview_Tax_ref.R $OUTCOMP >> Ranalyse.log
echo "Taxonomic table generation completed"

## Launch Calcul_Rarecurve
echo "Start rarecurve generation"
if [ $RARECURVE == "yes" ]
then
    if [ $(ls ../dataDADA2/result/$OUTCOMP/ | grep "Rarecurve" | wc -l) != 1 ]
    then
        nohup Rscript 7_Calcul_Rarecurve.R $OUTCOMP $REGION >> Ranalyse.log
    fi
    echo "Rarecurve generation completed"
fi

ELAPSED=$((($SECONDS-$BEFORE)/60))
echo "R process finished and takes "$ELAPSED" minutes !"
