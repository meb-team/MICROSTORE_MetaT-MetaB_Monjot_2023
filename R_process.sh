#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 25/11/2020
#
#!/bin/bash

if [ $(echo $1 | grep "." |wc -l ) == 1 ]
then
inifile=$1
else
echo 'Enter initialization file path (.ini) : '
read inifile
fi

## Set variable
INPUT=$(cat $inifile | grep "INPUT" | awk -F"\t" '{ print $2 }')
echo "Input file: "$INPUT
INPUT=$( echo "../"$INPUT)
OUTPUT=$(cat $inifile | grep "OUTPUT" | awk -F"\t" '{ print $2 }')
REGION=$(cat $inifile | grep "REGION" | awk -F"\t" '{ print $2 }')
echo "Region: "$REGION
MODE=$(cat $inifile | grep "MODE" | awk -F"\t" '{ print $2 }')
echo "Taxonomic level: "$MODE
if [ $MODE == "Superphylum" ]
then
DIVISION=$(echo "Eukaryota")
else
DIVISION=$(cat $inifile | grep "DIVISION" | awk -F"\t" '{ print $2 }')
fi
echo "Division: "$DIVISION
FILTER=$(cat $inifile | grep "FILTER" | awk -F"\t" '{ print $2 }')
echo "Bockulich filter: "$FILTER
TAXONOMY=$(cat $inifile | grep "TAXONOMY" | awk -F"\t" '{ print $2 }')
echo "Taxonomy mode: "$TAXONOMY
UNIFY=$(cat $inifile | grep "UNIFY" | awk -F"\t" '{ print $2 }')
echo "Unify duplicat: "$UNIFY
RAREFY=$(cat $inifile | grep "RAREFY" | awk -F"\t" '{ print $2 }')
echo "Rarefy global data: "$RAREFY
RARECURVE=$(cat $inifile | grep "RARECURVE" | awk -F"\t" '{ print $2 }')
echo "Calcul rarecurve: "$RARECURVE
IDENTIFIER=$(cat $inifile | grep "IDENTIFIER" | awk -F"\t" '{ print $2 }')
echo "ID: "$IDENTIFIER

## set auto output
if [ $OUTPUT == "auto" ]
then
OUTPUT=$(echo "Composition-"$REGION"-"$IDENTIFIER"-Filter"$FILTER"-Unify"$UNIFY"-Rarefy"$RAREFY"-"$TAXONOMY"-"$DIVISION)
fi
echo "Output file: "$OUTPUT

# set script directory
cd script/
## Launch Analyse_Duplicat
nohup Rscript 8_Analyse_Duplicat_AFC.R $INPUT $OUTPUT $REGION $FILTER $RAREFY >> Ranalyse.log
## Launch Analyse_composition
if [ $MODE == "Superphylum" ]
then
nohup Rscript 9_Analyse_Composition.R $INPUT $OUTPUT $REGION $FILTER $TAXONOMY $UNIFY $RAREFY $MODE >> Ranalyse.log
else
nohup Rscript 9_Analyse_Composition.R $INPUT $OUTPUT $REGION $FILTER $TAXONOMY $UNIFY $RAREFY $MODE $DIVISION >> Ranalyse.log
fi
## Launch Calcul_Rarecurve
if [ $RARECURVE == "yes" ]
then
nohup Rscript 10_Calcul_Rarecurve.R $INPUT $OUTPUT $REGION $FILTER $UNIFY >> Ranalyse.log
fi
