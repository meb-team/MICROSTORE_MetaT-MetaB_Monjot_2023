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
#!/bin/bash
BEFORE=$SECONDS
PATHCONDA=$(conda info | grep -i 'base environment' | awk -F" " '{print $4}')
source $PATHCONDA'/etc/profile.d/conda.sh'
conda activate REnv4

if [ $(echo $1 | grep "." |wc -l ) == 1 ]
then
inifile=$1
else
echo 'Enter initialization file path (.ini) : '
read inifile
fi

## Set variable
### DADA2 variables
INPUT=$(cat $inifile | grep "INPUT" | awk -F"\t" '{ print $2 }')
echo "Input file: "$INPUT
DATABASE=$(cat $inifile | grep "DATABASE" | awk -F"\t" '{ print $2 }')
echo "Database: "$DATABASE
RESULT=$(cat $inifile | grep "RESULT" | awk -F"\t" '{ print $2 }')
echo "Result: "$RESULT
MINLEN=$(cat $inifile | grep "MINLEN" | awk -F"\t" '{ print $2 }')
echo "minLen: "$MINLEN
MAXLEN=$(cat $inifile | grep "MAXLEN" | awk -F"\t" '{ print $2 }')
echo "maxLen: "$MAXLEN
MAXN=$(cat $inifile | grep "MAXN" | awk -F"\t" '{ print $2 }')
echo "maxN : "$MAXN
MINOVERLAP=$(cat $inifile | grep "MINOVERLAP" | awk -F"\t" '{ print $2 }')
echo "minOverlap : "$MINOVERLAP
MAXMISMATCH=$(cat $inifile | grep "MAXMISMATCH" | awk -F"\t" '{ print $2 }')
echo "maxMismatch : "$MAXMISMATCH
NTHREADS=$(cat $inifile | grep "NTHREADS" | awk -F"\t" '{ print $2 }')
echo "Number of threads : "$NTHREADS

### Composition variables
REGION=$(cat $inifile | grep "REGION" | awk -F"\t" '{ print $2 }')
echo "Region: "$REGION
MODE=$(cat $inifile | grep "MODE" | awk -F"\t" '{ print $2 }')
echo "Taxonomic level: "$MODE
if [ $MODE == "Superphylum" ]
then
DIVISION=$(echo "Eukaryota")
DIVISION=Eukaryota
else
DIVISION=$(cat $inifile | grep "DIVISION" | awk -F"\t" '{ print $2 }')
fi
echo "Division: "$DIVISION
FILTER=$(cat $inifile | grep "FILTER" | awk -F"\t" '{ print $2 }')
echo "Bokulich filter: "$FILTER
UNIFY=$(cat $inifile | grep "UNIFY" | awk -F"\t" '{ print $2 }')
echo "Unify duplicat: "$UNIFY
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

## Launch DADA2 script

if [ $(ls ../dataDADA2/result/ | grep $RESULT | wc -l) != 1 ]
then
echo "Start DADA2"
nohup Rscript 8_Dada2.R $INPUT $DATABASE $RESULT $MINLEN $MAXLEN $MAXN $MINOVERLAP $MAXMISMATCH $NTHREADS >> Ranalyse.log
fi

## Launch Composition script
if [ $COMPOSITION == "yes" ]
then
echo "Start Composition analysis"
nohup Rscript 9_Analyse_Composition.R  $RESULT $REGION $FILTER $MODE $DIVISION $RAREFY $UNIFY >> Ranalyse.log
fi

## Launch Calcul_Rarecurve
if [ $RARECURVE == "yes" ]
then
nohup Rscript 10_Calcul_Rarecurve.R $INPUT $OUTPUT $REGION $FILTER $UNIFY >> Ranalyse.log
fi

ELAPSED=$((($SECONDS-$BEFORE)/60))
echo "R process finished and takes "$ELAPSED" minutes !"
