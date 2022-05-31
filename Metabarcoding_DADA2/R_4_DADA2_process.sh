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
FWD=$(cat $inifile | grep "FWD" | awk -F"\t" '{ print $2 }')
echo "Forward primer : "$FWD
REV=$(cat $inifile | grep "REV" | awk -F"\t" '{ print $2 }')
echo "Reverse primer : "$REV
CUTADAPT=$(which cutadapt)
echo "CUTADAPT locality : "$CUTADAPT


# set script directory
cd script/
## Launch DADA2 script
if [ $(ls ../dataDADA2/rawTable/ | grep $RESULT | wc -l) != 1 ]
then
    echo "Start DADA2"
    nohup Rscript 4_Dada2.R $INPUT $DATABASE $RESULT $MINLEN $MAXLEN $MAXN $MINOVERLAP $MAXMISMATCH $NTHREADS $FWD $REV $CUTADAPT >> Ranalyse.log
    echo "DADA2 processing completed"
fi


ELAPSED=$((($SECONDS-$BEFORE)/60))
echo "R process finished and takes "$ELAPSED" minutes !"
