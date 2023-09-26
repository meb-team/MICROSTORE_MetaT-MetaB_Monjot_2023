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
# SET variables
## INPUT
if [ $(echo $1 | grep "." |wc -l ) -eq 1 ]
then
INPUT=$1
else
echo 'Enter INPUT file ("V4" or "V4-unified"): '
read INPUT
fi
## Number of threads to process data
if [ $(echo $2 | grep "." | wc -l ) -eq 1 ]
then
NTHREADS=$2
else
echo 'Enter the Number of threads to process data: '
read NTHREADS
fi
## Forward primer
if [ $(echo $3 | grep "." | wc -l ) -eq 1 ]
then
FWD=$3
else
echo 'Enter forward primer (for V4 : "GTG[CT]CAGC[AC]GCCGCGGTA"): '
read FWD
fi
## Reverse primer
if [ $(echo $4 | grep "." | wc -l ) -eq 1 ]
then
REV=$4
else
echo 'Enter reverse primer (for V4 : "TTGG[CT][AG]AATGCTTTCGC"): '
read REV
fi

# Orientation of the reads for ASV processing
cd script/
echo "Orient_reads script on $INPUT start"
echo "with $NTHREADS threads"
echo "Forward primer : $FWD"
echo "Reverse primer : $REV"


#Start Orient
cp 2_Orient_reads_parallel.sh ../dataDADA2/reads/$INPUT
cp 2_Orient_reads.py ../dataDADA2/reads/$INPUT
# cp 2_Orient_reads.sh../dataDADA2/reads/$INPUT #AM300522
cd ../dataDADA2/reads/$INPUT
cat 2_Orient_reads.py | sed "s/XFWDX/$FWD/g" | sed "s/XREVX/$REV/g" > 2_Orient_readsx.py
# cat 2_Orient_reads.sh | sed "s/XFWDX/$FWD/g" | sed "s/XREVX/$REV/g" > 2_Orient_readsx.sh #AM300522
mv 2_Orient_readsx.py 2_Orient_reads.py
#mv 2_Orient_readsx.sh 2_Orient_reads.sh #AM300522
bash 2_Orient_reads_parallel.sh $NTHREADS

for file in $(ls | grep -v "^Reads_Clean$")
do
rm -rf $file
done
for fastq in $(ls Reads_Clean)
do
mv Reads_Clean/$fastq ../$INPUT/$fastq
done
rm -r Reads_Clean
mv ../$INPUT"/" ../$INPUT"-correct-paired/"
#

ELAPSED=$((($SECONDS-$BEFORE)/60))
echo "Orient_reads stage is completed and takes $ELAPSED minutes"

