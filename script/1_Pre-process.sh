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
cd ../rawdata
echo "metadata_metaB table is used for preprocessing"
#
# Detection
echo "Prepare Data 0/2"
if [ $(ls reads/ | grep ".fastq" | wc -l) -eq 0 ]
then
for file in $(ls reads/)
do
mv reads/$file/* reads/
rm -r reads/$file/
done
fi

# Sort Reads
cd reads
mkdir ../V4/
for i in $(ls)
do
    j=$(echo $i | cut -d"_" -f2 | sed 's/OSTA//')
    f=$(grep ^$j ../metadata_metaB | cut -d"_" -f4)
    if [ $f = "V4" ]
    then
        if [ $(ls ../V4/ | grep $j | wc -l) -eq 0 ]
        then
        mkdir "../V4/"$j"/"
        fi
        mv $i "../V4/"$j"/"
    fi
done
echo "Data available 0/2"

## Concate amplicons generated through differents flow-cells
echo "Start Processing of Amplicons from multiple flow-cells 1/2"
cd ..
for h in V4
do
    echo $h
    for f in $(ls $h)
    do
        echo "   "$f

        i=$(ls $h/$f | wc -l)
        for x in $(ls $h/$f)
        do
            echo "      "$x
            if [ $i -gt 2 ]
            then
            n=$(echo $h/$f/$x | awk -F"_" '{ print $2 }')
            o=$(echo $h/$f/$x | awk -F"_" '{ print $3 }')
            g=$(echo $h/$f/$x | awk -F"_" '{ print $4 }')
            p=$(echo $h/$f/$x | awk -F"." '{ print $2 }' | cut -d"_" -f1)
            cellid=$(echo $h/$f/$x | awk -F"." '{ print $1 }' | cut -d"_" -f5)
            cat $h/$f/"CIN_"$n"_"$o"_"$g"_"$cellid"."$p"_clean.fastq.gz" >> $h/$f/"CIN_"$n"_"$o"_"$g"_CONCATE-"$p"_clean.fastq.gz"
            rm $h/$f/"CIN_"$n"_"$o"_"$g"_"$cellid"."$p"_clean.fastq.gz"
            fi
        done
    done
done
echo "Processing of Amplicons from multiple flow-cells completed 1/2"

# Create dataDADA2 directory and sort reads
cd ..
mkdir dataDADA2
mkdir dataDADA2/V4
for V in V4
do
    for sample in $(ls rawdata/$V)
    do
        mv rawdata/$V/$sample/*.fastq.gz dataDADA2/$V/
    done
rm -r -f rawdata/$V
done

# Rename reads
for x in $(echo dataDADA2/V*)
do
    for f in $(ls $x/*.fastq.gz)
    do
        rename=$(echo "$f" | awk -F"_" '{ print $2"_"$3"_R"$4"_"$6 }')
        mv $f $x/$rename
    done
done

# sort R1 and R2
echo "Start Processing of duplicats 2/2"
cd dataDADA2/
for V in V4
do
    mkdir $V'-unified'
    mkdir $V'-unified'/R1
    mkdir $V'-unified'/R2
    cp $V/*R1*.fastq.gz $V'-unified'/R1/
    cp $V/*R2*.fastq.gz $V'-unified'/R2/
done
## detect and unify duplicat
for V in V4-unified
do
    echo $V
    if [ $V = "V4-unified" ]
    then
    ID=2
    else
    ID=1
    fi
    for R in R1 R2
    do
        for file in $(ls $V/$R/ | awk -F"OSTA" '{ print $1"_" }')
        do
            if [ $(cat ../rawdata/metadata_metaB | grep $file | awk -F"_" '{ print $5 }') -eq 2 ]
            then
            duplicat=$(cat ../rawdata/metadata_metaB | grep $file | awk -F"_" '{ print $2"_"$3"_"$4 }')
            duplicatx1=$(cat ../rawdata/metadata_metaB |  grep $duplicat | grep "_1_" | awk -F"_" '{ print $1 }')
            duplicatx2=$(cat ../rawdata/metadata_metaB |  grep $duplicat | grep "_2_" | awk -F"_" '{ print $1 }')
            echo "   "$duplicatx1" combine with "$duplicatx2
            cat $V/$R/$duplicatx2"OSTA_"$ID"_"$R"_clean.fastq.gz" >> $V/$R/$duplicatx1"OSTA_"$ID"_"$R"_clean.fastq.gz"
            rm $V/$R/$duplicatx2"OSTA_"$ID"_"$R"_clean.fastq.gz"
            fi
        done
    done
done

echo "Processing of duplicats completed 2/2"

## bind R1 and R2
for V in V4-unified
do
    for R in R1 R2
    do
        mv $V/$R/* $V/
        rm -rf $V/$R/
    done
done

# Move and compress Reads in DADA2 directory
mkdir reads
for file in V4 V4-unified
do
    mv $file reads/$file
done
echo Pre-process completed !
    

