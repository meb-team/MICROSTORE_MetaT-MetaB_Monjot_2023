#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 20/10/2020
#
#!/bin/bash

cd ..
# Format table "data-inf.txt" :
cd rawdata ; cat data-inf.txt | awk -F"\t" '{ print $2"_"$3"_1_"$5 }' | awk -F"_" '{ print $2"_"$3"_"$4"_"$5"_"$NF }' | sed 's/_01_/_1_/g' | sed 's/_02_/_2_/g' | cut -c1-15 | awk -F"_" '{ print $1"_"$2"_"$3"_"$5"_"$4"_" }' > sortV4-V9

# Modify FR sample (V4 -> V9) :
sed -i 's/FR_DJAG_09_V4_1_/FR_DJAG_09_V9_1_/g' sortV4-V9
    ## Sur mac : sed -i "" 's/FR_DJAG_09_1_V4_/FR_DJAG_09_1_V9_/g' sortV4-V9

# 1_Sort_V4_V9
## Sort and Compress Reads
cd reads
mkdir ../V9/
mkdir ../V4/
for i in $(ls)
do
    f=$(grep ^"$i" ../sortV4-V9 | cut -d"_" -f4)
    tar -cvf $i.tar $i/
    if [ $f == "V9" ]
    then
    mv $i.tar ../V9/
    else
    mv $i.tar ../V4/
    fi
done
## Decompress
cd ../V9/
for g in $(ls)
do
    tar -xvf ../V9/$g ; rm $g
done
cd ../V4/
for h in $(ls)
do
    tar -xvf ../V4/$h ; rm $h
done

# 2_Concate amplicons generated through differents flow-cells
## Decompress and concate
cd ..
for h in V4 V9
do
    echo $h
    for f in $(ls $h)
    do
        echo $f
        gunzip $h/$f/*.fastq.gz
        i=$(ls $h/$f | wc -l)
        echo $i
        for x in $(ls $h/$f)
        do
            echo $x
            if [ $i \> 2 ]
            then
            n=$(echo $h/$f/$x | awk -F"_" '{ print $2 }')
            o=$(echo $h/$f/$x | awk -F"_" '{ print $3 }')
            g=$(echo $h/$f/$x | awk -F"_" '{ print $4 }')
            p=$(echo $h/$f/$x | awk -F"." '{ print $2 }' | cut -d"_" -f1)
            cellid=$(echo $h/$f/$x | awk -F"." '{ print $1 }' | cut -d"_" -f5)
            cat $h/$f/"CIN_"$n"_"$o"_"$g"_"$cellid"."$p"_clean.fastq" >> $h/$f/"CIN_"$n"_"$o"_"$g"_CONCATE-"$p"_clean.fastq"
            rm $h/$f/"CIN_"$n"_"$o"_"$g"_"$cellid"."$p"_clean.fastq"
            fi
        done
    done
done

# 3_Sort reads for DADA2 analyses
## Create dataDADA2 directory and sort reads
cd ..
mkdir dataDADA2
mkdir dataDADA2/V4
mkdir dataDADA2/V9
for V in V4 V9
do
    for sample in $(ls rawdata/$V)
    do
        mv rawdata/$V/$sample/*.fastq dataDADA2/$V/
    done
rm -r -f rawdata/$V
done

# 4_Rename rawdata
## Rename reads
for x in $(echo dataDADA2/V*)
do
    for f in $(ls $x/*.fastq)
    do
        rename=$(echo "$f" | awk -F"_" '{ print $2"_"$3"_R"$4"_"$6 }')
        mv $f $x/$rename
        echo $rename
    done
done

# 5_ Unify
## sort R1 and R2
cd dataDADA2/
for V in V4 V9
do
    mkdir $V'-unified'
    mkdir $V'-unified'/R1
    mkdir $V'-unified'/R2
    cp $V/*R1*.fastq $V'-unified'/R1/
    cp $V/*R2*.fastq $V'-unified'/R2/
done
## detect and unify duplicat
for V in V4-unified V9-unified
do
    if [ $V == "V4-unified" ]
    then
    ID=2
    else
    ID=1
    fi
    for R in R1 R2
    do
        for file in $(ls $V/$R/ | awk -F"OSTA" '{ print $1"_" }')
        do
            if [ $(cat ../rawdata/sortV4-V9 | grep $file | awk -F"_" '{ print $5 }') == 2 ]
            then
            duplicat=$(cat ../rawdata/sortV4-V9 | grep $file | awk -F"_" '{ print $2"_"$3"_"$4 }')
            duplicatx1=$(cat ../rawdata/sortV4-V9 |  grep $duplicat | grep "_1_" | awk -F"_" '{ print $1 }')
            duplicatx2=$(cat ../rawdata/sortV4-V9 |  grep $duplicat | grep "_2_" | awk -F"_" '{ print $1 }')
            #duplicatx1=$(echo $duplicatx | awk -F" " '{ print $1 }')
            #duplicatx2=$(echo $duplicatx | awk -F" " '{ print $2 }')
            echo $duplicatx1
            echo $duplicatx2
            #echo $duplicatx
            cat $V/$R/$duplicatx2"OSTA_"$ID"_"$R"_clean.fastq" >> $V/$R/$duplicatx1"OSTA_"$ID"_"$R"_clean.fastq"
            rm $V/$R/$duplicatx2"OSTA_"$ID"_"$R"_clean.fastq"
            fi
        done
    done
done
## bind R1 and R2
for V in V4-unified V9-unified
do
    for R in R1 R2
    do
        mv $V/$R/* $V/
        rm -rf $V/$R/
    done
done

# Move and compress Reads in DADA2 directory
mkdir reads
for file in V4 V9 V4-unified V9-unified
do
    echo $file
    mv $file reads/$file
done

echo Pre-process completed !
    
    
