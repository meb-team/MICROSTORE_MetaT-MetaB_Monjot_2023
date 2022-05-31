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

ARG=$1
Input=$(echo $ARG | cut -d"_" -f1)
Hash=40000
FWD=XFWDX
REV=XREVX
## Fusionner les R1 et R2

mkdir $Input"_Fuse/"
for R1 in $(ls | grep ".fastq" | grep $Input"_[12]_R1_")
do
    echo $R1
    R2=$(echo $R1 | sed 's/_R1_/_R2_/g' )
    echo $R2
    sample=$(echo $R1 | awk -F"_" '{ print $1 }')
    echo $sample
    cat $R1 $R2 > $Input"_Fuse/"$sample"_unify.fastq"
done

## Sort by primer
cd $Input"_Fuse/"
mkdir $Input"_Temp/"
for sample in $(ls | grep "_unify.fastq")
do
    ID=$(echo $sample | cut -d"_" -f1)
    LineN=$(cat $sample | wc -l)
    testN=1
    while [ $(echo $testN $Hash | awk ' {print $1 * $2} ' ) -lt $LineN ]
    do
    testN=$(expr $testN + 1)
    done
    testN=$(expr $testN - 1)
    echo $testN
    for slice in $(seq 1 1 $testN)
    do
        echo $slice
        head -n$(expr $slice \* $Hash) $sample | tail -n$Hash > $ID"_"$slice"_unify.fastq"
        LineNslice=$(cat $ID"_"$slice"_unify.fastq" | wc -l)
        for i in $(seq 4 4 $LineNslice)
        do
            lecture=$(head -n$i $ID"_"$slice"_unify.fastq" | tail -n4)
            if [ $(echo $lecture | head -n2 | tail -n1 | grep $FWD | wc -l) -eq 1 ]
            then
                if [ $(echo $lecture | head -n1 | grep "/2" | wc -l) -eq 1 ]
                then
                correctR1=$(echo $lecture | sed 's/\/2/\/1/')
                else
                correctR1=$(echo $lecture)
                fi
            echo $correctR1 >> $Input"_Temp/"$ID"_R1.fastq"
            elif [ $(echo $lecture | head -n2 | tail -n1 | grep $REV | wc -l) -eq 1 ]
            then
                if [ $(echo $lecture | head -n1 | grep "/1" | wc -l) -eq 1 ]
                then
                correctR2=$(echo $lecture | sed 's/\/1/\/2/')
                else
                correctR2=$(echo $lecture)
                fi
            echo $correctR2 >> $Input"_Temp/"$ID"_R2.fastq"
            else
            echo $lecture >> $Input"_Temp/"$ID"_NoPrimer.fastq"
            fi
        done
        rm $ID"_"$slice"_unify.fastq"
    done
    for slice in $(expr $testN + 1)
    do
        echo $slice
        rest=$(expr $LineN - $(echo $testN $Hash | awk ' {print $1 * $2} ' ))
        echo $rest
        tail -n$rest $sample > $ID"_"$slice"_unify.fastq"
        LineNslice=$(cat $ID"_"$slice"_unify.fastq" | wc -l)
        for i in $(seq 4 4 $LineNslice)
        do
            lecture=$(head -n$i $ID"_"$slice"_unify.fastq" | tail -n4)
            if [ $(echo $lecture | head -n2 | tail -n1 | grep $FWD | wc -l) -eq 1 ]
            then
                if [ $(echo $lecture | head -n1 | grep "/2" | wc -l) -eq 1 ]
                then
                correctR1=$(echo $lecture | sed 's/\/2/\/1/')
                else
                correctR1=$(echo $lecture)
                fi
            echo $correctR1 >> $Input"_Temp/"$ID"_R1.fastq"
            elif [ $(echo $lecture | head -n2 | tail -n1 | grep $REV | wc -l) -eq 1 ]
            then
                if [ $(echo $lecture | head -n1 | grep "/1" | wc -l) -eq 1 ]
                then
                correctR2=$(echo $lecture | sed 's/\/1/\/2/')
                else
                correctR2=$(echo $lecture)
                fi
            echo $correctR2 >> $Input"_Temp/"$ID"_R2.fastq"
            else
            echo $lecture >> $Input"_Temp/"$ID"_NoPrimer.fastq"
            fi
        done
        rm $ID"_"$slice"_unify.fastq"
    done
done

## change space by tab
cd "../"$Input"_Fuse/"$Input"_Temp/"
for samples in $(ls | grep ".fastq")
do
    ID=$(echo $samples | cut -d"." -f1)
    cat $samples | tr " " "\n" > $ID"_correct.fastq"
    rm $samples
done

## remove unpaired reads
#for R1reads in $(ls | grep "R1_correct.fastq")
#do
#    echo $R1reads
#    R2reads=$(echo $R1reads | sed 's/R1/R2/')
#    echo $R2reads
#    ID=$(echo $R1reads | cut -d"_" -f1)
#    cat $R1reads | grep "^@H7" > $ID"_Access.txt"
#    for identifierR1 in $(cat $ID"_Access.txt")
#    do
#        identifierR2=$(echo $identifierR1 | sed 's/\/1/\/2/')
#        if [ $(cat $R2reads | grep $identifierR2 | wc -l) -eq 1 ]
#        then
#        cat $R1reads | grep -A 3 $identifierR1 >> $ID"_R1_paired_correct.fastq"
#        cat $R2reads | grep -A 3 $identifierR2 >> $ID"_R2_paired_correct.fastq"
#        fi
#    done
#    rm $R1reads
#    rm $R2reads
#done

## remove unpaired reads
for R1reads in $(ls | grep "R1_correct.fastq")
do
    echo $R1reads
    R2reads=$(echo $R1reads | sed 's/R1/R2/')
    echo $R2reads
    ID=$(echo $R1reads | cut -d"_" -f1)
    R1=$(echo $R1reads | cut -d"_" -f2)
    ../../bin/fastq-pair/fastq_pair $R1reads $R2reads
    if [ $(ls ../../ | grep "^Reads_Clean$" | wc -l) -eq 0 ]
    then
        mkdir ../../Reads_Clean
    fi
    for file in $(ls | grep ".fastq.paired.fq")
    do
        if [ $(echo $file | grep "R1" | wc -l) == 1 ]
        then
        mv $file "../../Reads_Clean/"$ID"_R1_correct_paired.fastq"
        fi
        if [ $(echo $file | grep "R2" | wc -l) == 1 ]
        then
        mv $file "../../Reads_Clean/"$ID"_R2_correct_paired.fastq"
        fi
    done
    rm $R1reads
    rm $R2reads
done

## remove temp file
rm -r "../../"$Input"_Fuse"
