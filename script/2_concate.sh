#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 24/01/2020
#
#!/bin/sh

## Decompress and Concate
cd ../rawdata
for h in $(echo V*)
do
    echo $h
    for f in $(ls $h)
    do
        gunzip $h/$f/*.fastq.gz
        i=$(ls $h/$f | wc -l)
        if [ $i \> 2 ]
        then
            
            echo $f
            echo $i
            for x in $(ls $h/$f)
            do
                n=$(echo $h/$f/$x | awk -F"_" '{ print $2 }')
                echo $n
                o=$(echo $h/$f/$x | awk -F"_" '{ print $3 }')
                echo $o
                g=$(echo $h/$f/$x | awk -F"_" '{ print $4 }')
                echo $g
                p=$(echo $h/$f/$x | awk -F"." '{ print $2 }' | cut -d"_" -f1)
                echo $p
                cat $h/$f/CIN_"$n"_"$o"_1_*.*_clean.fastq > $h/$f/CIN_"$n"_"$o"_1_CONCATE-"$p"_clean.fastq
                echo "$h/$f/CIN_"$n"_"$o"_1_CONCATE-"$p"_clean.fastq"
                cat $h/$f/CIN_"$n"_"$o"_2_*.*_clean.fastq > $h/$f/CIN_"$n"_"$o"_2_CONCATE-"$p"_clean.fastq
                echo "$h/$f/CIN_"$n"_"$o"_2_CONCATE-"$p"_clean.fastq"
            done
            for x in $(ls $h/$f)
            do
                rm $h/$f/CIN_*_1_*.*_clean.fastq
                rm $h/$f/CIN_*_2_*.*_clean.fastq
            done
        fi
    done
done
