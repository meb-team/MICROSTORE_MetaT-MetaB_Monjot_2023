#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 27/02/2020
#
#!/bin/bash

#Ask V4 or V9
if [ $(echo $1 | grep "." |wc -l ) == 1 ]
then
region=$1
else
echo 'Enter ribosomal region (V4 or V9) : '
read region
fi
#Ask assembler
if [ $(echo $2 | grep "." |wc -l ) == 1 ]
then
assembler=$2
else
echo 'Enter assembler (vsearch, NGmerge, pear, flash2 or fastxtend) : '
read assembler
fi

if [ $region == "V4" ]
then
## Sort primers
    cd ../dataPANAM
    result=$(echo 'V4-testPrimer_'$assembler)
    echo $result
    for f in $(ls $result/*-*.fastq.fasta)
    do
        x=$(grep "^TTGG[CT][AG]AATGCTTTCGC.*TACCGCGGC[GT]GCTG[AG]CAC$" $f | wc -l)
        y=$(grep "^GTG[CT]CAGC[AC]GCCGCGGTA.*GCGAAAGCATT[CT][AG]CCAA$" $f | wc -l)
        w=$(cat $f | grep -v ">" | grep "N" | wc -l)
        z=$(($x+$y))
        m=$(grep "^>" $f | wc -l)
        n=$(($z*100/$m))
        echo $f
        amp=$(echo $f | awk -F"-" '{ print $3 }' | awk -F"_" '{ print $1 }')
        echo $amp
        forward=$(ls $result/$amp'_2_R1_clean.fastq')
        echo $forward
        h=0
        reverse=$(echo $forward | sed 's/R1/R2/g')
        echo $reverse
        for i in $(cat $forward | grep "^@" | grep "1$" | awk -F"/" '{ print $1 }')
        do
            j=$(grep $i $reverse | wc -l)
            h=$(($h+$j))
        done
        merged=$(($m*100/$h))
        echo $amp" Total Sequences : "$h" - Sequences merged : "$m" ("$merged"%) - Sequences with primers : "$z" ("$n"%) and "$w" sequences with Ns" | tee -a $result/V4-primer.out
        
    done
fi
if [ $region == "V9" ]
then
    cd ../dataPANAM
    result=$(echo 'V9-testPrimer_'$assembler)
    echo $result
    for f in $(ls $result/*-*.fastq.fasta)
    do
        x=$(grep "^TTGTACACACCGCCC.*GTAGGTGAACCTGC[AG]GAAGG$" $f | wc -l)
        y=$(grep "^CCTTC[CT]GCAGGTTCACCTAC.*GGGCGGTGTGTACAA$" $f | wc -l)
        w=$(cat $f | grep -v ">" | grep "N" | wc -l)
        z=$(($x+$y))
        m=$(grep "^>" $f | wc -l)
        n=$(($z*100/$m))
        echo $f
        amp=$(echo $f | awk -F"-" '{ print $3 }' | awk -F"_" '{ print $1 }')
        echo $amp
        forward=$(ls $result/$amp'_1_R1_clean.fastq')
        echo $forward
        h=0
        reverse=$(echo $forward | sed 's/R1/R2/g')
        echo $reverse
        for i in $(cat $forward | grep "^@" | grep "1$" | awk -F"/" '{ print $1 }')
        do
            j=$(grep $i $reverse | wc -l)
            h=$(($h+$j))
        done
        merged=$(($m*100/$h))
        echo $amp" Total Sequences : "$h " - Sequences merged : "$m" ("$merged"%) - Sequences with primers : "$z" ("$n"%) and" $w "sequences with Ns" | tee -a $result/V9-primer.out
    done
fi
    

