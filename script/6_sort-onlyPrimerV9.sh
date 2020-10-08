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

## Sort primers
cd ../dataPANAM
for f in $(ls V9-testPrimer/nonchimera-*.fastq.fasta)
do
    x=$(grep "^TTGTACACACCGCCC.*GTAGGTGAACCTGC[AG]GAAGG$" $f | wc -l)
    y=$(grep "^CCTTC[CT]GCAGGTTCACCTAC.*GGGCGGTGTGTACAA$" $f | wc -l)
    w=$(cat $f | grep -v ">" | grep "N" | wc -l)
    z=$(echo $(($x+$y)))
    m=$(grep "^>" $f | wc -l)
    n=$(($z*100/$m))
    echo $f "- Sequences with primers : "$z "-" $n"%" "of total sequences -" $m "- and" $w "sequences with Ns" | cut -d"/" -f2 | sed 's/nonchimera-//g' |sed 's/_1_R1_clean.fastq.fasta//g' | tee -a V9-testPrimer/V9-primer.out
done
    
