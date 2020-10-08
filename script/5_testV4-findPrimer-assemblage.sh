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
#!/bin/bash

## Boucle for
cd ..
cp -r dataPANAM/V4 dataPANAM/V4-testPrimer

for f in $(echo dataPANAM/V4-testPrimer)
do
    for fichier in $(ls $f/*_R1_*.fastq)
    do
        echo $fichier
        reverse=$(echo $fichier | sed 's/R1/R2/g')
        echo $reverse
    
        label=$(echo $fichier | awk -F"/" '{ print $3}' |sed -E 's/(F.*_)S.*.fastq$/\1/')
        echo $label
        
#Assemblage
        vsearch -fastq_mergepairs $fichier -reverse $reverse -fastq_maxmergelen 500 -fastq_minlen 200 -fastq_minovlen 50 -fastq_maxdiffs 0 -fastqout $f/assemble.fastq
#filter
        vsearch -fastq_filter $f/assemble.fastq -fastaout $f/assemble.fasta -relabel $label\_ -fasta_width 0
#chimère
        vsearch -uchime_denovo $f/assemble.fasta -nonchimeras $f/nonchimera-$label.fasta -fasta_width 0
    done
done

