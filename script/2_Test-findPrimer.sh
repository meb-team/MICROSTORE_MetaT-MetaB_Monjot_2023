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

if [ $(echo $3 | grep "." | wc -l ) == 1 ]
then
chimere=$3
else
echo 'Check chimeras ? (yes or no) : '
read chimere
fi

if [ $region == "V4" ]
then
## Boucle for
    cd ..
    if [ $(ls dataPANAM/ | grep "V4-testPrimer" | wc -l) == 1 ]
    then
    rm -r dataPANAM/V4-testPrimer
    fi
    cp -r dataPANAM/PANAM2/V4 dataPANAM/V4-testPrimer
    gunzip dataPANAM/V4-testPrimer/*
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
            if [ $assembler == "vsearch" ]
            then
            bin/vsearch -fastq_mergepairs $fichier -reverse $reverse -fastq_maxmergelen 500 -fastq_minlen 200 -fastq_allowmergestagger -fastq_minovlen 50 -fastq_maxdiffs 0 -fastqout $f/assemble.fastq
            fi
            
            if [ $assembler == "NGmerge" ]
            then
            bin/NGmerge -1 $fichier -2 $reverse -m 50 -p 0 -n 16 -o $f/assemble.fastq
            fi
            
            if [ $assembler == "pear" ]
            then
            bin/pear -f $fichier -r $reverse -p 0.0001 -v 50 -m 500 -n 200 -j 16 -o $f/assemble.fastq
            fi
            
            if [ $assembler == "flash2" ]
            then
            bin/flash2 $fichier $reverse -m 50 -r 250 -f 400 -t 16 -o $f/assemble.fastq
            fi
            
            if [ $assembler == "fastxtend" ]
            then
            bin/fastx_mergepairs -a $fichier -b $reverse -l 50 -m 0 -L 200 -i 100 -o $f/assemble.fastq
            fi
            
    #filter
            bin/vsearch -fastq_filter $f/assemble.fastq -fastaout $f/assemble.fasta -relabel $label\_ -fasta_width 0
    #chimère
            if [ $chimere ==  "yes" ]
            then
            bin/vsearch -uchime_denovo $f/assemble.fasta -nonchimeras $f/nonchimera-$label.fasta -fasta_width 0
            else
            mv $f/assemble.fasta $f/assembly-$label.fasta
            fi
        done
    done
fi
if [ $region == "V9" ]
then
    ## Boucle for
    cd ..
    if [ $(ls dataPANAM/ | grep "V9-testPrimer" | wc -l) == 1 ]
    then
    rm -r dataPANAM/V9-testPrimer
    fi
    cp -r dataPANAM/PANAM2/V9 dataPANAM/V9-testPrimer
    gunzip dataPANAM/V9-testPrimer/*
    for f in $(echo dataPANAM/V9-testPrimer)
    do
        for fichier in $(ls $f/*_R1_*.fastq)
        do
            echo $fichier
            reverse=$(echo $fichier | sed 's/R1/R2/g')
            echo $reverse
        
            label=$(echo $fichier | awk -F"/" '{ print $3}' |sed -E 's/(F.*_)S.*.fastq$/\1/')
            echo $label
            
    #Assemblage
            if [ $assembler == "vsearch" ]
            then
            bin/vsearch -fastq_mergepairs $fichier -reverse $reverse -fastq_maxmergelen 300 -fastq_minlen 100  -fastq_allowmergestagger  -fastq_minovlen 50 -fastq_maxdiffs 0 -fastqout $f/assemble.fastq
            fi
            
            if [ $assembler == "NGmerge" ]
            then
            bin/NGmerge -1 $fichier -2 $reverse -m 50 -p 0 -n 16 -o $f/assemble.fastq
            fi
            
            if [ $assembler == "pear" ]
            then
            bin/pear -f $fichier -r $reverse -p 0.0001 -v 50 -m 300 -n 100 -j 16 -o $f/assemble.fastq
            fi
            
            if [ $assembler == "flash2" ]
            then
            bin/flash2 $fichier $reverse -m 50 -r 150 -f 200 -t 16 -o $f/assemble.fastq
            fi
            
            if [ $assembler == "fastxtend" ]
            then
            bin/fastx_mergepairs -a $fichier -b $reverse -l 50 -m 0 -L 100 -i 100 -o $f/assemble.fastq
            fi
            
    #filter
            bin/vsearch -fastq_filter $f/assemble.fastq -fastaout $f/assemble.fasta -relabel $label\_ -fasta_width 0
    #chimère
            if [ $chimere ==  "yes" ]
            then
            bin/vsearch -uchime_denovo $f/assemble.fasta -nonchimeras $f/nonchimera-$label.fasta -fasta_width 0
            else
            mv $f/assemble.fasta $f/assembly-$label.fasta
            fi
        done
    done
fi
