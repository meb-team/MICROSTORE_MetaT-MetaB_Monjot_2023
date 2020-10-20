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
cd rawdata ; cat data-inf.txt | awk -F"\t" '{ print $2"_"$3"_"$5 }' | awk -F"_" '{ print $2"_"$3"_"$NF }' | cut -c1-10 | awk -F"_" '{ print $1"_"$2"_"$3"_" }' > sortV4-V9

# Modify FR sample (V4 -> V9) :
sed -i 's/FR_DJAG_V4/FR_DJAG_V9/g' sortV4-V9
    ## Sur mac : sed -i "" 's/FR_DJAG_V4/FR_DJAG_V9/g' sortV4-V9

# 1_Sort_V4_V9
## Sort and Compress Reads
cd reads
mkdir ../V9/
mkdir ../V4/
for i in $(ls)
do
    f=$(grep ^"$i" ../sortV4-V9 | cut -d"_" -f3)
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

# 3_Sort reads for PANAM analyses
## Create dataPANAM directory and sort reads
cd ..
mkdir dataPANAM
mkdir dataPANAM/V4
mkdir dataPANAM/V9
for V in V4 V9
do
    for sample in $(ls rawdata/$V)
    do
        mv rawdata/$V/$sample/*.fastq dataPANAM/$V/
    done
rm -r -f rawdata/$V
done



# 4_Rename rawdata
## Rename reads
for x in $(echo dataPANAM/V*)
do
    for f in $(ls $x/*.fastq)
    do
        rename=$(echo "$f" | awk -F"_" '{ print $2"_"$3"_R"$4"_"$6 }')
        mv $f $x/$rename
        echo $rename
    done
done

# 5_Install PANAM2
cd dataPANAM ; git clone https://github.com/panammeb/PANAM2.git
cd PANAM2 ; perl setup.pl
chmod 777 bd_ssrna

# Move and compress Reads in PANAM2 directory
cd .. ; mv V4 PANAM2/V4 ; gzip PANAM2/V4/*.fastq ; mv V9 PANAM2/V9 ; gzip PANAM2/V9/*.fastq

echo Pre-process completed !
