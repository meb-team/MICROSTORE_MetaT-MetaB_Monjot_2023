#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 26/02/2020
#
#!/bin/sh

## Create dataPANAM directory and sort reads
cd ..
mkdir dataPANAM
mkdir dataPANAM/V4
mkdir dataPANAM/V9
for V in V4 V9
do
    for sample in $(ls rawdata/$V)
    do
        cp rawdata/$V/$sample/*.fastq dataPANAM/$V
    done
done
