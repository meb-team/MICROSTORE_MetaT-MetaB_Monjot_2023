#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 21/10/2020
#
#!/bin/bash

echo Enter result directory name : 
read result
cd ..
if [ $(ls | grep "quality_test" | wc -l) == 1 ]
then
rm -r quality_test/
fi
mkdir quality_test
touch quality_test/seqAll.fastq
for file in $(ls dataPANAM/PANAM2/$result/quality_output/ | grep -v ".fasta")
do
    echo $file
    cat dataPANAM/PANAM2/$result/quality_output/$file/quality_output/pairs/tmp.assembled.fastq >> quality_test/seqAll.fastq
done

#Quality test
mkdir quality_test/result/
fastqc quality_test/seqAll.fastq -o quality_test/result/

if [ $(ls quality_test/result/ | wc -l) != 0 ]
then
echo Fastqc quality test completed ! 
else
echo Fastqc quality test aborted ! 
fi
