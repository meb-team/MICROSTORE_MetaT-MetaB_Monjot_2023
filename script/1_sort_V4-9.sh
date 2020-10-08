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

## Sort & Compress
cd ../rawdata/reads
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
##Decompress
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
