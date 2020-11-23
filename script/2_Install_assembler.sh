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

# Intsall assembler
## Install vsearch -v 2.15.0
cd ../bin
tar xzf v2.15.0.tar.gz
cd vsearch-2.15.0
./autogen.sh
./configure --prefix=$(pwd)
make
mv bin/vsearch ../

## Install pear -v 0.9.11
cd ..
tar -xvf pear-0.9.11-linux-x86_64.tar
mv pear-0.9.11-linux-x86_64/bin/pear ./

## Install FLASH -v 2.2.00
tar -xvf FLASH2-2.2.00.tar
cd FLASH2-2.2.00
make
mv flash2 ../

## Install NGmerge -v 0.3
cd ..
tar -xvf NGmerge-0.3.tar
cd NGmerge-0.3
make
mv NGmerge ../

## Install fastxtend
cd ..
tar -xvf fastxtend-0.0.13.1.tar
mv fastxtend-0.0.13.1/bin/fastx_mergepairs ./

##
echo 'Assembler Installed !'
