#!/bin/bash
#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 24/05/2022
#
# Install fastq_pair
mkdir bin
cd bin
git clone https://github.com/linsalrob/fastq-pair.git
cd fastq-pair/
cmake -DCMAKE_INSTALL_PREFIX=$(pwd) .
make
make install
cd ../..
# Launch parallelisation
NTHREADS=$1
parallel -j $NTHREADS -k 'gunzip {}' ::: *.fastq.gz
parallel -j $NTHREADS -k 'python3 2_Orient_reads.py -i {}' ::: *_R1_clean.fastq
#parallel -j $NTHREADS -k 'bash 2_Orient_reads.sh -i {}' ::: *_R1_clean.fastq #AM300522
