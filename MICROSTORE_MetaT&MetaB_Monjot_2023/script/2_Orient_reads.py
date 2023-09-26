#!/usr/bin/env python3
#
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
# Import
import  os, argparse, sys, re, subprocess
from collections import defaultdict
#
# Argument
arg=dict()
## Input
for i in range(1,len(sys.argv)):
    if sys.argv[i]=="-i":
        Input="".join(["\"",sys.argv[i+1],"\""])
        print("input:",Input)
Input=Input.lstrip("\"").split("_")[0]
## Primers
FWD="XFWDX"
REV="XREVX"
## Create directory
if not os.path.exists(Input+'_Fuse/'):
    os.makedirs(Input+'_Fuse/')
    
for R1 in os.listdir():
    reg=re.search('^'+Input+'_[12]_R1',R1)
    if reg:
        print(R1)
        R2=R1.replace('R1','R2')
        print(R2)
        os.system(" ".join(['cat',R1,R2,'>',Input+'_Fuse/'+Input+'_unify.fastq']))

## Create fastq Dict for FWD and REV
###FWD
ID_seq_F=defaultdict(str)
ID_score_F=defaultdict(str)
###REV
ID_seq_R=defaultdict(str)
ID_score_R=defaultdict(str)
#
T_sens=""
with open(Input+'_Fuse/'+Input+'_unify.fastq') as f1:
    for li in f1:
        li=li.rstrip("\n")
        if li.startswith("@H7:"):
            scoreT=0
            ID=li
        regFWD=re.search('^'+FWD,li) #[A]
        #regFWD=re.search(FWD,li) #replace by [A]
        if regFWD:
            ID=ID.replace("/2","/1")
            ID_seq_F[ID]=li
            T_sens="F"
        regREV=re.search('^'+REV,li) # [B]
        if regREV: # [C]
        #regREV=re.search(REV,li) #replace by [B]
        #if regREV and not regFWD: #replace by [C]
            ID=ID.replace("/1","/2")
            ID_seq_R[ID]=li
            T_sens="R"
        if li.startswith("+"):
            scoreT=1
        if scoreT==1:
            if not li.startswith("+"):
                if T_sens=="F":
                    ID=ID.replace("/2","/1")
                    ID_score_F[ID]=li
                if T_sens=="R":
                    ID=ID.replace("/1","/2")
                    ID_score_R[ID]=li
                    
## Write correct reads
with open(Input+'_Fuse/'+Input+'_R1_correct.fastq',"w") as F:
    for ID in ID_seq_F:
        F.write(ID+"\n"+ID_seq_F[ID]+"\n+\n"+ID_score_F[ID])
        F.write("\n")
with open(Input+'_Fuse/'+Input+'_R2_correct.fastq',"w") as R:
    for ID in ID_seq_R:
        R.write(ID+"\n"+ID_seq_R[ID]+"\n+\n"+ID_score_R[ID])
        R.write("\n")

## remove unpaired reads
for R1 in os.listdir(Input+'_Fuse/'):
    reg=re.search('^'+Input+'_R1_correct',R1)
    if reg:
        print(R1)
        R2=R1.replace('R1','R2')
        print(R2)
        os.system(" ".join(['bin/fastq-pair/fastq_pair',Input+'_Fuse/'+R1,Input+'_Fuse/'+R2]))
    if not os.path.exists('Reads_Clean/'):
        os.makedirs('Reads_Clean/')
    for file in os.listdir(Input+'_Fuse/'):
        Region=file.split("_")[1]
        if file.endswith('.fastq.paired.fq'):
            os.system(" ".join(['mv',Input+'_Fuse/'+file,'Reads_Clean/'+Input+'_'+Region+'_correct_paired.fastq']))
os.system('rm -rf '+Input+'_Fuse/')
