#!/bin/bash
#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 22/09/2023
#
mkdir ../rawdata/reads
cd ../rawdata/reads
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259103/CIN_FFOSTA_2_1_HV2Y7BCX2.12BA354_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259103/CIN_FFOSTA_2_2_HV2Y7BCX2.12BA354_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259065/CIN_DTOSTA_2_1_HV2Y7BCX2.12BA373_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259065/CIN_DTOSTA_2_2_HV2Y7BCX2.12BA373_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259075/CIN_EDOSTA_2_1_HV2Y7BCX2.12BA303_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259094/CIN_EWOSTA_2_1_HV2Y7BCX2.12BA341_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259107/CIN_FJOSTA_2_1_HV2Y7BCX2.12BA307_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259107/CIN_FJOSTA_2_2_HV2Y7BCX2.12BA307_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259090/CIN_ESOSTA_2_1_HV2Y7BCX2.12BA293_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259090/CIN_ESOSTA_2_2_HV2Y7BCX2.12BA293_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259059/CIN_DNOSTA_2_1_HV2Y7BCX2.12BA301_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259094/CIN_EWOSTA_2_2_HV2Y7BCX2.12BA341_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259108/CIN_FKOSTA_2_2_HV2Y7BCX2.12BA319_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259070/CIN_DYOSTA_2_2_HV2Y7BCX2.12BA338_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259059/CIN_DNOSTA_2_2_HV2Y7BCX2.12BA301_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259070/CIN_DYOSTA_2_1_HV2Y7BCX2.12BA338_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259076/CIN_EEOSTA_2_2_HV2Y7BCX2.12BA315_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259069/CIN_DXOSTA_2_2_HV2Y7BCX2.12BA326_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259076/CIN_EEOSTA_2_1_HV2Y7BCX2.12BA315_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259060/CIN_DOOSTA_2_2_HV2Y7BCX2.12BA313_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259069/CIN_DXOSTA_2_1_HV2Y7BCX2.12BA326_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259060/CIN_DOOSTA_2_1_HV2Y7BCX2.12BA313_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259085/CIN_ENOSTA_2_2_HV2Y7BCX2.12BA328_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259064/CIN_DSOSTA_2_2_HV2Y7BCX2.12BA361_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259072/CIN_EAOSTA_2_2_HV2Y7BCX2.12BA362_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259102/CIN_FEOSTA_2_2_HV2Y7BCX2.12BA342_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259061/CIN_DPOSTA_2_2_HV2Y7BCX2.12BA325_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259068/CIN_DWOSTA_2_1_HV2Y7BCX2.12BA314_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259079/CIN_EHOSTA_2_1_HV2Y7BCX2.12BA351_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259081/CIN_EJOSTA_2_1_HV2Y7BCX2.12BA375_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259104/CIN_FGOSTA_2_1_HV2Y7BCX2.12BA366_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259113/CIN_FPOSTA_2_2_HV2Y7BCX2.12BA379_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259075/CIN_EDOSTA_2_2_HV2Y7BCX2.12BA303_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259095/CIN_EXOSTA_2_1_HV2Y7BCX2.12BA353_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259077/CIN_EFOSTA_2_1_HV2Y7BCX2.12BA327_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259089/CIN_EROSTA_2_1_HV2Y7BCX2.12BA376_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259096/CIN_EYOSTA_2_1_HV2Y7BCX2.12BA365_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259096/CIN_EYOSTA_2_2_HV2Y7BCX2.12BA365_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259105/CIN_FHOSTA_2_2_HV2Y7BCX2.12BA378_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259105/CIN_FHOSTA_2_1_HV2Y7BCX2.12BA378_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259084/CIN_EMOSTA_2_2_HV2Y7BCX2.12BA316_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259084/CIN_EMOSTA_2_1_HV2Y7BCX2.12BA316_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259062/CIN_DQOSTA_2_1_HV2Y7BCX2.12BA337_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259067/CIN_DVOSTA_2_1_HV2Y7BCX2.12BA302_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259099/CIN_FBOSTA_2_1_HV2Y7BCX2.12BA306_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259080/CIN_EIOSTA_2_2_HV2Y7BCX2.12BA363_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259091/CIN_ETOSTA_2_2_HV2Y7BCX2.12BA305_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259106/CIN_FIOSTA_2_1_HV2Y7BCX2.12BA295_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259082/CIN_EKOSTA_2_2_HV2Y7BCX2.12BA292_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259102/CIN_FEOSTA_2_1_HV2Y7BCX2.12BA342_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259078/CIN_EGOSTA_2_2_HV2Y7BCX2.12BA339_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259111/CIN_FNOSTA_2_2_HV2Y7BCX2.12BA355_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259079/CIN_EHOSTA_2_2_HV2Y7BCX2.12BA351_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259063/CIN_DROSTA_2_2_HV2Y7BCX2.12BA349_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259113/CIN_FPOSTA_2_1_HV2Y7BCX2.12BA379_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259073/CIN_EBOSTA_2_1_HV2Y7BCX2.12BA374_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259068/CIN_DWOSTA_2_2_HV2Y7BCX2.12BA314_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259058/CIN_DMOSTA_2_1_HV2Y7BCX2.12BA289_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259087/CIN_EPOSTA_2_1_HV2Y7BCX2.12BA352_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259110/CIN_FMOSTA_2_1_HV2Y7BCX2.12BA343_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259110/CIN_FMOSTA_2_2_HV2Y7BCX2.12BA343_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259093/CIN_EVOSTA_2_1_HV2Y7BCX2.12BA329_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259093/CIN_EVOSTA_2_2_HV2Y7BCX2.12BA329_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259092/CIN_EUOSTA_2_1_HV2Y7BCX2.12BA317_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259071/CIN_DZOSTA_2_1_HV2Y7BCX2.12BA350_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259088/CIN_EQOSTA_2_2_HV2Y7BCX2.12BA364_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259078/CIN_EGOSTA_2_1_HV2Y7BCX2.12BA339_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259111/CIN_FNOSTA_2_1_HV2Y7BCX2.12BA355_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259097/CIN_EZOSTA_2_2_HV2Y7BCX2.12BA377_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259082/CIN_EKOSTA_2_1_HV2Y7BCX2.12BA292_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259083/CIN_ELOSTA_2_1_HV2Y7BCX2.12BA304_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259063/CIN_DROSTA_2_1_HV2Y7BCX2.12BA349_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259098/CIN_FAOSTA_2_2_HV2Y7BCX2.12BA294_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259114/CIN_FQOSTA_2_2_HV2Y7BCX2.12BA296_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259066/CIN_DUOSTA_2_2_HV2Y7BCX2.12BA290_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259109/CIN_FLOSTA_2_1_HV2Y7BCX2.12BA331_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259058/CIN_DMOSTA_2_2_HV2Y7BCX2.12BA289_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259073/CIN_EBOSTA_2_2_HV2Y7BCX2.12BA374_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259087/CIN_EPOSTA_2_2_HV2Y7BCX2.12BA352_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259100/CIN_FCOSTA_2_2_HV2Y7BCX2.12BA318_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259089/CIN_EROSTA_2_2_HV2Y7BCX2.12BA376_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259091/CIN_ETOSTA_2_1_HV2Y7BCX2.12BA305_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259077/CIN_EFOSTA_2_2_HV2Y7BCX2.12BA327_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259062/CIN_DQOSTA_2_2_HV2Y7BCX2.12BA337_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259067/CIN_DVOSTA_2_2_HV2Y7BCX2.12BA302_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259112/CIN_FOOSTA_2_1_HV2Y7BCX2.12BA367_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259112/CIN_FOOSTA_2_2_HV2Y7BCX2.12BA367_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259099/CIN_FBOSTA_2_2_HV2Y7BCX2.12BA306_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259080/CIN_EIOSTA_2_1_HV2Y7BCX2.12BA363_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259074/CIN_ECOSTA_2_2_HV2Y7BCX2.12BA291_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259106/CIN_FIOSTA_2_2_HV2Y7BCX2.12BA295_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259074/CIN_ECOSTA_2_1_HV2Y7BCX2.12BA291_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259086/CIN_EOOSTA_2_1_HV2Y7BCX2.12BA340_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259101/CIN_FDOSTA_2_2_HV2Y7BCX2.12BA330_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259101/CIN_FDOSTA_2_1_HV2Y7BCX2.12BA330_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259086/CIN_EOOSTA_2_2_HV2Y7BCX2.12BA340_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259092/CIN_EUOSTA_2_2_HV2Y7BCX2.12BA317_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259071/CIN_DZOSTA_2_2_HV2Y7BCX2.12BA350_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259085/CIN_ENOSTA_2_1_HV2Y7BCX2.12BA328_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259104/CIN_FGOSTA_2_2_HV2Y7BCX2.12BA366_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259083/CIN_ELOSTA_2_2_HV2Y7BCX2.12BA304_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259108/CIN_FKOSTA_2_1_HV2Y7BCX2.12BA319_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259072/CIN_EAOSTA_2_1_HV2Y7BCX2.12BA362_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259097/CIN_EZOSTA_2_1_HV2Y7BCX2.12BA377_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259088/CIN_EQOSTA_2_1_HV2Y7BCX2.12BA364_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259081/CIN_EJOSTA_2_2_HV2Y7BCX2.12BA375_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259064/CIN_DSOSTA_2_1_HV2Y7BCX2.12BA361_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259098/CIN_FAOSTA_2_1_HV2Y7BCX2.12BA294_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259114/CIN_FQOSTA_2_1_HV2Y7BCX2.12BA296_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259109/CIN_FLOSTA_2_2_HV2Y7BCX2.12BA331_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259066/CIN_DUOSTA_2_1_HV2Y7BCX2.12BA290_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259061/CIN_DPOSTA_2_1_HV2Y7BCX2.12BA325_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259095/CIN_EXOSTA_2_2_HV2Y7BCX2.12BA353_clean.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR112/ERR11259100/CIN_FCOSTA_2_1_HV2Y7BCX2.12BA318_clean.fastq.gz
#
echo Downloading completed !
