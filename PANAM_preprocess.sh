#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 30/11/2020
#
#!/bin/bash

# Preprocess and Panam intallation
cd script/
bash 1_Pre-process.sh

# Install Assembler
bash 2_Install_assembler.sh

# Use Vsearch v2.15.1
cd ../dataPANAM
mv PANAM2/bin/vsearch PANAM2/bin/vsearch-1.9.9/bin/
cp ../bin/vsearch-2.15.1/ PANAM2/bin/
mv PANAM2/bin/vsearch-2.15.1/bin/vsearch PANAM2/bin/

# End
echo "Preprocess and assembler installation finish !"
