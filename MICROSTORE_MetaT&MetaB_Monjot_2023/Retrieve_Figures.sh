#!/bin/bash
#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 25/09/2023
#
BEFORE=$SECONDS
# R CStack
ulimit -s unlimited
# Activate conda environment
PATHCONDA=$(conda info | grep -i 'base environment' | awk -F" " '{print $4}')
source $PATHCONDA'/etc/profile.d/conda.sh'
conda activate REnv_Monjot_2023A
if [ $(echo $1 | grep "." |wc -l ) == 1 ]
then
OUTPUT=$1
else
echo 'Enter result file path : '
read OUTPUT
fi

echo "Output: "$OUTPUT


# set script directory
mkdir Monjot_etal_2023
mkdir Monjot_etal_2023/Figures
mkdir Monjot_etal_2023/Figures_Supp
#Fig 2
cp dataDADA2/result/$OUTPUT/AFC-Distribution/AFC-Sequence-90-article.svg Monjot_etal_2023/Figures/Figure_2.svg
#Fig 3
cp dataDADA2/result/$OUTPUT/Composition/Table-Total-Full.svg Monjot_etal_2023/Figures/Figure_3.svg
#Fig 4
cp dataDADA2/result/$OUTPUT/Functional-Analyse/Composition_Function/Totalcompo.svg Monjot_etal_2023/Figures/Figure_4.svg
#Fig 5
cp dataDADA2/result/$OUTPUT/Functional-Analyse/Tree2_Function/ZoneXFractionXGenus_ASV_Tree.svg Monjot_etal_2023/Figures/Figure_5A.svg
cp dataDADA2/result/$OUTPUT/Functional-Analyse/Tree2_Function/ZoneXFractionXGenus_Sequence_Tree.svg Monjot_etal_2023/Figures/Figure_5B.svg
#Fig 6
cp dataDADA2/result/$OUTPUT/Functional-Analyse/Hist_Function/ZoneXPeriods-Full.svg Monjot_etal_2023/Figures/Figure_6.svg
#Fig 7
cp dataDADA2/result/$OUTPUT/Metatrans/DESeq2/Quaternaryplot/Quaterneryplot_KO_1e-02_log2_2_withlabel_only_10first.svg Monjot_etal_2023/Figures/Figure_7.svg
#Fig 8
cp dataDADA2/result/$OUTPUT/Metatrans/DESeq2/Heatmap_raw/Heatmap_Zone_KO_5e-02_log2_2_up10_morethan5sharingko.svg Monjot_etal_2023/Figures/Figure_8.svg
#
#Fig S1
cp dataDADA2/result/$OUTPUT/Stat-Analyse/Analyse-SumAvRare-ASV-art.svg Monjot_etal_2023/Figures_Supp/Figure_S1A.svg
cp dataDADA2/result/$OUTPUT/Stat-Analyse/Analyse-SumAvRare-Sequences-art.svg Monjot_etal_2023/Figures_Supp/Figure_S1B.svg
#Fig S2
cp dataDADA2/result/$OUTPUT/Functional-Analyse/Group_FUNCTION/PCoA_All.svg Monjot_etal_2023/Figures_Supp/Figure_S2.svg
#Fig S3
cp dataDADA2/result/$OUTPUT/Functional-Analyse/Group_FUNCTION/Correlation.svg Monjot_etal_2023/Figures_Supp/Figure_S3.svg
#Fig S4
cp dataDADA2/result/$OUTPUT/Functional-Analyse/Group_FUNCTION/Kmeans-ssi.svg Monjot_etal_2023/Figures_Supp/Figure_S4.svg
#Fig S5
cp dataDADA2/result/$OUTPUT/Functional-Analyse/Group_FUNCTION/PCoA_Kmean-ssi-DIM12.svg Monjot_etal_2023/Figures_Supp/Figure_S5.svg
#Fig S6
cp dataDADA2/result/$OUTPUT/Functional-Analyse/Group_FUNCTION/Def_group_by_factor.svg Monjot_etal_2023/Figures_Supp/Figure_S6.svg
#Fig S7
cp dataDADA2/result/$OUTPUT/Rarecurve/Rarecurve-Raw.pdf Monjot_etal_2023/Figures_Supp/Figure_S7.pdf
#Fig S8
cp dataDADA2/result/$OUTPUT/Diversity/All_article_Diversity.svg Monjot_etal_2023/Figures_Supp/Figure_S8.svg
#Fig S9
cp dataDADA2/result/$OUTPUT/Stat-Analyse/ASV-Stat-FractionXZone.svg Monjot_etal_2023/Figures_Supp/Figure_S9.svg
#Fig S10
cp dataDADA2/result/$OUTPUT/Stat-Analyse/Stat-Taxo.svg Monjot_etal_2023/Figures_Supp/Figure_S10.svg
#Fig S11
cp dataDADA2/result/$OUTPUT/Functional-Analyse/Composition_Function/Polar-Total-asv-grp.svg Monjot_etal_2023/Figures_Supp/Figure_S11.svg
#Fig S12
cp dataDADA2/result/$OUTPUT/Functional-Analyse/Tree2_Function/PeriodsXGenus_All_Tree.svg Monjot_etal_2023/Figures_Supp/Figure_S12.svg
#Fig S13
cp dataDADA2/result/$OUTPUT/Metatrans/DESeq2/Heatmap_raw/Heatmap_Zone_KO_5e-02_log2_2.svg Monjot_etal_2023/Figures_Supp/Figure_S13.svg
#Fig S14
cp dataDADA2/result/$OUTPUT/Metatrans/DESeq2/Heatmap_raw/Heatmap_Fraction_KO_5e-02_log2_2.svg Monjot_etal_2023/Figures_Supp/Figure_S14.svg
#Fig S15
cp dataDADA2/result/$OUTPUT/Metatrans/DESeq2/Heatmap_raw/Heatmap_Periods_04vs06_KO_5e-02_log2_2.svg Monjot_etal_2023/Figures_Supp/Figure_S15A.svg
cp dataDADA2/result/$OUTPUT/Metatrans/DESeq2/Heatmap_raw/Heatmap_Periods_04vs09_KO_5e-02_log2_2.svg Monjot_etal_2023/Figures_Supp/Figure_S15B.svg
cp dataDADA2/result/$OUTPUT/Metatrans/DESeq2/Heatmap_raw/Heatmap_Periods_04vs11_KO_5e-02_log2_2.svg Monjot_etal_2023/Figures_Supp/Figure_S15C.svg
cp dataDADA2/result/$OUTPUT/Metatrans/DESeq2/Heatmap_raw/Heatmap_Periods_06vs09_KO_5e-02_log2_2.svg Monjot_etal_2023/Figures_Supp/Figure_S15D.svg
cp dataDADA2/result/$OUTPUT/Metatrans/DESeq2/Heatmap_raw/Heatmap_Periods_06vs11_KO_5e-02_log2_2.svg Monjot_etal_2023/Figures_Supp/Figure_S15E.svg
cp dataDADA2/result/$OUTPUT/Metatrans/DESeq2/Heatmap_raw/Heatmap_Periods_09vs11_KO_5e-02_log2_2.svg Monjot_etal_2023/Figures_Supp/Figure_S15F.svg
#Fig S16
cp dataDADA2/result/$OUTPUT/Metatrans/DESeq2/Heatmap_raw/Bubble_Mix_KO_5e-02_log2_2.svg Monjot_etal_2023/Figures_Supp/Figure_S16.svg
#Fig S17
cp dataDADA2/result/$OUTPUT/Metatrans/DESeq2/Heatmap_raw/Bubble_Periods_KO_5e-02_log2_2.svg Monjot_etal_2023/Figures_Supp/Figure_S17.svg
#

ELAPSED=$((($SECONDS-$BEFORE)/60))
echo "R process finished and takes "$ELAPSED" minutes !"
