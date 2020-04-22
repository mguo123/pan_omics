#!/bin/bash
# Margaret Guo
# pan-omics
# 04/08/2020 - deprecated 04/20/2020
# process_anchors.sh
# for merging anchor files
# Sort bed files with bedtools sort
# start in src ##TODO FIX
# cd ../..
# parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
# echo $parent_path
# ATAC_DIR="$parent_path"/data/raw/atac
# FOOT_DIR="$parent_path"/data/raw/footprinting
#
# OUT_DIR="$parent_path"/data/interim
# OUT_ATAC_DIR="$OUT_DIR"/atac
# OUT_FOOT_DIR="$OUT_DIR"/dfootprinting
# OUT_MERGE_DIR="$OUT_DIR"/merged



# cd /Users/mguo123/Google Drive/1_khavari/omics_project-LD/diffloop_data/bedpe_files_csvs_rna
# for f in `ls *bed`; do
# echo $f;
# tissue=${f%%-*};
# echo $tissue;
# bedtools sort -i $f > ${tissue}_anchors_sort.bed;
# done
#
# cd /Users/mguo123/Google Drive/1_khavari/omics_project-LD/atac_footprinting/merged_into_anchors_hoco`
# cp ../../diffloop_data/bedpe_files_csvs_rna/*sort.bed .`

# # Remove tab from end of footprinting files and then bedtools sort
# cd /Users/mguo123/Google Drive/1_khavari/omics_project-LD/atac_footprinting/match_hoco
# #### sed 's/[[:blank:]]*$//'`
# for f in `ls *bed`; do
# echo $f; tissue=${f%%_*};echo $tissue;
# f2=${tissue}_footprinting_hoco.bed;
# f3=${tissue}_footprinting_hoco_sort.bed;
# sed 's/[[:blank:]]*$//' $f > $f2;
# bedtools sort -i $f2 > $f3;
# bedtools merge ...
# done
#
# mv *footprinting_hoco_sort.bed ../merged_into_anchors_hoco


# # Intersect tissue specific files
# cd /Users/mguo123/Google Drive/1_khavari/omics_project-LD/atac_footprinting/merged_into_anchors_hoco
#   b. `mv Melanocytes_footprinting_hoco_sort.bed MC_footprinting_hoco_sort.bed`
#   c. `mv GDS-D0_footprinting_hoco_sort.bed GDSD0_footprinting_hoco_sort.bed`
#   d. `mv GDS-D3_footprinting_hoco_sort.bed GDSD3_footprinting_hoco_sort.bed`
#   e. `mv GDS-D6_footprinting_hoco_sort.bed GDSD6_footprinting_hoco_sort.bed`
#   f. `####bedtools intersect -wa -wb -a Air_anchors_sort.bed -b Airway_footprinting_sorted.bed -sorted -names footprinting -f 2.5E-6 > Airway_anchors_footprinting.bed`
# for f in `ls *anchors_sort.bed`; do
# echo $f; tissue=${f%%_*}; echo $tissue;
# bfile=`find ${tissue}*footprinting_hoco_sort.bed`; echo $bfile;
# newfile=${tissue}_anchors_footprinting.bed;
# bedtools intersect -wa -wb -a $f -b $bfile -sorted -names footprinting -f 2.5E-6 > $newfile;
# done
