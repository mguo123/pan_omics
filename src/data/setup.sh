#!/bin/bash
# Margaret Guo
# pan-omics
# 04/08/2020
# how to setup the data

# 1. Figure out which tissues you want
# edit `tissues.txt` make sure the hyphens and underscores, capitalization is correct

# 2. make directories in data raw folder
cd ../..
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
echo $parent_path
DATA_DIR="$parent_path"/data/raw
HICHIP_DIR=$DATA_DIR/hichip
ATAC_DIR=$DATA_DIR/atac
FOOT_DIR=$DATA_DIR/footprinting
RNA_DIR=$DATA_DIR/rna

mkdir -p "$HICHIP_DIR"
mkdir "$ATAC_DIR"
mkdir "$FOOT_DIR"
mkdir "$RNA_DIR"

# 3. print tissue
echo "tissue folders to make... "
#cat tissues.txt
cd src/data
for tissue in `cat tissues.txt`;
do echo $tissue;
mkdir "$HICHIP_DIR/$tissue"
mkdir "$ATAC_DIR/$tissue"
mkdir "$FOOT_DIR/$tissue"
mkdir "$RNA_DIR/$tissue"
done
