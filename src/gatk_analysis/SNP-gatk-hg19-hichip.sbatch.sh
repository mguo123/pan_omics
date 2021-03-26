#!/bin/bash
##SNP-gatk-b37-hichip.sbatch
#SBATCH --job-name=gatk-hichip-Colon
#SBATCH --output=/oak/stanford/groups/khavari/users/mguo123/pan_omics/data/processed/hichip_output_snp/GDSD3_gatk.out
#SBATCH --error=/oak/stanford/groups/khavari/users/mguo123/pan_omics/data/processed/hichip_output_snp/GDSD3_gatk.err
#SBATCH --time=24:00:00
#SBATCH -p rbaltman,khavari
#SBATCH --nodes=1
#SBATCH --mem=100000
#SBATCH -c 16
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mguo123@stanford.edu
# source ~/.bashrc
module load biology
ml gatk
ml bwa
module load bowtie2
module load samtools
module load bedtools
module load py-macs2
module load tbb
module load java
module load python
module load perl
module load R

BWA_fasta_ref_file="/oak/stanford/groups/khavari/users/mguo123/reference_genome/b37/human_g1k_v37.fasta"
# BWA_fasta_ref_file="/oak/stanford/groups/khavari/users/mguo123/reference_genome/hg38/Homo_sapiens_assembly38.fasta"
# BWA_fasta_ref_dir="/oak/stanford/groups/khavari/users/mguo123/reference_genome/hg19_UCSC"
BWA_fasta_ref_dir="/oak/stanford/groups/khavari/users/mguo123/reference_genome/b37"
picard_path="/oak/stanford/groups/rbaltman/mguo123/picard/build/libs"

refSNP_path="/oak/stanford/groups/khavari/users/mguo123/refsnps/b37"

# refSNP_path="/oak/stanford/groups/khavari/users/mguo123/refsnps/hg19"
# # make a merged hichip bed file of all the tissue anchor peaks (from FitHiChIP so a bit stringent but that's okay) see OneNote for more details
interval_file="/oak/stanford/groups/khavari/users/mguo123/pan_omics/data/processed/hichip_output_snp/all_epithelial_pks_b37.bed"
# consider using only atac peak regions???
# interval_file="/oak/stanford/groups/khavari/users/mguo123/pan_omics/data/processed/atac_output_snp/all_epithelial_pks.bed"


# SET SAMPLE/celltype INFO
# Celltypes: Airway Astrocytes Bladder Colon Esophageal GDS-D0 GDS-D3 GDS-D6 GM12878 HMEC Melanocytes Ovarian Pancreas Prostate Renal Thyroid Uterine

cell_type="Airway"
sample_list=(Air_B1 Air_B3)
data_path="/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/air_output/bowtie_results/bwt2"
data_path="/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs_022520/bowtie_results/bwt2"

cell_type="Astrocytes"
sample_list=(Astro_B1 Astro_B2)
data_path="/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/blad_output/bowtie_results/bwt2"

sample_list=(Blad_B1 Blad_B2)
cell_type="Bladder"
data_path="/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/astro_output/bowtie_results/bwt2"

sample_list=(Colon_B1 Colon_B2)
cell_type="Colon"
data_path="/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/colon_output/bowtie_results/bwt2"

sample_list=(Eso_B1 Eso_B2)
cell_type="Esophageal"
data_path="/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/eso_output/bowtie_results/bwt2"

sample_list=(GDSD0_B1 GDSD0_B2)
cell_type="GDSD0"
data_path="/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/GDSD0_output/bowtie_results/bwt2"

sample_list=(GDSD3_B1 GDSD3_B2)
cell_type="GDSD3"
data_path="/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/GDSD3_output/bowtie_results/bwt2"

sample_list=(GDSD6_B1 GDSD6_B2)
cell_type="GDSD6"
data_path="/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/GDSD6_output/bowtie_results/bwt2"

sample_list=(GM12878_B1 GM12878_B2)
cell_type="GM12878"
data_path="/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/GM12878_output/bowtie_results/bwt2"

sample_list=(HMEC0_B1 HMEC4_B2)
cell_type="HMEC"
data_path="/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/HMEC_output/bowtie_results/bwt2"

sample_list=(MC_B1 MC_B2)
cell_type="Melanocytes"
data_path="/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/MC_output/bowtie_results/bwt2"

sample_list=(Ova_B1 Ova_B2)
cell_type="Ovarian"
data_path="/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/Ova_output/bowtie_results/bwt2"

#### REVISE
sample_list=(Panc_B2 Panc_B3)
cell_type="Pancreas"
data_path="/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/panc_output/bowtie_results/bwt2"
data_path="/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs_022520/bowtie_results/bwt2"


sample_list=(Pros_B1 Pros_B2)
cell_type="Prostate"
data_path="/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/P$ros_output/bowtie_results/bwt2"


sample_list=(Renal_B1 Renal_B2)
cell_type="Renal"
data_path="/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/Renal_output/bowtie_results/bwt2"

sample_list=(Thy_B1 Thy_B2)
cell_type="Thyroid"
data_path="/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/Thy_output/bowtie_results/bwt2"


## REVISE
sample_list=(Uter_B2 Uter_B3)
cell_type="Uterine"
data_path="/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/uter_output/bowtie_results/bwt2"
data_path="/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/uter_b3_out/bowtie_results/bwt2"

 


### STEP 1A: setting output path
output_path="/oak/stanford/groups/khavari/users/mguo123/pan_omics/data/processed/hichip_output_snp"
output_path=${output_path}/$cell_type
echo $output_path
mkdir -p $output_path

cd $output_path

# for airway premerge samples Manually generated command 
java -jar /oak/stanford/groups/rbaltman/mguo123/picard/build/libs/picard.jar MergeSamFiles O=/oak/stanford/groups/khavari/users/mguo123/pan_omics/data/processed/hichip_output_snp/Airway/Air_B3_merged.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs_022520/bowtie_results/bwt2/Air_B3/00_Air_B3_CKDL200147132-1a-AK1750-AK19316_HTFTCDSXX_L3_hg19.bwt2pairs.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs_022520/bowtie_results/bwt2/Air_B3/01_Air_B3_CKDL200147132-1a-AK1750-AK19316_HTFTCDSXX_L3_hg19.bwt2pairs.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs_022520/bowtie_results/bwt2/Air_B3/02_Air_B3_CKDL200147132-1a-AK1750-AK19316_HTFTCDSXX_L3_hg19.bwt2pairs.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs_022520/bowtie_results/bwt2/Air_B3/03_Air_B3_CKDL200147132-1a-AK1750-AK19316_HTFTCDSXX_L3_hg19.bwt2pairs.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs_022520/bowtie_results/bwt2/Air_B3/04_Air_B3_CKDL200147132-1a-AK1750-AK19316_HTFTCDSXX_L3_hg19.bwt2pairs.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs_022520/bowtie_results/bwt2/Air_B3/05_Air_B3_CKDL200147132-1a-AK1750-AK19316_HTFTCDSXX_L3_hg19.bwt2pairs.bam
java -jar /oak/stanford/groups/rbaltman/mguo123/picard/build/libs/picard.jar MergeSamFiles O=/oak/stanford/groups/khavari/users/mguo123/pan_omics/data/processed/hichip_output_snp/Airway/Air_B1_merged.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/air_output/bowtie_results/bwt2/Air_B1/Air_B1_USPD16099870-N704-AK402_HL7H3DSXX_L1_hg19.bwt2pairs.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/air_output/bowtie_results/bwt2/Air_B1/Air_B1_USPD16099870-N704-AK402_HL7H3DSXX_L2_hg19.bwt2pairs.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/air_output/bowtie_results/bwt2/Air_B1/Air_B1_USPD16099870-N704-AK402_HL7H3DSXX_L3_hg19.bwt2pairs.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/air_output/bowtie_results/bwt2/Air_B1/Air_B1_USPD16099870-N704-AK402_HL7H3DSXX_L4_hg19.bwt2pairs.bam

# for panc premerge samples Manually generated command 
java -jar /oak/stanford/groups/rbaltman/mguo123/picard/build/libs/picard.jar MergeSamFiles O=/oak/stanford/groups/khavari/users/mguo123/pan_omics/data/processed/hichip_output_snp/Pancreas/Panc_B2_merged.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/panc_output/bowtie_results/bwt2/Panc_B2/Panc_B2_USPD16099870-AK11822-AK13225_HL7H3DSXX_L1_hg19.bwt2pairs.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/panc_output/bowtie_results/bwt2/Panc_B2/Panc_B2_USPD16099870-AK11822-AK13225_HL7H3DSXX_L2_hg19.bwt2pairs.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/panc_output/bowtie_results/bwt2/Panc_B2/Panc_B2_USPD16099870-AK11822-AK13225_HL7H3DSXX_L3_hg19.bwt2pairs.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/panc_output/bowtie_results/bwt2/Panc_B2/Panc_B2_USPD16099870-AK11822-AK13225_HL7H3DSXX_L4_hg19.bwt2pairs.bam
java -jar /oak/stanford/groups/rbaltman/mguo123/picard/build/libs/picard.jar MergeSamFiles O=/oak/stanford/groups/khavari/users/mguo123/pan_omics/data/processed/hichip_output_snp/Pancreas/Panc_B3_merged.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs_022520/bowtie_results/bwt2/Panc_B3/00_Panc_B3_CKDL200147132-1a-AK11823-AK19316_HTFTCDSXX_L3_hg19.bwt2pairs.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs_022520/bowtie_results/bwt2/Panc_B3/01_Panc_B3_CKDL200147132-1a-AK11823-AK19316_HTFTCDSXX_L3_hg19.bwt2pairs.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs_022520/bowtie_results/bwt2/Panc_B3/02_Panc_B3_CKDL200147132-1a-AK11823-AK19316_HTFTCDSXX_L3_hg19.bwt2pairs.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs_022520/bowtie_results/bwt2/Panc_B3/03_Panc_B3_CKDL200147132-1a-AK11823-AK19316_HTFTCDSXX_L3_hg19.bwt2pairs.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs_022520/bowtie_results/bwt2/Panc_B3/04_Panc_B3_CKDL200147132-1a-AK11823-AK19316_HTFTCDSXX_L3_hg19.bwt2pairs.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs_022520/bowtie_results/bwt2/Panc_B3/05_Panc_B3_CKDL200147132-1a-AK11823-AK19316_HTFTCDSXX_L3_hg19.bwt2pairs.bam

# for uter premerge samples Manually generated command
java -jar /oak/stanford/groups/rbaltman/mguo123/picard/build/libs/picard.jar MergeSamFiles O=/oak/stanford/groups/khavari/users/mguo123/pan_omics/data/processed/hichip_output_snp/Uterine/Uter_B2_merged.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/uter_output/bowtie_results/bwt2/Uter_B2/Uter_B2_USPD16099870-N708-AK403_HL7H3DSXX_L1_hg19.bwt2pairs.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/uter_output/bowtie_results/bwt2/Uter_B2/Uter_B2_USPD16099870-N708-AK403_HL7H3DSXX_L2_hg19.bwt2pairs.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/uter_output/bowtie_results/bwt2/Uter_B2/Uter_B2_USPD16099870-N708-AK403_HL7H3DSXX_L3_hg19.bwt2pairs.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/uter_output/bowtie_results/bwt2/Uter_B2/Uter_B2_USPD16099870-N708-AK403_HL7H3DSXX_L4_hg19.bwt2pairs.bam
java -jar /oak/stanford/groups/rbaltman/mguo123/picard/build/libs/picard.jar MergeSamFiles O=/oak/stanford/groups/khavari/users/mguo123/pan_omics/data/processed/hichip_output_snp/Uterine/Uter_B3_merged.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/uter_b3_out/bowtie_results/bwt2/Uter_B3/00_Uter_B3_CKDL200147132-1a-AK11822-AK19316_HTFTCDSXX_L3_hg19.bwt2pairs.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/uter_b3_out/bowtie_results/bwt2/Uter_B3/01_Uter_B3_CKDL200147132-1a-AK11822-AK19316_HTFTCDSXX_L3_hg19.bwt2pairs.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/uter_b3_out/bowtie_results/bwt2/Uter_B3/02_Uter_B3_CKDL200147132-1a-AK11822-AK19316_HTFTCDSXX_L3_hg19.bwt2pairs.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/uter_b3_out/bowtie_results/bwt2/Uter_B3/03_Uter_B3_CKDL200147132-1a-AK11822-AK19316_HTFTCDSXX_L3_hg19.bwt2pairs.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/uter_b3_out/bowtie_results/bwt2/Uter_B3/04_Uter_B3_CKDL200147132-1a-AK11822-AK19316_HTFTCDSXX_L3_hg19.bwt2pairs.bam I=/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/uter_b3_out/bowtie_results/bwt2/Uter_B3/05_Uter_B3_CKDL200147132-1a-AK11822-AK19316_HTFTCDSXX_L3_hg19.bwt2pairs.bam

for sample in ${sample_list[@]}; do 
   ## FOR EVERYTHING EXCEPT ASTROCYTES
# sample=${sample_list[${idx}]}
echo $sample




### STEP 2 merge the bam files per sample of a tissue (only done if data_path exists so not Airway, Pancreas, Uterine)
# bam_files=`ls ${data_path}/${sample}/*bwt2pairs.bam`

# merge_cmd="java -jar ${picard_path}/picard.jar MergeSamFiles O=${output_path}/${sample}_merged.bam"
# for f in $bam_files; do 
# echo $f
# merge_cmd="${merge_cmd} I=${f}"
# done
# echo $merge_cmd
# eval $merge_cmd
# echo 'done merging bam files of sample: ${sample}...'
###### done merging bam file (for not Airway, Pancreas, Uterine)

 
# java -jar $picard_path/picard.jar SortSam INPUT=${output_path}/${sample}_merged.bam OUTPUT=$output_path/${sample}_sorted.bam SORT_ORDER=coordinate
# java -jar $picard_path/picard.jar MarkDuplicates \
#       I=$output_path/${sample}_sorted.bam \
#       O=$output_path/${sample}_sorted_mkdup.bam \
#       M=$output_path/${sample}_marked_dup_metrics.txt
# java -jar $picard_path/picard.jar BuildBamIndex INPUT=$output_path/${sample}_sorted_mkdup.bam

java -jar $picard_path/picard.jar AddOrReplaceReadGroups \
    I=$output_path/${sample}_sorted_mkdup.bam   \
    O=$output_path/${sample}_sorted_mkdup_rg.bam  \
    SORT_ORDER=coordinate \
    RGID=psych \
    RGLB=sample \
    RGPU=sample \
    RGPL=illumina \
    RGSM=${sample}
java -jar $picard_path/picard.jar BuildBamIndex INPUT=$output_path/${sample}_sorted_mkdup_rg.bam



## DEBUG CHECK Readgroups
# samtools view -H $output_path/${sample}_sorted_mkdup.bam  | grep '@RG'




# STEP 3: base recalibration takes 5-10 min
#First pass of the Base Quality Score Recalibration (BQSR) -- Generates recalibration table based on various user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context).
gatk BaseRecalibrator \
-R $BWA_fasta_ref_file \
-I $output_path/${sample}_sorted_mkdup_rg.bam \
-O $output_path/${sample}_recal.table \
--known-sites $refSNP_path/dbsnp_138.b37.vcf \
--known-sites $refSNP_path/1000G_phase1.snps.high_confidence.b37.vcf \
--known-sites $refSNP_path/1000G_phase1.indels.b37.vcf \
--known-sites $refSNP_path/Mills_and_1000G_gold_standard.indels.b37.vcf \
-L $interval_file


#  STEP 4: Apply base quality score recalibration (takes ~3-5 min)
gatk ApplyBQSR \
-I $output_path/${sample}_sorted_mkdup.bam \
-bqsr $output_path/${sample}_recal.table \
--reference $BWA_fasta_ref_file \
-L $interval_file \
-O $output_path/${sample}_recal.bam 
# java -jar $picard_path/picard.jar AddOrReplaceReadGroups \
#     I=$output_path/${sample}_recal.bam  \
#     O=$output_path/${sample}_recal_rg.bam \
#     SORT_ORDER=coordinate \
#     RGID=psych \
#     RGLB=sample \
#     RGPU=sample \
#     RGPL=illumina \
#     RGSM=${sample} 
# java -jar $picard_path/picard.jar BuildBamIndex INPUT=$output_path/${sample}_recal_rg.bam



# STEP 5: snp calling for single bam file sample
# allele specific version : https://gatk.broadinstitute.org/hc/en-us/articles/360035890551-Allele-specific-annotation-and-filtering-of-germline-short-variants
# 1. Variant  (20 min) get gvcf
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531812?id=4017 - GVCF vs VCF

gatk HaplotypeCaller  \
--reference $BWA_fasta_ref_file \
-I $output_path/${sample}_recal.bam \
-O $output_path/${sample}.g.vcf.gz \
-L $interval_file \
-ERC GVCF  \
   -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation
 #
# --annotation-group,-G:String  One or more groups of annotations to apply to variant calls  This argument may be
#                               specified 0 or more times. Default value: null. Possible Values:
#                               {AlleleSpecificAnnotation, AS_StandardAnnotation, ReducibleAnnotation, StandardAnnotation,
#                               StandardHCAnnotation, StandardMutectAnnotation}


done




# STEP6. Data aggregation step ( for hichip)
# combine gvcfs
gatk CombineGVCFs \
--reference $BWA_fasta_ref_file \
--variant $output_path/${sample_list[0]}.g.vcf.gz \
--variant $output_path/${sample_list[1]}.g.vcf.gz \
-O $output_path/${cell_type}.g.vcf.gz
 -G StandardAnnotation -G AS_StandardAnnotation



# STEP7. Joint genotyping
# Take the outputs from step 2 (or step 1 if dealing with fewer samples) and run GenotypeGVCFs on all of them together to create the raw SNP and indel VCFs that are usually emitted by the callers.
gatk GenotypeGVCFs \
   -R $BWA_fasta_ref_file \
   -V $output_path/${cell_type}.g.vcf.gz \
   -O $output_path/${cell_type}.vcf.gz\
   -G StandardAnnotation -G AS_StandardAnnotation
 

# STEP 8 . Variant recalibration
# vcf spec https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
# https://gatk.broadinstitute.org/hc/en-us/articles/360037594511-VariantRecalibrator#top
gatk VariantRecalibrator \
   -R $BWA_fasta_ref_file \
   -V $output_path/${cell_type}.vcf.gz \
   -AS \
  --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${refSNP_path}/hapmap_3.3.b37.vcf \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 ${refSNP_path}/1000G_phase1.indels.b37.vcf \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${refSNP_path}/1000G_phase1.snps.high_confidence.b37.vcf \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${refSNP_path}/dbsnp_138.b37.vcf \
   -an DP -an QD -an FS -an SOR -an MQ -an ReadPosRankSum \
   -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
   -mode SNP \
   -O $output_path/${cell_type}.AS.recal \
   --tranches-file $output_path/${cell_type}.AS.tranches \
   --rscript-file $output_path/${cell_type}.plots.AS.R
 
 gatk ApplyVQSR \
   -R $BWA_fasta_ref_file \
   -V $output_path/${cell_type}.vcf.gz \
   -O $output_path/${cell_type}_postvqsr.vcf.gz \
   --truth-sensitivity-filter-level 99.0 \
   --tranches-file $output_path/${cell_type}.AS.tranches \
   --recal-file $output_path/${cell_type}.AS.recal \
   -mode SNP



