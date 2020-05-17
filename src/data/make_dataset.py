# -*- coding: utf-8 -*-
import click
import logging
from pathlib import Path
import os,glob
# from dotenv import find_dotenv, load_dotenv
import process_rna as process_rna
import process_hichip as process_hichip
import process_atac as process_atac
import process_bedtools as process_bedtools
import process_crms as process_crms

##TODO: make workflow for single tissue
###TODO: snakemake this workflow
@click.command()
@click.argument('input_filepath', type=click.Path(exists=True))
@click.argument('interim_filepath', type=click.Path())
@click.argument('output_filepath', type=click.Path())
@click.argument('external_filepath', type=click.Path())
def main(input_filepath, interim_filepath, output_filepath, external_filepath):
    """ Runs data processing scripts to turn raw data from (../raw) into
        cleaned data ready to be analyzed (saved in ../processed).
    """
    logger = logging.getLogger(__name__)
    logger.info('Data Preprocessing')

    logger.info('processing interim RNA')
    # process_rna.run(os.path.join(input_filepath,'rna'), os.path.join(interim_filepath,'rna'))

    logger.info('processing interim ATAC')
    # process_atac.merge_samples(os.path.join(input_filepath,'atac'), os.path.join(interim_filepath,'merged/atac'),
    #                         split_delim = '.narrowPeak.bed', extension=".bed", type_prefix='', verbose=True)

    logger.info('processing interim Footprinting')
    # process_atac.preprocess_footprinting(os.path.join(input_filepath,'footprinting'), os.path.join(interim_filepath,'footprinting'))
    # process_atac.merge_samples(os.path.join(interim_filepath,'footprinting'), os.path.join(interim_filepath,'merged/footprinting'),
    #                         split_delim = '_output_mpbs_pre.bed', extension=".bed", type_prefix='', verbose=True)


    logger.info('processing interim HICHIP')
    # process_hichip.make_bedpe(os.path.join(input_filepath,'hichip'), os.path.join(interim_filepath,'hichip'))
    # process_hichip.make_csvs(os.path.join(interim_filepath,'hichip'), os.path.join(interim_filepath,'merged'))
    # process_bedtools.sort_bed_dir(os.path.join(interim_filepath,'merged/anchors_bed'),os.path.join(interim_filepath,'merged/anchors_bed_sort'))

    logger.info('annotating promoters')
    # promoter_file = os.path.join(external_filepath,'promoter_hg19_2000_500_sort.bed')
    # process_bedtools.annotate_batch(promoter_file, os.path.join(interim_filepath,'merged/atac'),
    #                 os.path.join(interim_filepath,'annon/promoter_atac'),
    #                 extensions=['_hg19_2000_500_sort.bed','_merged.bed'],f=1E-6, wo=False)
    # process_bedtools.annotate_batch(promoter_file, os.path.join(interim_filepath,'merged/footprinting'),
    #                 os.path.join(interim_filepath,'annon/promoter_footprinting'),
    #                 extensions=['_hg19_2000_500_sort.bed','_merged.bed'],f=1E-6, wo=False)
    # logger.info('annotating promoters with anchors that are close by')
    region_search_file = os.path.join(external_filepath,'promoter_hg19_5000_5000_sort.bed')
    process_bedtools.annotate_batch(region_search_file, os.path.join(interim_filepath,'merged/anchors_bed_sort'),
                    os.path.join(interim_filepath,'annon/promoter_anchors'),
                    extensions=['_hg19_5000_5000_sort.bed','_sort.bed'],f=1E-9, wo=True)

    logger.info('annotating anchors')
    # process_bedtools.annotate_batch(os.path.join(interim_filepath,'merged/anchors_bed_sort'),
    #                 os.path.join(interim_filepath,'merged/atac'),
    #                 os.path.join(interim_filepath,'annon/anchor_atac'),
    #                 extensions=['_sort.bed','_merged.bed'],f=2.5E-6)
    # process_bedtools.annotate_batch(os.path.join(interim_filepath,'merged/anchors_bed_sort'),
    #                 os.path.join(interim_filepath,'merged/footprinting'),
    #                 os.path.join(interim_filepath,'annon/anchor_footprinting'),
    #                 extensions=['_sort.bed','_merged.bed'],f=2.5E-6)

    logger.info('creating crms')
    rna_file = os.path.join(interim_filepath, 'rna', 'tissue_tpm_sym.csv')
    tf_annon_file = os.path.join(external_filepath, 'HOCOMOCOv11_annotation.csv')
    annon_file_path = os.path.join(interim_filepath,'annon')
    loop_file_path = os.path.join(interim_filepath, 'merged/loops')
    output_dir = os.path.join(output_filepath,'tissue_crms')
    process_crms.create_crm_batch(rna_file, tf_annon_file, annon_file_path, loop_file_path,
                         output_dir,tissues_sel = ['Airway', 'Pancreas', 'Uterine'],
                         type='all',THRES=1)
    # process_crms.create_crm_batch(rna_file, tf_annon_file, annon_file_path, loop_file_path,
    #                      output_dir,type='all',THRES=1)



if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'


    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]
    print(project_dir)
    # # find .env automagically by walking up directories until it's found, then
    # # load up the .env entries as environment variables
    # load_dotenv(find_dotenv())

    main()
