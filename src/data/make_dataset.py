# -*- coding: utf-8 -*-
import click
import logging
from pathlib import Path
import os
# from dotenv import find_dotenv, load_dotenv
import process_rna as process_rna
import process_hichip_1 as process_hichip_1

@click.command()
@click.argument('input_filepath', type=click.Path(exists=True))
@click.argument('interim_filepath', type=click.Path())
@click.argument('output_filepath', type=click.Path())
def main(input_filepath, interim_filepath, output_filepath):
    """ Runs data processing scripts to turn raw data from (../raw) into
        cleaned data ready to be analyzed (saved in ../processed).
    """
    logger = logging.getLogger(__name__)
    logger.info('making final data set from raw data')

    logger.info('processing interim RNA')
    # process_rna.run(os.path.join(input_filepath,'rna'), os.path.join(interim_filepath,'rna'))

    logger.info('processing interim HICHIP_DIR')
    # process_hichip_1.make_bedpe(os.path.join(input_filepath,'hichip'), os.path.join(interim_filepath,'hichip'))
    process_hichip_1.make_csvs(os.path.join(interim_filepath,'hichip'), os.path.join(interim_filepath,'merged'))






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
