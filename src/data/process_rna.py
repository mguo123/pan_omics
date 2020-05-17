"""
process_rna.py

Make TPM/FPKM matrix

4/09/2020 - from `make_rna_mat.ipynb` 12/19/2019

"""
import subprocess
import os, glob
import pandas as pd
from collections import defaultdict, Counter

def run(input_filepath, output_filepath, extension=".genes.results"):
    """
    Make TPM/FPKM matrix

    """
    print('processing RNA')
    if not os.path.exists(output_filepath):
        os.makedirs(output_filepath)

    tissue_to_sample_tpm = defaultdict(dict)
    tissue_to_sample_fpkm = defaultdict(dict)

    tissue_tpm_dict = {}
    tissue_fpkm_dict = {}
    sample_tpm_dict = {}
    sample_fpkm_dict = {}


    # read file
    for subdir, dirs, files in os.walk(input_filepath):
        for filename in files:
            tissue = os.path.basename(subdir)
            filepath = subdir + os.sep + filename

            if filepath.endswith(extension) :
                sample = filename.split(extension)[0]
                print ('extracted', tissue, sample)

                df = pd.read_table(filepath, header=0)
                df['gene_id_short']=df['gene_id'].str.split(".",expand=True).loc[:,0]

                gene_to_tpm = pd.Series(df.TPM.values, index=df.gene_id_short.values).to_dict()
                gene_to_fpkm = pd.Series(df.FPKM.values, index=df.gene_id_short.values).to_dict()
                tissue_to_sample_tpm[tissue][sample] = gene_to_tpm
                tissue_to_sample_fpkm[tissue][sample] = gene_to_fpkm
                sample_tpm_dict[sample] = gene_to_tpm
                sample_fpkm_dict[sample] = gene_to_fpkm

    # tissue_to_sample_tpm
    for tissue in tissue_to_sample_tpm.keys():
        mean_gene_tpm = pd.Series(pd.DataFrame(tissue_to_sample_tpm[tissue]).mean(axis=1)).to_dict()
        mean_gene_fpkm = pd.Series(pd.DataFrame(tissue_to_sample_fpkm[tissue]).mean(axis=1)).to_dict()
        tissue_tpm_dict[tissue] = mean_gene_tpm
        tissue_fpkm_dict[tissue] = mean_gene_fpkm

    # tpm_df = pd.DataFrame(tissue_tpm_dict)
    # tpm_df = tpm_df.reindex(sorted(tpm_df.columns), axis=1)
    #
    # fpkm_df = pd.DataFrame(tissue_fpkm_dict)
    # fpkm_df = fpkm_df.reindex(sorted(fpkm_df.columns), axis=1)
    #
    # tpm_df.to_csv(os.path.join(output_filepath,'tissue_tpm.csv'))
    # fpkm_df.to_csv(os.path.join(output_filepath,'tissue_fpkm.csv'))

    make_rna_df(tissue_tpm_dict, os.path.join(output_filepath,'tissue_tpm.csv'))
    make_rna_df(tissue_fpkm_dict, os.path.join(output_filepath,'tissue_fpkm.csv'))
    make_rna_df(sample_tpm_dict, os.path.join(output_filepath,'sample_tpm.csv'))
    make_rna_df(sample_fpkm_dict, os.path.join(output_filepath,'sample_fpkm.csv'))

    ###TODO: figure out how to make this silent
    subprocess.call (["Rscript", "--vanilla",
            "src/data/process_rna.R",
            output_filepath])

def make_rna_df(key_rnaval_dict, save_path=None):
    df = pd.DataFrame(key_rnaval_dict)
    df = df.reindex(sorted(df.columns), axis=1)
    if save_path is not None:
        df.to_csv(save_path)
