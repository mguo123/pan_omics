"""
Margaret Guo

Statistical Functions for Fisher's exact to get the pairwise TF values


"""

import pandas as pd
import glob, os
from collections import Counter,defaultdict
import seaborn as sns
import numpy as np
import scipy.stats as stats

def eval_tf(tf, genes, original_df, pseudo=1,verbose=False):
    """
    tf = TF column name, i.e. 'TBL1XR1'
    genes = list of genes significant in CRISPR experiment
    original_df = gene (rows) x experiment/TFs (columns), # TF footprints

    matrix:
        foots_0_0 = observed # of TF footprints in crispr sig genes for specified TF
        m_0_2 = # of TF footprints in crispr sig genes across all chipseq TF KO experiments (with TF KO's shared btn crispr and chipseq)
        foots_0_1 = # of TF footprints in crispr sig genes across all other chipseq TF KO experiments
        foots_1_0 = # of TF footprints for specified chipseq TF KO experiments, excluding those found in crispr sig genes
        foots_1_1 = # of TF footprints not found in specifed chip seq experiment, and not found in



    """
    # get marginals
    row_sum = original_df.sum(axis=1).to_dict()
    col_sum = original_df.sum(axis=0).to_dict()
    tot_sum = sum(row_sum.values())

    # get number of footprints in crispr sig genes for specificied TF
    foots_0_0 =  sum(original_df[tf][original_df.index.isin(genes)])

    # get number of footprints for all genes which were crispr sig, over all experiments
    m_0_2 = 0
    for g in genes:
        try:
            m_0_2 +=row_sum.get(g)
        except:
            continue # can't find gene

    foots_0_1 = m_0_2 - foots_0_0
    foots_1_0 = col_sum[tf] - foots_0_0
    t_min_m_1_2 = tot_sum - m_0_2
    foots_1_1 = t_min_m_1_2 - foots_1_0

    observed_num = foots_0_0
    expected_num = m_0_2*col_sum[tf]/tot_sum


    mat_genes = np.array([[foots_0_0, foots_0_1],[foots_1_0,foots_1_1]]).reshape((2,2))
#     oddsratio, pvalue = stats.fisher_exact(mat_genes,alternative='greater')
    mat_genes_pseudo = mat_genes+pseudo
    _, pvalue_pseudo = stats.fisher_exact(mat_genes_pseudo,alternative='greater')
    oddsratio_pseudo = mat_genes_pseudo[0][0]*mat_genes_pseudo[1][1]/(mat_genes_pseudo[1][0]*mat_genes_pseudo[0][1])
    if verbose:
        if oddsratio_pseudo>1:
            print(mat_genes_pseudo)
            print(tf, oddsratio_pseudo, pvalue_pseudo)
    return observed_num, expected_num, oddsratio_pseudo, pvalue_pseudo


def calculate_enrichment(foot_count_df,type_count_df, seed_TF=None, pseudo=1,verbose=False):
    """
    calls eval_tf

    ###TODO fill in
    Arguments:

    Outputs
    """
    foot_count_df_wseedTF = foot_count_df[foot_count_df[seed_gene]>0]
    foot_count_df_wseedTF = foot_count_df_wseedTF.loc[:,foot_count_df_wseedTF.sum(axis=0)>0]
    genes = foot_count_df_wseedTF.index

    print('count table shape pre seed TF filter', seed_TF, foot_count_df.shape)
    print('count table shape post seed TF filter', seed_TF, foot_count_df_wseedTF.shape)

    results = pd.DataFrame()
    results['num_sig_genes']=len(genes)
    results['obs_num_footprints']=0
    results['exp_num_footprints']=0
    results['odds_ratio']=0.0
    results['pval']=0.0

    for tf in foot_count_df_wseedTF.columns:
        observed_num, expected_num, oddsratio, pvalue = eval_tf(tf, genes, foot_count_df,
                                                                pseudo=pseudo, verbose=verbose)
#         results.at[tf, 'TF'] = tf
        results.at[tf, 'num_sig_genes'] = len(genes)
        results.at[tf, 'obs_num_footprints'] = observed_num
        results.at[tf, 'exp_num_footprints'] = expected_num
        results.at[tf, 'odds_ratio'] = oddsratio
        results.at[tf, 'pval'] = pvalue

    results.sort_index(inplace=True)
    results['log_odds_ratio'] = np.log2(results.odds_ratio)
    results['pval_bonf'] =  results.pval.apply(lambda x: min(1, x* results.shape[0]))
    results['type']=type_count_df
    results['seed_TF']=seed_TF
    print('num sig tfs', results[results.pval_bonf<0.05].shape[0],results[results.pval<0.05].shape[0])
    return results
