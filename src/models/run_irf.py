# basic packages
import os, glob
# import argparse
import pandas as pd
import numpy as np; np.random.seed(0)
import itertools
import time
from collections import Counter

# machine learning packages from sklearn
from sklearn.preprocessing import MinMaxScaler #StandardScaler
from sklearn import preprocessing, metrics
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import train_test_split#, KFold, cross_validate, cross_val_score, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor, VotingClassifier, AdaBoostClassifier
from sklearn.metrics import roc_auc_score, auc, roc_curve, plot_roc_curve, confusion_matrix, accuracy_score
from scipy import interp

# for IRF
from functools import reduce
# Needed for the scikit-learn wrapper function
import irf
from irf import irf_utils, utils, irf_jupyter_utils
from irf.ensemble.wrf import RandomForestClassifierWithWeights
from math import ceil

# Import our custom utilities
from imp import reload
# Import tools needed for visualization
import seaborn as sns; sns.set()
import matplotlib
import matplotlib.pyplot as plt

tissues = ['Airway','Astrocytes','Bladder','Colon', 'Esophageal', 'GDSD0', 'GDSD3', 'GDSD6', 'GM12878', 'HMEC',
 'Melanocytes', 'Ovarian', 'Pancreas', 'Prostate', 'Renal', 'Thyroid', 'Uterine']

groups = ['purple', 'blue', 'purple', 'green', 'green', 'purple', 'purple', 'purple', 'grey', 'purple', 'blue', 'green', 'green', 'purple', 'green', 'green', 'purple']

tissue_to_group = dict(zip(tissues, groups))

def load_data(tissues):
    # import data
    data_all = pd.read_csv('/Users/mguo123/Google Drive/1_khavari/omics_project-LD/pan_omics/data/processed/tissue_crms/all_count_comb_overall.csv',index_col=0,header=0)
    data_all = data_all[data_all.tissue.isin(tissues)]
    data_all = data_all[data_all.iloc[:,2:].sum(axis=1)>1e-1]

    # expression labels
    # exp_label = list(data_all.exp.values)
    labels_all  = data_all.tissue.values#np.array((data_all.exp>THRES).values)

    genes_all = data_all.index.values
    gene_to_num_dict = dict(zip(np.unique(genes_all),range(len(np.unique(genes_all)))))
    genes_num_all = np.vectorize(gene_to_num_dict.get)(genes_all)


    # select tf columns and normalize
    data_all.drop(['tissue','exp','num_loop_counts','num_loops','num_atac_regions_pro','num_atac_regions_loop'],axis=1,inplace=True)
    selector = VarianceThreshold()
    data_all_varfilt = selector.fit_transform(data_all)
    data_all_varfilt_cols = data_all.columns[selector.get_support()]
    print(data_all.shape, data_all_varfilt.shape, len(data_all_varfilt_cols))
    scaler = MinMaxScaler()
    data_all_norm = scaler.fit_transform(data_all_varfilt)
    data_all_norm = pd.DataFrame(data_all_norm, columns = data_all_varfilt_cols)
    print('files_loaded', data_all.shape)

    return data_all_norm, labels_all, genes_all

def get_tissue_labels(labels_all,tissue):
    # select out section data containing on
    tissue_bool_labels = np.array(labels_all==tissue).flatten()
    tissue_bool_idx = np.where(tissue_bool_labels)[0]
    print('loaded tissue labels', len(tissue_bool_idx))
    return tissue_bool_labels, tissue_bool_idx


def get_genes(gene_filepath, genes_all):
    cell_type_genes = pd.read_csv(gene_filepath,header=None).loc[:,0]

    genes_bool = np.isin(genes_all, cell_type_genes)
    genes_bool_idx = np.where(genes_bool)[0]
    print('loaded genes', len(genes_bool),len(genes_bool_idx))

    return genes_bool,genes_bool_idx



def run_rf(data, labels, select_idx, save_dir, save_prefix, num_tries=1):
    """
    genes_bool_idx = select_idx
    """
    best_model = None
    best_roc = 0
    best_data = []

    for tr in range(num_tries):
        random_state = np.random.randint(0,1000)

        train_features, test_features, train_labels, test_labels = train_test_split(np.array(data)[select_idx,:],
                                                                                labels[select_idx],
                                                                                test_size = 0.25, random_state = random_state)



        print('Training Features Shape:', train_features.shape)
        print('Training Labels Shape:', train_labels.shape)
        print('Testing Features Shape:', test_features.shape)
        print('Testing Labels Shape:', test_labels.shape)

        print(Counter(labels[select_idx]))
        # Create the model with 100 trees
        model = RandomForestClassifier(n_estimators=100,#50,
                                       max_features = 'sqrt',
                                       n_jobs=-1, verbose = 0)
        # Fit on training data
        model.fit(train_features, train_labels)
        # Actual class predictions
        rf_predictions = model.predict(test_features)

        # evaluation
        acc = accuracy_score(test_labels, rf_predictions)
        print(acc)
    
        # Probabilities for each class
        rf_probs = model.predict_proba(test_features)[:, 1]
        # Calculate roc auc
        print(confusion_matrix(test_labels, rf_predictions))
        roc_value = roc_auc_score(test_labels, rf_probs)
        print(roc_value)
        # plot_roc_curve(model, test_features, test_labels)

        if roc_value>best_roc:
            best_model = model
            best_data = train_features, test_features, train_labels, test_labels
                
        
    if best_model is not None:
        fi = pd.DataFrame({'feature': list(data.columns),
                           'importance': best_model.feature_importances_}).\
                            sort_values('importance', ascending = False)
        output_file = os.path.join(save_dir, save_prefix+'_feature_importances.csv')

        fi.to_csv(output_file)
    return best_model, best_data

def run_irf(model, data,save_dir,save_prefix):
    print('running irf ...')
    train_features, test_features, train_labels, test_labels = data

    print('Get all Random Forest and Decision Tree Data ...')
    all_rf_tree_data = utils.get_rf_tree_data(rf=model, X_train=train_features, X_test=test_features, y_test=test_labels)
    print('Get the RIT data and produce RIT...')
    all_rit_tree_data = irf_utils.get_rit_tree_data(
                                                all_rf_tree_data=all_rf_tree_data,
                                                bin_class_type=1,
                                                M=100,
                                                max_depth=2,
                                                noisy_split=False,
                                                num_splits=2)


    print("Run the iRF function ....")
    try:

        all_rf_weights, all_K_iter_rf_data, \
        all_rf_bootstrap_output, all_rit_bootstrap_output, \
        stability_score = irf_utils.run_iRF(X_train=train_features,
                                        X_test=test_features,
                                        y_train=train_labels,
                                        y_test=test_labels,
                                        K=10,
                                        rf=RandomForestClassifierWithWeights(n_estimators=100),
                                        B=30,
                                        random_state_classifier=2018,
                                        propn_n_samples=.2,
                                        bin_class_type=1,
                                        M=25,
                                        max_depth=5,
                                        noisy_split=False,
                                        num_splits=2,
                                        n_estimators_bootstrap=5)

        print('writing stability scores ...')
        stability_score_names = {}
        for feat_idx, score in sorted(stability_score.items(),key=lambda x:x[1],reverse=True):
            feat_names_arr = [data_all_varfilt_cols_1[int(x)] for x in feat_idx.split('_')]
            feat_names = '::'.join(feat_names_arr)
            stability_score_names[feat_names] = score
        stability_df = pd.Series(stability_score_names)
        stability_df.to_csv(os.path.join(save_dir, save_prefix,'_stability_score.csv'))
    except:
        print('ERROR: could not run iRF no stability scores generated')
        


if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    # parser.add_argument('--tissue', dest='tissue',
 #                     help='tissue')
    # parser.add_argument('--save_dir', dest='save_dir', default = 'networks/representations_irf_cca'
 #                     help='tissue')

    save_dir = '/Users/mguo123/Google Drive/1_khavari/omics_project-LD/pan_omics/data/processed/fig4_modelling/irf'

    genes_file_paths = glob.glob('/Users/mguo123/Google Drive/1_khavari/omics_project-LD/pan_omics/data/processed/fig1/rna/*_genes.csv')
    cluster_names = [os.path.basename(f).split('_')[0] for f in genes_file_paths]

    if not os.path.exists(save_dir):
           os.makedirs(save_dir)

    data_all, labels_all, genes_all= load_data(tissues)


    for tissue_rf_label in tissues:

        for (cluster_name, gene_filepath) in zip(cluster_names, genes_file_paths):
            group = tissue_to_group[tissue_rf_label]
            if cluster_name == group: # only run on the group of interest
                save_prefix = tissue_rf_label + '_rf_' + cluster_name + '_genes'
                print(save_prefix)
                
                if os.path.exists(os.path.join(save_dir, save_prefix+'_feature_importances.csv')):                
                # if os.path.exists(os.path.join(save_dir, save_prefix,'_stability_score.csv')):
                    print('already ran', os.path.join(save_dir, save_prefix+'_feature_importances.csv'))
                else:
                    tissue_bool_labels, tissue_bool_idx = get_tissue_labels(labels_all, tissue_rf_label)
                    # gene_filepath = '/Users/mguo123/Google Drive/1_khavari/omics_project-LD/rnaseq/unique_gene_lists/'+tissue_geneset+'_genes.txt'
                    genes_bool, genes_bool_idx = get_genes(gene_filepath, genes_all)
                    try:
                        best_model, best_data = run_rf(data_all, tissue_bool_labels, genes_bool_idx, save_dir, save_prefix)
                    except:
                        print('failed rf ', save_prefix)
                        continue

                    try:
                        run_irf(best_model, best_data,save_dir,save_prefix)
                    except:
                        print('failed irf ', save_prefix)
                        continue                   
