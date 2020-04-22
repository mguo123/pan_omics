"""
process_crms.py

04/20/2020

how to make crms from the annotated data

need to combine the data/interim/annon files

"""
import os, glob
import pandas as pd
from collections import Counter, defaultdict

files_dict = {
    'rna_file':'pass',
    'tf_annon_file' : 'data/external/HOCOMOCOv11_annotation.csv',
    'tissue':'pass',
    'atac_annon_file_promoter':'pass',
    'atac_annon_file_anchor':'pass',
    'footprinting_annon_file_promoter':'pass',
    'footprinting_annon_file_anchor':'pass',
    'anchor_promoter_file': 'pass',
    'loop_file':'pass',
    'output_dir' : 'data/processed/tissue_crms'
}

def create_crm_batch(rna_file, tf_annon_file, annon_file_path, loop_file_path,
                     output_dir,files_dict=files_dict, type='all',
                     THRES=1, verbose=False):
    """
    takes in a couple file paths, for each tissue
    makes sure all necessary files exist and then runs create_crm_per_tissue

    """

    if os.path.isfile(rna_file):
        files_dict['rna_file'] = rna_file
    else:
        raise ValueError("Need to run RNA processing to get tpm file")

    if os.path.isfile(tf_annon_file):
        files_dict['tf_annon_file'] = tf_annon_file
    else:
        raise ValueError("Need the TF HOCOMOCOv11_annotation file")

    files_dict['output_dir'] = output_dir
    if not os.path.exists(files_dict['output_dir']):
        os.makedirs(files_dict['output_dir'])
        os.makedirs(os.path.join(files_dict['output_dir'],'promoter_only_tissue'))
        os.makedirs(os.path.join(files_dict['output_dir'],'pro_anc_loop_tissue'))
        os.makedirs(os.path.join(files_dict['output_dir'],'pro_loop_tissue'))
        os.makedirs(os.path.join(files_dict['output_dir'],'combined_tissue'))

    if not os.path.isdir(annon_file_path):
        raise ValueError("Need to run annotation processing (bedtools) to get *_annon.bed files")
    if not os.path.isdir(loop_file_path):
        raise ValueError("Need to run hichip preprocessing to get *.loops.csv files")

    # setup global variables
    if type=='all':
        all_count_pro_overall = pd.DataFrame()
        all_count_sep_overall = pd.DataFrame()
        all_count_comb_overall = pd.DataFrame()
        all_count_sum_overall = pd.DataFrame()
    else:
        all_count_overall = pd.DataFrame()

    # loop through tissue
    for loop_file in sorted(glob.glob(os.path.join(loop_file_path, "*.loops.csv"))):
        tissue = os.path.basename(loop_file).split('.loops.csv')[0]
        # # #### DEBUG:
        # if tissue !='WM_shMITF_PLX':
        #     continue
        print('creating crms for tissue', tissue)
        files_dict['tissue'] = tissue
        files_dict['loop_file'] = loop_file

        atac_annon_file_promoter = os.path.join(annon_file_path, 'promoter_atac', 'promoter_'+tissue+'_annon.bed')
        if os.path.isfile(atac_annon_file_promoter):
            files_dict['atac_annon_file_promoter'] = atac_annon_file_promoter
        else:
            print("atac_annon_file_promoter file not found", atac_annon_file_promoter)
            continue

        atac_annon_file_anchor = os.path.join(annon_file_path, 'anchor_atac', tissue+'_annon.bed')
        if os.path.isfile(atac_annon_file_anchor):
            files_dict['atac_annon_file_anchor'] = atac_annon_file_anchor
        else:
            print("atac_annon_file_anchor file not found", atac_annon_file_anchor)
            continue

        footprinting_annon_file_promoter = os.path.join(annon_file_path, 'promoter_footprinting', 'promoter_'+tissue+'_annon.bed')
        if os.path.isfile(footprinting_annon_file_promoter):
            files_dict['footprinting_annon_file_promoter'] = footprinting_annon_file_promoter
        else:
            print("footprinting_annon_file_promoter file not found", footprinting_annon_file_promoter)
            continue

        footprinting_annon_file_anchor = os.path.join(annon_file_path, 'anchor_footprinting', tissue+'_annon.bed')
        if os.path.isfile(footprinting_annon_file_anchor):
            files_dict['footprinting_annon_file_anchor'] = footprinting_annon_file_anchor
        else:
            print("footprinting_annon_file_anchor file not found", footprinting_annon_file_anchor)
            continue

        anchor_promoter_file = os.path.join(annon_file_path, 'promoter_anchors', 'promoter_'+tissue+'_annon.bed')
        if os.path.isfile(anchor_promoter_file):
            files_dict['anchor_promoter_file'] = anchor_promoter_file
        else:
            print("anchor_promoter_file file not found", anchor_promoter_file)
            continue

        #
        print('file check passed')
        print(files_dict)
        # run crms
        crm_result = create_crm_per_tissue(files_dict,type=type, THRES=THRES, verbose=verbose )
        if crm_result is None:
            continue

        if type=='all':
            all_count_pro, all_count_sep, all_count_comb, all_count_sum = crm_result
            all_count_pro_overall = pd.concat([all_count_pro_overall, all_count_pro],sort=False).fillna(0)
            all_count_sep_overall = pd.concat([all_count_sep_overall, all_count_sep],sort=False).fillna(0)
            all_count_comb_overall = pd.concat([all_count_comb_overall, all_count_comb],sort=False).fillna(0)
            all_count_sum_overall = pd.concat([all_count_sum_overall, all_count_sum],sort=False).fillna(0)
        else:
            tissue_all_count = crm_result
            all_count_overall = pd.concat([all_count_overall, all_count_sum],sort=False).fillna(0)

    # end tissue loop
    # saving
    if type=='all':
        all_count_pro_overall.to_csv(os.path.join(output_dir,'all_count_pro_overall.csv'))
        all_count_sep_overall.to_csv(os.path.join(output_dir,'all_count_sep_overall.csv'))
        all_count_comb_overall.to_csv(os.path.join(output_dir,'all_count_comb_overall.csv'))
        all_count_sum_overall.to_csv(os.path.join(output_dir,'all_count_sum_overall.csv'))
    else:
        tissue_all_count = create_crm_per_tissue(files_dict,type=type, THRES=THRES, verbose=verbose )
        all_count_overall.to_csv(os.path.join(output_dir,'all_count_overall_{}.csv'.format(type)))


def create_crm_per_tissue(files_dict,type='all', THRES=1, verbose=False ):
    """
    create crm information per each tissue (see files_dict) example above
    ###TODO: add verbose stuff

    outputs  and saves a <DataFrame>
        index: TSSs with RNA values
        'tissue': tissue <str> (for later ease of combining)
        'exp': RNA TPM expression value
        'num_atac_regions<suffix>' : number of atac peaks per crm over designed region (pro, anc ,loop, or all)
        'num_loop_counts' and 'num_loops': number of counts for looped regions and number of looped regions respectively
        '<TF><suffix>': footprinting information depending on type, should be lots of values

    type:
        options:
            'all': does all the following:
            'promoter': only include promoter information in the result, will include "_pro" suffixes
            'anchor_loop_sep': include anchor and loop information separately, will include "_pro", "_anc" and "_loop" suffixes
            'anchor_loop_comb': include anchor and loop information summed, will include "_pro" and "_loop" suffixes
            'sum': agg all "_pro" "_anc" and "_loop", no suffixes
    """
    tissue = files_dict['tissue']
    rna_df = read_rna_file(files_dict['rna_file'])
    annon_dict =get_tf_dict(files_dict['tf_annon_file'])

    if tissue not in list(rna_df.columns.values):
        print ('tissue not found in rna matrix', tissue)
        return None

    results_promoter = crm_annotation(files_dict['atac_annon_file_promoter'],
                            files_dict['footprinting_annon_file_promoter'],
                            rna_df,annon_dict,files_dict['tissue'])

    results_loops = crm_annotation(files_dict['atac_annon_file_anchor'],
                            files_dict['footprinting_annon_file_anchor'],
                            rna_df,annon_dict,files_dict['tissue'],
                            anchor_promoter_file = files_dict['anchor_promoter_file'],
                            loop_file = files_dict['loop_file'])

    # create file DataFrame
    all_tfs = annon_dict.values()
    all_tf_all_tss = pd.DataFrame( index=rna_df.index)
    all_tf_all_tss['tissue'] = tissue
    all_tf_all_tss['exp'] = rna_df[tissue] # shape, (24811, 1)

    # get loop results
    loop_counts_df = results_loops['loop']

    # get footprinting results
    foot_count_pro = results_promoter['foot']
    foot_count_anchor, foot_count_loop = results_loops['foot']
    foot_count_sum = (foot_count_pro + foot_count_anchor + foot_count_loop).fillna(0)
    foot_count_comb = foot_count_anchor.add(foot_count_loop,fill_value=0)
    foot_count_pro.columns = [x+'_pro' for x in foot_count_pro
                              .columns.values]
    foot_count_anchor.columns = [x+'_anc' for x in foot_count_anchor.columns.values]
    foot_count_loop.columns = [x+'_loop' for x in foot_count_loop.columns.values]
    foot_count_comb.columns = [x+'_loop' for x in foot_count_comb.columns.values]

    # get atac results
    atac_count_pro = results_promoter['atac']
    atac_count_anchor, atac_count_loop = results_loops['atac']
    atac_count_sum = (atac_count_pro + atac_count_anchor + atac_count_loop).fillna(0)
    atac_count_pro.columns = ['num_atac_regions_pro']
    atac_count_comb = atac_count_anchor.add(atac_count_loop,fill_value=0)
    atac_count_anchor.columns = ['num_atac_regions_anc']
    atac_count_loop.columns = ['num_atac_regions_loop']
    atac_count_comb.columns = ['num_atac_regions_loop']


    # type = 'promoter': only include promoter information in the result, will include "_pro" suffixes
    all_count_pro = pd.concat([all_tf_all_tss, atac_count_pro, foot_count_pro],axis=1,sort=True).fillna(0)
    all_count_pro = all_count_pro[all_count_pro.index.isin(rna_df.index)]
    # type = 'anchor_loop_sep': include anchor and loop information separately, will include "_pro", "_anc" and "_loop" suffixes
    all_count_sep = pd.concat([all_tf_all_tss, loop_counts_df, atac_count_pro, atac_count_anchor, atac_count_loop,
                               foot_count_pro, foot_count_anchor,foot_count_loop],axis=1,sort=True).fillna(0)
    all_count_sep = all_count_sep[all_count_sep.index.isin(rna_df.index)]
    # type = 'anchor_loop_comb': include anchor and loop information summed, will include "_pro" and "_loop" suffixes
    all_count_comb = pd.concat([all_tf_all_tss, loop_counts_df, atac_count_pro, atac_count_comb,
                               foot_count_pro, foot_count_comb],axis=1,sort=True).fillna(0)
    all_count_comb = all_count_comb[all_count_comb.index.isin(rna_df.index)]
    # type = 'sum': agg all "_pro" "_anc" and "_loop", no suffixes
    all_count_sum = pd.concat([all_tf_all_tss, loop_counts_df, atac_count_sum,
                               foot_count_sum],axis=1,sort=True).fillna(0)
    all_count_sum = all_count_sum[all_count_sum.index.isin(rna_df.index)]

    # saving information
    if type=='all':
        all_count_pro.to_csv(os.path.join(files_dict['output_dir'],'promoter_only_tissue',tissue+'_crm.csv'))
        all_count_sep.to_csv(os.path.join(files_dict['output_dir'],'pro_anc_loop_tissue',tissue+'_crm.csv'))
        all_count_comb.to_csv(os.path.join(files_dict['output_dir'],'pro_loop_tissue',tissue+'_crm.csv'))
        all_count_sum.to_csv(os.path.join(files_dict['output_dir'],'combined_tissue',tissue+'_crm.csv'))
        return all_count_pro, all_count_sep, all_count_comb, all_count_sum
    elif type=='promoter':
        all_count_pro.to_csv(os.path.join(files_dict['output_dir'],'promoter_only_tissue',tissue+'_crm.csv'))
        return all_count_pro
    elif type=='anchor_loop_sep':
        all_count_sep.to_csv(os.path.join(files_dict['output_dir'],'pro_anc_loop_tissue',tissue+'_crm.csv'))
        return all_count_sep
    elif type=='anchor_loop_comb':
        all_count_comb.to_csv(os.path.join(files_dict['output_dir'],'pro_loop_tissue',tissue+'_crm.csv'))
        return all_count_comb
    elif type=='sum':
        all_count_sum.to_csv(os.path.join(files_dict['output_dir'],'combined_tissue',tissue+'_crm.csv'))
        return all_count_sum



def crm_annotation(atac_annon_file, footprinting_annon_file,rna_df,annon_dict,tissue,
                        anchor_promoter_file = None, loop_file = None,
                        THRES=1, verbose=False ):
    """
    get  annotation

    Arguments:
        footprinting_annon_file = <str> path, i.e. '../data/interim/annon/anchor_footprinting/GDSD3_annon.bed'
        rna_df <DataFrame> from `read_rna_file`
        annon_dict <Dict> motif<str> --> TF<str> from `get_tf_dict`
        anchor_promoter_file = <str> path, i.e. '../data/interim/annon/promoter_anchors/promoter_GDSD3_annon.bed'
        loop_file = <str> path, i.e. '../data/interim/merged/loops/GDSD3.loops.csv'
        THRES = <int> TPM threshold for considering something expressed
    Returns dictionary with results
        if only promoter region, return only foot_count  DataFrame
        o.w. return both foot_count_anchor,foot_count_loop DataFrame
    """

    foot_df = get_foot_df(footprinting_annon_file, annon_dict, rna_df, tissue, thres=THRES,verbose=verbose)
    atac_count = read_atac_file(atac_annon_file)
    if anchor_promoter_file is None: # them we are with the proximal/promoter region
        # promoter region (Case 1)
        # foot_df is already setup

        foot_count = get_foot_counts(foot_df)
        return {'atac':atac_count, 'foot':foot_count}

    else:
        # Case 2: getting anchor TFs and loops
        # foot_file = '../data/interim/annon/anchor_footprinting/GDSD3_annon.bed'
        foot_df = foot_df[['TF','p_loc']]

        # read in anchors
        # anchor_promoter_file = '../data/interim/annon/promoter_anchors/promoter_GDSD3_annon.bed'
        anchor_annon_df = pd.read_csv(anchor_promoter_file,sep='\t',header=None)
        anchor_annon_df.columns = ['chr_p','start_p','stop_p','TSS','chr_f','start_f','stop_f','anchor','overlap']
        anchor_annon_df_filt = anchor_annon_df.groupby('TSS').apply(lambda x: x.nlargest(1, "overlap")).reset_index(drop=True)
        anchor_annon_df_filt = anchor_annon_df_filt[['anchor','TSS']]
        # get footprints of the "annotated-with-TSS" anchor regions, then make the counts DataFrame
        foot_df_anchor = foot_df.merge(anchor_annon_df_filt,how='inner', left_on='p_loc',right_on='anchor')
        foot_count_anchor = get_foot_counts(foot_df_anchor)

        # read in loops, then "annotate" the source TSS regions the loop anchors
        loop_df_bi = make_loop_df(loop_file)
        loop_df_bi = loop_df_bi.merge(anchor_annon_df_filt,how='inner',left_on='source',right_on='anchor' )
        # get footprints of the "annotated-with-TSS" looped regions, then make the counts DataFrame
        foot_df_loop = foot_df.merge(loop_df_bi,how='inner', left_on='p_loc',right_on='target')
        foot_count_loop = get_foot_counts(foot_df_loop)


        ### for atac_num_regions, anchor and loop
        atac_count_anchor = atac_count.merge(anchor_annon_df_filt,how='inner', left_index=True,right_on='anchor')
        atac_count_anchor = pd.DataFrame(atac_count_anchor.groupby('TSS')['num_atac_regions'].sum())
        atac_count_anchor.index.name = None

        atac_count_loop = atac_count.merge(loop_df_bi,how='inner',left_index=True,right_on='target')
        atac_count_loop = pd.DataFrame(atac_count_loop.groupby('TSS')['num_atac_regions'].sum())
        atac_count_loop.index.name = None

        # loop_counts
        loop_counts_df = loop_df_bi.groupby('TSS').agg({"count":"sum","target":"count"})
        loop_counts_df.columns = ['num_loop_counts','num_loops']
        loop_counts_df.index.name = None

        return {'atac':[atac_count_anchor, atac_count_loop], 'foot':[foot_count_anchor,foot_count_loop],'loop':loop_counts_df }


    # # foot_file = '../data/interim/annon/anchor_footprinting/GDSD3_annon.bed'
    # foot_df = get_foot_df(foot_file, annon_dict, rna_df, thres=1,verbose=False)
    #
    # pass



######## HELPER FUNCTIONS #########
def get_exp(rna_df, sym, tissue,verbose=True):
    try:
        return rna_df.loc[sym,tissue]
    except KeyError:
        if verbose:
            print(sym, 'not found')
        return 0

def read_rna_file(rna_file):
    rna_df = pd.read_csv(rna_file, index_col=0,header=0)
    rna_df.index = [x.upper() for x in rna_df.index.values]
    return rna_df

def read_atac_file(atac_file):
    """
    read in to DataFrame
    atac_file = '../data/interim/annon/promoter_atac/promoter_GDSD3_annon.bed'
    and get # of atac regions per tss

    output: atac_count <DataFrame> with `TSS` and index and  `num_atac_regions` as columns
    """
    filename = os.path.basename(atac_file)
    tissue = filename.split('_annon.bed')[0].split("_")[-1]
    atac_df = pd.read_csv(atac_file,sep='\t', header=None)
    atac_df.columns = ['chr_p','start_p','stop_p','TSS','chr_a','start_a','stop_a']
    atac_df['p_loc'] = atac_df.chr_p + '_' + atac_df.start_p.map(str) + '_' + atac_df.stop_p.map(str)
    atac_df['a_loc'] = atac_df.chr_a + '_' + atac_df.start_a.map(str) + '_' + atac_df.stop_a.map(str)
    atac_count = pd.DataFrame(atac_df.groupby('TSS')['p_loc'].count())
    atac_count.columns=['num_atac_regions']
    atac_count.index.name = None

    return atac_count


def get_foot_df(foot_file, annon_dict, rna_df, tissue, thres=1,verbose=False):
    """
    from footprinting annotation file, annotate with motif/TF (rna and transcription factor info)

    Arguments:
        foot_file  = <str> to filepath i.e. '../data/interim/annon/promoter_footprinting/promoter_GDSD0_annon.bed'
        annon_dict = dict from `get_tf_dict`
        rna_df = DataFrame from `read_rna_file`, filter by expression motif --> TF
        tissue = <str>
        thres = <int> for considering what is expressed

    Return:
        foot_df_filt = <DataFrame> with columns, ['TSS','TF', 'p_loc']
            indicating the TSS (gene) the TF is annotated for and the region annotated in`
                p_loc = <str> with <chr>_<start>_<stop>, can correspond to anchor region
    """
    # read in footprinting file
    filename = os.path.basename(foot_file)
    foot_df = pd.read_csv(foot_file,sep='\t',
                header=None)
    foot_df.columns = ['chr_p','start_p','stop_p','TSS','chr_f','start_f','stop_f','motif_id','score','-']
    foot_df['p_loc'] = foot_df.chr_p + '_' + foot_df.start_p.map(str) + '_' + foot_df.stop_p.map(str)
    foot_df['f_loc'] = foot_df.chr_f + '_' + foot_df.start_f.map(str) + '_' + foot_df.stop_f.map(str)
    foot_df[['motif_abbr','motif_info']]=foot_df['motif_id'].str.split("_",expand=True)
    foot_df['animal']=foot_df['motif_info'].str.split(".",expand=True).loc[:,0]
    # human tfs only
    foot_df = foot_df[foot_df.animal=='HUMAN']

    # mapping motif --> TF with the annotation dictionary
    foot_df['TF'] = foot_df.motif_abbr.map(annon_dict)
    foot_df = foot_df[~foot_df.TF.isna()]

    # get expression from rna_df
    foot_df['TF_exp']= foot_df.TF.apply(lambda x: get_exp(rna_df, x, tissue,verbose=verbose))
    foot_df_filt = foot_df[foot_df.TF_exp>thres]
    foot_df_filt = foot_df_filt[['TSS','TF', 'p_loc']]


    return foot_df_filt


def get_foot_counts(foot_df):
    """
    from a foot_df from `get_foot_df`
    get `foot_df_count_pivot` that is a mxn DataFrame
        with m = # of genes annotated in foot file, and
            n = # of tfs annotated from `annon_dict`
            filtered through expression,

    Arguments:
        foot_df: <DataFrame> with at least required columns TSS, TF (all other columns are ignored)


    """
    # pivoting table to create counts DataFrame
    foot_df_filt_count = pd.DataFrame({'count':foot_df.groupby('TSS')['TF'].value_counts()}).reset_index()
    foot_df_filt_count_pivot = foot_df_filt_count.pivot(index='TSS',columns='TF',values='count').fillna(0)
    return foot_df_filt_count_pivot

def get_tf_dict(annon_file):
    # return dictionary mapping motif to tf (for HOCOMOCOv11_annotation)
    # annon_file = '../data/external/HOCOMOCOv11_annotation.csv'
    annon_df = pd.read_csv(annon_file)
    annon_df[['motif_abbr','motif_info']]=annon_df['id'].str.split("_",expand=True)
    annon_dict = pd.Series(annon_df.tf.values, index=annon_df.motif_abbr.values).to_dict()
    return annon_dict

def make_loop_df(loop_file):
    """
    from hichip_df make hichip df that's bidirectionary target <-> source


    Arguments
        loop_file : <str> like '../data/interim/merged/loops/GDSD3.loops.csv'

    return: loop_df_bi <DataFrame> with colummns ['target','source','count']

    """
    # loop_file = '../data/interim/merged/loops/GDSD3.loops.csv'
    loop_df = pd.read_csv(loop_file,index_col=0)
    loop_df_rev = loop_df.copy()
    loop_df_rev.columns = ['target','source','count']
    loop_df_bi = pd.concat([loop_df, loop_df_rev],sort=False)
    return loop_df_bi
