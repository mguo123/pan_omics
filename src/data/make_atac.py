"""
process_atac.py

Margaret Guo

04/15/2020


footprinting (.bed) --> csv)



notebooks
- make_hichip_df-withsnps-cancer-only
"""
#### FIX FOR HOCO
def footprinting_to_df(footprinting_file):
    footprinting_df = pd.read_table(footprinting_file,header=None)
    footprinting_df.columns = ['chr', 'start', 'end', 'motif', 'score', 'strand', 'other']
    temp = footprinting_df.motif.str.split('(') #i.e. remove (var.2)
    temp = temp.apply(lambda x: x[0])
    temp = temp.str.split('.')
    footprinting_df['motif_abbr'] = temp.apply(lambda x: x[-1])
    return footprinting_df

def atac_to_df(atac_narrowPeak_file):
    atac_df = pd.read_table(atac_narrowPeak_file,header=None)
    atac_df.columns = ['chr', 'start', 'end']
    atac_df = atac_df.groupby(atac_df.columns.tolist()).size().reset_index().rename(columns={0:'count'})
    atac_df['region'] = atac_df.apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
    return atac_df
