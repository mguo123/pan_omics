"""
process_hichip.py

Margaret Guo

04/09/2020

takes output of Fithichip (*.bed files) and output *.loop_counts.bedpe files that can be used for diffloop analysis
(involves removing a column)

make loop_counts and anchor files

merge files

for anchor.csv make bed file out of it


notebooks
- fitchip_bedpe_preprocessing
- make_hichip_df
- make_hichip_df-withsnps-cancer-only
- annotate_bedpe_csv_rna
"""

import subprocess
import os, glob
import pandas as pd


def make_bedpe(input_filepath, output_filepath, split_delim = '.', extension=".bed", verbose=True):
    """
    takes output of Fithichip (*.bed files) and output *.loop_counts.bedpe files that can be used for diffloop analysis


    """
    print('processing HiChIP 1.. bed to loop_counts.bedpe')

    command_template = """tail -n +2 test.bed | awk '$7=("loop_"FNR FS $7)' | sed 's/\t/ /g' |sed '/GCContent/d' > test.loop_counts.bedpe"""
    cmd_arr = command_template.split()
    # replace cmd_arr[3] and cmd_arr[-1]

    if not os.path.exists(output_filepath):
        os.makedirs(output_filepath)

    # read file
    for subdir, dirs, files in os.walk(input_filepath):
        for filename in sorted(files):
            tissue = os.path.basename(subdir)
            input_filepath = subdir + os.sep + filename
            output_filepath_tissue = os.path.join(output_filepath, tissue)


            if filename.endswith(extension) :
                if not os.path.exists(output_filepath_tissue):
                    os.makedirs(output_filepath_tissue)
                sample = filename.split(split_delim)[0]
                output_file = os.path.join(output_filepath_tissue, filename.split(extension)[0] + ".loop_counts.bedpe")
                cmd_arr[3] = input_filepath
                cmd_arr[-1] = output_file
                if verbose:
                    print ('converting to bedpe... ', tissue, sample)
                    print(' '.join(cmd_arr))
                cmd = subprocess.run(str(' '.join(cmd_arr)),shell=True, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


###TODO comment
def hichip_to_df(bedpe_file):
    """
    bedpe to hichip_df
    """
    with open(bedpe_file, 'r') as f:
        hichip_dict = {}
        for idx, line in enumerate(f.readlines()):
            line_arr = line.strip().split()
    #         print(line_arr)
            if line_arr[0].startswith('chr'):
                prefix=''
            else:
                prefix='chr'
            start_A, end_A = line_arr[1], line_arr[2]
            start_B, end_B = line_arr[4], line_arr[5]
            if all([x.isdigit() for x in [start_A, end_A , start_B, end_B ]]):
                start_node=prefix+ line_arr[0] + '_' + start_A + '_' + end_A
                end_node = prefix+ line_arr[3] + '_' + start_B + '_' + end_B
                edge_attr = line_arr[7]
                hichip_dict[idx] = {'source': start_node, 'target':end_node, 'count':edge_attr}
        hichip_df = pd.DataFrame.from_dict(hichip_dict, orient='index')
#
    return hichip_df

###TODO comment
def hichip_to_anchor(hichip_df, sample):
    """
    hichip_df<DataFrame> from `hichip_to_df`
    """
    anchors_list = sorted(set(hichip_df.source).union(set(hichip_df.target)))
    anchors_df = pd.DataFrame({"anchors":anchors_list})
    anchors_df_split = anchors_df.anchors.str.split('_', expand=True)
    anchors_df_split.columns = ['chr', 'start', 'end']
    anchors_df = pd.concat([anchors_df, anchors_df_split], axis=1)
    anchors_df['sample']= sample
    return anchors_df


def make_csvs(input_filepath, output_filepath, split_delim = '.', extension=".loop_counts.bedpe", type_prefix='', verbose=True):
    """
    takes *.loop_counts.bedpe files  make .loop_counts.csv and .anchors.csv (separate samples),
    saves the sample separated csvs in input folder


    make new folder in interim/ called
        merged/
            loops/
            anchors/
            anchors_bed/
    Arguments:
        type_prefix <str> if different parameter versions of inputs are used

    Return: none
    """
    print('processing HiChIP 2.. loop_counts.bedpe to csv files,  merge by tissue, create bed anchor files')

    # set up file paths
    if not os.path.exists(output_filepath):
        os.makedirs(output_filepath)

    loops_dir = os.path.join(output_filepath, 'loops')
    anchors_dir = os.path.join(output_filepath, 'anchors')
    anchors_bed_dir = os.path.join(output_filepath, 'anchors_bed')
    if not os.path.exists(loops_dir):
        os.makedirs(loops_dir)
        os.makedirs(anchors_dir)
        os.makedirs(anchors_bed_dir)


    # read file
    for subdir, dirs, files in os.walk(input_filepath):
        tissue = os.path.basename(subdir)
        if verbose:
            print('merging and csv...', tissue)
        merged_loop_df = pd.DataFrame()
        merged_anchor_df = pd.DataFrame()

        # loop through tissue folder
        for file_idx, filename in enumerate(sorted(files)):

            input_filepath = os.path.join(subdir, filename )# bedpe file
            output_filepath_tissue = subdir

            if filename.endswith(extension) :

                sample = filename.split(split_delim)[0]
                output_loop_file = os.path.join(output_filepath_tissue, sample + ".loop_counts.csv")
                output_anchor_file = os.path.join(output_filepath_tissue, sample + ".anchors.csv")
                hichip_df = hichip_to_df(input_filepath)
                # loop_dict[sample] = hichip_df
                # save loop counts
                hichip_df.to_csv(output_loop_file)
                if verbose:
                    print ('converting to loop and anchor csvs... ', tissue, sample)
                    print('loops', hichip_df.shape)
                    print('wrote', output_loop_file)

                # save anchors
                #     hichip_df = hichip_to_df(bedpe_file)
                anchors_df = hichip_to_anchor(hichip_df, sample)
                anchors_df.to_csv(output_anchor_file)
                if verbose:
                    print('anchors', anchors_df.shape)
                    print('wrote', output_anchor_file)

                # merging loops and anchors per tissue
                if file_idx==0:
                    merged_loop_df = hichip_df
                    merged_anchor_df = anchors_df
                else:
                    merged_loop_df = pd.concat([merged_loop_df, hichip_df],ignore_index=True)
                    merged_loop_df = merged_loop_df.groupby(['source', 'target']).sum()
                    merged_loop_df.reset_index(inplace=True)

                    merged_anchor_df = pd.concat([merged_anchor_df, anchors_df],ignore_index=True)
                    merged_anchor_df.drop_duplicates(subset=['chr','start', 'end'],inplace=True)
                    merged_anchor_df.reset_index(drop=True,inplace=True)

        # save tissue files post merge + save bedfile
        if verbose:
            print('merged loops', tissue, merged_loop_df.shape)
            print('merged anchors', tissue, merged_anchor_df.shape)
        if merged_loop_df.shape[0]>0:
            merged_loop_df.to_csv(os.path.join(loops_dir, tissue+type_prefix+'.loops.csv'))
            merged_anchor_df.to_csv(os.path.join(anchors_dir, tissue+type_prefix+'.anchors.csv'))
            merged_anchor_df = merged_anchor_df[['chr','start','end', 'anchors']]
            merged_anchor_df.to_csv(os.path.join(anchors_bed_dir, tissue+type_prefix+'.anchors.bed'),sep='\t', index=False, header=False)






# def sort_bed_file(input_bedfile,ouput_bedfile):
