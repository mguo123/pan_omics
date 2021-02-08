"""
preprocess_atac.py


# Remove tab from end of footprinting files and then bedtools sort and merge each tissue folder store all results in merged folder
cd /Users/mguo123/Google Drive/1_khavari/omics_project-LD/atac_footprinting/match_hoco
#### sed 's/[[:blank:]]*$//'`
for f in `ls *bed`; do
echo $f; tissue=${f%%_*};echo $tissue;
f2=${tissue}_footprinting_hoco.bed;
f3=${tissue}_footprinting_hoco_sort.bed;
sed 's/[[:blank:]]*$//' $f > $f2;
bedtools sort -i $f2 > $f3;
bedtools merge ...
done

"""
import os,sys, glob
import subprocess
import pybedtools
import pandas as pd

def preprocess_footprinting(input_filepath, output_filepath,split_delim = '.', type_prefix='', extension=".bed", verbose=True):
    """
    # Remove tab from end of footprinting files, and saving in interim files, same directory structure as raw

    """
    print('processing Footprinting .. ')

    # cmd_template1 = """sed 's/[[:blank:]]*$//' $f > $f2"""
    # cmd_arr1 = cmd_template1.split()
    # replace cmd_arr1[2] and cmd_arr1[-1]

    if not os.path.exists(output_filepath):
        os.makedirs(output_filepath)

    # read file
    for subdir, dirs, files in os.walk(input_filepath):
        tissue = os.path.basename(subdir)
        bedtool_obj = None

        for filename in sorted(files):
            input_filepath = os.path.join(subdir, filename )#  file full path

            output_filepath_tissue = os.path.join(output_filepath, tissue)


            if filename.endswith(extension) :
                # if not os.path.exists(output_filepath_tissue):
                #     os.makedirs(output_filepath_tissue)
                sample = filename.split(split_delim)[0]
                bedtool_df = pd.read_csv(input_filepath, sep='\t', header=None)

                if bedtool_obj is None:

                    bedtool_obj = pybedtools.BedTool.from_dataframe(bedtool_df.loc[:,:6])
                else:
                    bedtool_obj = bedtool_obj.cat(pybedtools.BedTool.from_dataframe(bedtool_df.loc[:,:6]))

                # output_file1 = os.path.join(output_filepath_tissue, filename.split(extension)[0] + "_pre.bed")
                # cmd_arr1[2] = input_filepath
                # cmd_arr1[-1] = output_file1
                # cmd1 = subprocess.run(str(' '.join(cmd_arr1)),shell=True, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                # # sorting
                # output_file2 = os.path.join(output_filepath_tissue, filename.split(extension)[0] + "_sort.bed")
                # bedtool_obj = pybedtools.BedTool(output_file1)
                # bedtool_obj.sort().saveas(output_file2)
                if verbose:
                    print ('preprocessing footprinting... ', tissue, sample)
                    # print(' '.join(cmd_arr1))
        # save tissue files post merge + save bedfile

        if bedtool_obj is not None:
            if verbose:
                print('num peaks', tissue, bedtool_obj.count())
            bedtool_obj.sort().saveas(os.path.join(output_filepath, tissue+type_prefix+'_merged.bed'))



def merge_samples(input_filepath, output_filepath, split_delim = '.bed', extension=".bed", type_prefix='', verbose=True):
    """
    takes all sample bed files within a tissue file with a certain `extension` merges them


    make new folder <output_filepath> (usually in in interim/merged/*)

    Arguments:
        type_prefix <str> if different parameter versions of inputs are used

    Return: none
    """
    print('processing files 2.. loop_counts.bedpe to csv files,  merge by tissue, create bed anchor files')

    # set up file paths
    if not os.path.exists(output_filepath):
        os.makedirs(output_filepath)


    # read file
    for subdir, dirs, files in os.walk(input_filepath):
        tissue = os.path.basename(subdir)
        bedtool_obj = None
        if verbose:
            print('merging ...', tissue)
        # loop through tissue folder
        for file_idx, filename in enumerate(files):

            input_filepath = os.path.join(subdir, filename )#  file full path
            # output_filepath_tissue = subdir

            if filename.endswith(extension) :
                sample = filename.split(split_delim)[0]
                # merging loops and anchors per tissue
                if bedtool_obj is None:
                    bedtool_obj = pybedtools.BedTool(input_filepath)
                else:
                    bedtool_obj = bedtool_obj.cat(pybedtools.BedTool(input_filepath))

        # save tissue files post merge + save bedfile

        if bedtool_obj is not None:
            if verbose:
                print('num peaks', tissue, bedtool_obj.count())
            bedtool_obj.sort().saveas(os.path.join(output_filepath, tissue+type_prefix+'_merged.bed'))
