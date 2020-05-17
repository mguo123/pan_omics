"""
process_bedtools.py

Margaret Guo 04/20/2020

all the subprocess commands with bedtools


"""
import os,glob,sys
import subprocess
import pybedtools

def sort_bed_dir(input_beddir,output_beddir):
    """
    sort all bed files in directory
    cmd line call:
        bedtools sort -i ${FILENAME}.bed > ${FILENAME}_sort.bed
    """
    print("processing .. sorting  bed files ",input_beddir)
    if not os.path.exists(input_beddir):
        raise("input directory does not exist")
    if not os.path.exists(output_beddir):
        os.makedirs(output_beddir)

    cmd = ["bedtools", "sort", "-i", "${FILENAME}.bed"]
    for file in glob.glob(os.path.join(input_beddir, '*bed')):
        file_name = os.path.basename(file).split('.')[0]
        output_filename = os.path.join(output_beddir, file_name+"_sort.bed")
        bedtool_obj = pybedtools.BedTool(file)
        bedtool_obj.sort().saveas(output_filename)
        # cmd[-1] = file
        # with open(output_filename, "w") as file:
        #     subprocess.run(cmd, check=True, stdout=file)


def annotate_batch(input, db_dir, output_dir, extensions=['.bed','.bed'],f=1E-9,wo=False, verbose=True):
    """
    calls `annotate` for directory

    if <input>=`input_prefix`<extensions[0]> is a file, the run input on db and store output bedfile as:
        db = bed file prefix in db_dir: <db><extensions[1]>
         <output_dir>/`input_prefix`_<db>_annon.bed

    elif <input> is a directory with many bed files inside,
    then run annotation matched by name of input and db
        run annotate with below arguments
            input_bedfile = <input>/`tissue`<extensions[0]>
            db_bedfile = <db_dir>/`tissue`<extensions[1]>
            output_bedfile = <output_dir>/`tissue`_annon.bed

    return: None
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if os.path.isdir(input):
        for input_filepath in sorted(glob.glob(os.path.join(input, "*"+extensions[0]))):
            input_filename = os.path.basename(input_filepath)
            tissue = input_filename.split(extensions[0])[0]

            db_filepath = os.path.join(db_dir, tissue + extensions[1])
            output_filepath = os.path.join(output_dir, tissue+'_annon.bed')
            if os.path.isfile(db_filepath):
                if verbose:
                    print('annotating in ',output_dir, '...', tissue)
                annotate(input_filepath, db_filepath,output_filepath,f=f,wo=wo, verbose=verbose)
            else:
                print('WARNING:',tissue, 'not found in ', db_dir)

    elif os.path.isfile(input):
        input_filename = os.path.basename(input)
        input_prefix = input_filename.split(extensions[0])[0]
        for db_filepath in sorted(glob.glob(os.path.join(db_dir, "*"+extensions[1]))):
            db = os.path.basename(db_filepath).split(extensions[1])[0]
            output_filepath = os.path.join(output_dir, input_prefix+"_"+db+'_annon.bed')
            # print('input',input)
            # print('db_filepath',db_filepath)
            # print('output_filepath',output_filepath)
            if verbose:
                print('annotating in ',output_dir, '...', db)
            annotate(input, db_filepath,output_filepath,f=f, wo=wo,verbose=verbose)



def annotate(input_bedfile, db_bedfile,output_bedfile,f=1E-9, wo=False, verbose=True):
    """
    calls
        bedtools intersect -wa -wb -a $input_bedfile -b $db_bedfile -sorted -names footprinting -f $f > $output_bedfile;


    -f	Minimum overlap required as a fraction of A. Default is 1E-9 (i.e. 1bp)
        for anchors 5kb overlap, overlap should be 2.5E-6
        for promoter 2kb region, overlap should be 1E-6

    -wo ????

    return: None
    """
    input_bed = pybedtools.BedTool(input_bedfile)
    db_bed = pybedtools.BedTool(db_bedfile)
    if wo:
        output_bed = input_bed.intersect(db_bed, wo=True, sorted=True, names='db',f=f)
    else:
        output_bed = input_bed.intersect(db_bed, wa=True, wb=True, sorted=True, names='db',f=f)
    output_bed.saveas(output_bedfile)
    if verbose:
        print('annotated to', output_bedfile)
