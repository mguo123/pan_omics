"""
build_random_loops.py


Margaret Guo
June 22, 2020
"""

import pandas as pd
import os, glob
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pybedtools
from collections import defaultdict, Counter

###TODO: make super class
###TODO: clean up code
class Hichip_random(object):
    def __init__(self, merged_dir, tissue, num_anchors, num_loops, anchors_all, save_bool=True, allow_self_loop=False):
        self.tissue = tissue
        self.merged_dir = merged_dir
        self.save_bool = save_bool
        self.allow_self_loop = allow_self_loop
        if self.save_bool:
            self.random_anchor_dir = os.path.join(merged_dir, 'anchors_random')
            if not os.path.exists(self.random_anchor_dir):
                os.makedirs(self.random_anchor_dir)
            self.save_anchor_path = os.path.join(self.random_anchor_dir, self.tissue+'.anchors.csv')

            self.random_anchor_bed_dir = os.path.join(merged_dir, 'anchors_random_bed_sort')
            if not os.path.exists(self.random_anchor_bed_dir):
                os.makedirs(self.random_anchor_bed_dir)
            self.save_anchor_bed_path = os.path.join(self.random_anchor_bed_dir, self.tissue +'.sort.bed')

            self.random_loop_dir = os.path.join(merged_dir, 'loops_random')
            if not os.path.exists(self.random_loop_dir):
                os.makedirs(self.random_loop_dir)
            self.save_loop_path = os.path.join(self.random_loop_dir, self.tissue+'.loops.csv')


        # self.tissue_anchor_file = os.path.join(merged_dir, 'anchors', tissue+'.anchors.csv')
        # self.tissue_anchor_df = pd.read_csv(self.tissue_anchor_file,index_col=0)
        # self.tissue_loop_file = os.path.join(merged_dir  , 'loops', tissue +'.loops.csv' )
        # self.tissue_loop_df = pd.read_csv(self.tissue_loop_file,index_col=0)

        self.num_anchors = num_anchors#self.tissue_anchor_df.shape[0]
        self.num_loops = num_loops#self.tissue_loop_df.shape[0]
        self.anchors_all = anchors_all
        # self.num_loops_counts = sum(self.tissue_loop_df['count'])


        self.anchors, self.tissue_anchor_df = self.get_random_anchors()#self.tissue_anchor_df.anchors.values
        # self.prob_anchors_counter = Counter(self.tissue_loop_df.source)+Counter(self.tissue_loop_df.target)
        # self.prob_anchor_df = pd.DataFrame.from_dict({'anchor':list(self.prob_anchors_counter.keys()), 'count':list(self.prob_anchors_counter.values())})
        # self.prob_anchor_df['prob'] = self.prob_anchor_df['count']/sum(self.prob_anchor_df['count'])

        # self.chr_anchor_prob_dict = {} # dict: key chrom, value: dataframe of anchor,count, and probability (within chrom)
        self.chr_anchor_dict = {} # dict: key chrom, value: dataframe of anchor,count, and probability (within chrom)
        for chrom in self.tissue_anchor_df.chr.unique():
            anchor_df_chr = self.tissue_anchor_df[self.tissue_anchor_df.chr==chrom]
            self.chr_anchor_dict[chrom] = anchor_df_chr.anchors.values





    def get_random_anchors(self):
        random_anchors = np.random.choice(list(self.anchors_all), self.num_anchors,replace=False)
        random_anchors_df = pd.DataFrame(columns = ['anchors','chr','start','end'])
        random_anchors_df['anchors'] = random_anchors
        random_anchors_df[['chr','start','end']] = random_anchors_df['anchors'].str.split('_',expand=True)
        if self.save_bool:
            random_anchors_df.to_csv(self.save_anchor_path)
            print('saved anchors', self.save_anchor_path)
            random_anchors_df=random_anchors_df[['chr','start','end', 'anchors']]
            bedtool_obj = pybedtools.BedTool.from_dataframe(random_anchors_df)
            bedtool_obj.sort().saveas(self.save_anchor_bed_path)
            print('saved anchors bed', self.save_anchor_bed_path)

        return random_anchors, random_anchors_df


    def get_random_loops(self):
        source_arr = np.random.choice(self.anchors,size=self.num_loops,replace=True)#, p=self.prob_anchor_df['prob'])
        target_arr = []
        for anchor in source_arr:
            chrom = anchor.split('_')[0]
            target_choice = np.random.choice(self.chr_anchor_dict[chrom])#, p=self.chr_anchor_prob_dict[chrom].prob)
            if not self.allow_self_loop:
                while target_choice==anchor:
                    target_choice = np.random.choice(self.chr_anchor_dict[chrom])#np.random.choice(self.chr_anchor_prob_dict[chrom].anchor, p=self.chr_anchor_prob_dict[chrom].prob)

            target_arr.append(target_choice)
    #
        random_loop_df = pd.DataFrame()
        random_loop_df['source'] = source_arr
        random_loop_df['target'] = target_arr
        random_loop_df['count'] = 1#self.tissue_loop_df['count']

        if self.save_bool:
            random_loop_df.to_csv(self.save_loop_path)
            print('saved loop', self.save_loop_path)
        return random_loop_df
    #


        #
        # for loop_idx in range(self.num_loops):
            # np.random.choice(astro.anchors)


def make_random(merged_dir,sel_tissues = None):
    anchor_path=os.path.join(merged_dir,'anchors')
    # 2. get all anchors, counts for loops and anchors for each tissue
    anchors_all = set()
    anchor_counts = {}
    loop_counts = {}
    for anchor_file in glob.glob(os.path.join(anchor_path, '*csv')):
        tissue = os.path.basename(anchor_file).split('.')[0]
        if sel_tissues is not None:
            if tissue not in sel_tissues:
                continue
        print(tissue)
        loop_file = os.path.join(merged_dir, 'loops', tissue+'.loops.csv')
        tissue_anchor_df = pd.read_csv(anchor_file,index_col=0)
        tissue_loop_df = pd.read_csv(loop_file,index_col=0)

        num_anchors = tissue_anchor_df.shape[0]
        num_loops = tissue_loop_df.shape[0]
    #     num_loops_counts = sum(tissue_loop_df['count'])
        anchor_counts[tissue] = num_anchors
        loop_counts[tissue] = num_loops

        anchors = set(tissue_anchor_df.anchors.values)
        anchors_all = anchors_all.union(anchors)






    # anchor_files = glob.glob(os.path.join(merged_dir, 'anchors', '*.anchors.csv'))
    # all_anchors = set()
    # for anchor_file in anchor_files:
    #     tissue = os.path.basename(anchor_file).split('.')[0]
    for tissue in anchor_counts.keys():
        print(tissue)
        # tissue_anchor_file = os.path.join(merged_dir, 'anchors', tissue+'.anchors.csv')
        # tissue_anchor_df = pd.read_csv(self.tissue_anchor_file,index_col=0)
        # tissue_loop_file = os.path.join(merged_dir  , 'loops', tissue +'.loops.csv' )
        # tissue_loop_df = pd.read_csv(self.tissue_loop_file,index_col=0)
        obj = Hichip_random(merged_dir, tissue,num_anchors=anchor_counts[tissue], num_loops=loop_counts[tissue],
                            anchors_all = anchors_all, save_bool=True)
        obj.get_random_loops() # reformat object for better initiation
        # print('saved random loops..', obj.save_path)

"""
print(os.getcwd())
merged_dir = '/Users/mguo123/Google Drive/1_khavari/omics_project-LD/pan_omics/data/interim/merged'

os.listdir()
tissue = 'Astrocytes'
astro = Hichip(merged_dir, tissue)
astro.get_random_anchors()

"""



"""
anchors
0,chr10_100000000_100005000,chr10,100000000,100005000,Astro_B2
1,chr10_100005000_100010000,chr10,100005000,100010000,Astro_B2
2,chr10_100015000_100020000,chr10,100015000,100020000,Astro_B2
3,chr10_100020000_100025000,chr10,100020000,100025000,Astro_B2
4,chr10_100025000_100030000,chr10,100025000,100030000,Astro_B2
5,chr10_100030000_100035000,chr10,100030000,100035000,Astro_B2


,source,target,count
0,chr10_100020000_100025000,chr10_100055000_100060000,13
1,chr10_100020000_100025000,chr10_100065000_100070000,13
2,chr10_100025000_100030000,chr10_100055000_100060000,14
3,chr10_100025000_100030000,chr10_100065000_100070000,11
4,chr10_101150000_101155000,chr10_101180000_101185000,14
5,chr10_101150000_101155000,chr10_101190000_101195000,1529

"""
