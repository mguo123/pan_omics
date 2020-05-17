"""
Margaret Guo
05/15/2020

trim PWM by IC 


"""

import numpy as np
import os,glob
import pandas as pd
import argparse



class PCM(object):
	def __init__(self, pcm_file,thres_IC = 0.4,pseudo=1,verbose=True):
		self.pcm_file = pcm_file
		self.thres_IC = thres_IC
		self.pseudo = pseudo
		self.verbose = verbose

		# get number of bases and read pcm file (# bases in alphabet DNA=4, protein ~20, doesn't change per PCM)
		self.pcm = self.read_pcm_file(self.pcm_file)
		self.num_bases = self.pcm.shape[0]
		self.update_pcm_params()

	def read_pcm_file(self,pcm_file):
		"""
		read in pcm which is an BxN matrix with B = # bases in alphabet and N = length of sequence motif
		"""
		pcm =  np.array(pd.read_csv(pcm_file, sep=' ',header=None))#,dtype='float64')
		pcm += self.pseudo
		return pcm

	def update_pcm_params(self):
		"""
		update, PPM, IC, length sequence

		called in trim_pcm()
		"""
		self.len_motif = self.pcm.shape[1]
		self.ppm =self.calculate_ppm()
		self.ic = self.calculate_ic()

	def calculate_ppm(self):
		"""
		calculate position probability matrix
		PPM[b,i]
		"""
		if not type(self.pcm) == np.ndarray:
			raise ValueError('no position count matrix')
		return self.pcm/np.sum(self.pcm,axis=0)

	def calculate_ic(self):
		"""
		calculate informtion content for each position i
		I_i = 2 + sum_(over each base b) ppm[b,i] * log2 ppm[b,i]
		return IC 
		"""
		if not type(self.ppm) == np.ndarray:
			raise ValueError('no position probability matrix')
		I = 2.0 + np.sum(self.ppm * np.log2(self.ppm),axis=0)
		return I

	def trim_pcm(self):
		"""

		trim PCM for end positions > self.thres_IC, 

		update self.len_motif, self.pcm, self.ppm, self.ic
		"""
		good_idx = np.where(self.ic > self.thres_IC)[0]
		start_idx = good_idx[0]
		end_idx = good_idx[-1]
		new_pcm = self.pcm[:,start_idx:(end_idx+1)]
		if self.verbose:
			print('trimming pcm', self.len_motif, new_pcm.shape[1])
		self.pcm = new_pcm
		self.update_pcm_params()


	def write_pcm_file(self, file_out,sep=' '):
		pd.DataFrame(self.pcm).to_csv(file_out, sep=' ',header=None,index=None)
		# with open(file_out, 'w') as f:
		# 	for b in range(self.num_bases):
		# 		f.write(' '.join(pcm[b,:])+'\n')






if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Parsing position count matrices (PCMs)')
	parser.add_argument('-i', dest='in_dir', type=str, 
	                   help='input directory of pcms, with *pwm extensions')
	parser.add_argument('-o', dest='out_dir', type=str, 
	                   help='output directory of pcms, with *pwm extensions')
	parser.add_argument('-v', dest='verbose', type=bool, default=False,
	                   help='output directory of pcms, with *pwm extensions')	
	args = parser.parse_args()

	# in_dir = '../../data/raw/motifs/' # these are both testing dir
	# out_dir = '../../data/raw/motifs_trim'
	if not os.path.exists(args.out_dir):
		os.makedirs(args.out_dir)

	for pwm_file in glob.glob(os.path.join(args.in_dir, '*.pwm')):
		print(pwm_file)
		pwm_name = os.path.basename(pwm_file).split('*.pwm')[0]

		pcm = PCM(pwm_file,verbose=args.verbose)
		pcm.trim_pcm()
		pcm.write_pcm_file(os.path.join(args.out_dir, pwm_name + '.trim.pwm'))
