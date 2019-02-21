import argparse, random, subprocess, glob
import numpy as np
import pandas as pd
import re

def pi(SFS,per_site = True):
	if sum(SFS) ==0: 
		return -99
	N = len(SFS)-1
	binom = (N * (N -1))/2
	pi = sum([(1.0*i*(N-i)*(SFS[i]))/(binom) for i in xrange(N) if i != 0])
	if per_site == True:
		return pi/sum(SFS)
	else:
		return pi

def parseLabel(label):
	d = label.split(':')[0]
	if d == 'u':
		mult = -1.
	elif d == 'd':
		mult = 1.
	dist = sum(map(float,label.split(':')[1].split('-')))/2
	return dist * mult
		
def main():
	parser = argparse.ArgumentParser(description="Gets the SFS files from a directory and collates the site frequency spectra from a directory containing many entries with the suffix *sfs.txt.gz")
	parser.add_argument("--input","-i", 
			required = True,
			dest = "input",
			type =str, 
			help = "name the directory containing the site frequnecy spectra files")

	args = parser.parse_args()
	data = {}
	for i in glob.glob(args.input +'/*sfs.txt.gz'):
		for k in pd.read_csv(i, compression = 'gzip').values:
			
			try:
				data[k[0]].append(k[2:])
			except KeyError:
				data[k[0]] = [k[2:]]
	for key in sorted(data.keys()):
		data[key] = np.array(data[key])
		full_sfs = np.sum(data[key], axis=0)
		print parseLabel(key), full_sfs.sum(),pi(full_sfs)
		


main()
