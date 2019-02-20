## A script that is a Python rewrite of Dan's get_bed_flanks.R script.

import argparse, math
import numpy as np
import pandas as pd

def readBed(input_file):
	x = pd.read_csv(input_file,  names = ["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"], sep = '\t')
	if pd.isnull(x.strand[1]):
		x['strand'] = '+'
	return x
	
	
def getLwrWin(x, offset, length):
	x['chromEnd'] = np.ceil(np.maximum(x['chromStart'] - offset, x['lwr_lim']))
	x['chromStart'] = np.ceil(np.maximum(x['chromEnd'] - length, x['lwr_lim']))

	ind = x['chromEnd'] == x['chromStart']


	x.loc[ind,'chromEnd'] = np.nan
	x.loc[ind,'chromStart'] = np.nan
	x['direction'] = np.nan
	x.loc[x['strand'] == '-', 'direction'] = 'd'
	x.loc[x['strand'] == '+', 'direction'] = 'u'
	x['length'] = x['chromEnd'] - x['chromStart']
	return x
	
def getUprWin(x, offset, length):
	x['chromStart'] = np.floor(np.minimum(x['chromEnd'] + offset, x['upr_lim']))
	x['chromEnd'] = np.floor(np.minimum(x['chromStart'] + length, x['upr_lim']))

	ind = x['chromEnd'] == x['chromStart']
	x.loc[ind,'chromEnd'] = np.nan
	x.loc[ind,'chromStart'] = np.nan
	x['direction'] = np.nan
	x.loc[x['strand'] == '-', 'direction'] = 'u'
	x.loc[x['strand'] == '+', 'direction'] = 'd'
	x['length'] = x['chromEnd'] - x['chromStart']
	return x
	
	
	
def main():
	parser = argparse.ArgumentParser(description="PYTHON 3 !!!! This script takes a .trees file from SLiM 3.x and add neutral mutations to the simulated populations. ")
	parser.add_argument("--input","-i", 
			required = True,
			dest = "input",
			type =str, 
			help = "The bed file that you want to analyse")
	parser.add_argument("--wsize","-w", 
			required = True,
			dest = "wsize",
			type = int, 
			help = "The window size that you want to analyse")
	parser.add_argument("--maxLen","-m", 
			required = True,
			dest = "maxLen",
			type = int, 
			help = "The maximum lenggth that you want to extend to")
	parser.add_argument("--minLen", 
			required = False,
			dest = "minLen",
			type = int, 
			help = "The maximum lenggth that you want to extend to",
			default = 0)


	args = parser.parse_args()

	offsets = range(0, args.maxLen + args.wsize, args.wsize)
	#print offsets
	
	mouseChroms = [str(i+1) for i in range(19)]# + ['X']

	bed = readBed(args.input)

	clengths = {}
	for i in open('/home/booker/work/MouseAdap/testResources/mm10.fa.fai'):
		x  = i.strip().split()
		clengths[x[0]] = int(x[1])
	sbed = {} # A dictionary of the BED file that has been split by chromosome
	for mc in mouseChroms:
		sbed[mc] = bed[bed['chrom'] == 'chr'+mc].copy()
	
	for mc in mouseChroms:
		chrom = 'chr' + mc
		cbed = sbed[mc]
		n = cbed.shape[0] # the number of genomic elements corresponding to this particular chromosome  
		cend = clengths[chrom]
		if n == 0:
			continue ## If this condition is met, it means that there are no entries in the BED file for that chromosome

	## Let's get LOOOW - the lower threshold for each element
		cbed['lwr_lim'] =  cbed['chromStart'] - (cbed['chromStart'] - np.insert(np.array(cbed['chromEnd'][[i for i in range(0,(n-1))]]), 0, 0, axis=0))/2
	## Let's get HIGH - the lower threshold for each element
		cbed['upr_lim'] = cbed['chromEnd'] + (np.append(np.array(cbed['chromStart'][[i for i in range(1,(n))]]), cend) - cbed['chromEnd'])/2

	## Make all entries below the size threshold NaNs
		cbed[cbed['chromStart'] - cbed['lwr_lim'] < args.minLen] = np.nan
		cbed[cbed['upr_lim'] - cbed['chromEnd'] < args.minLen] = np.nan
		
	## Now let's drop all invalid entries from the dataframe
		cbed = cbed[cbed.chromStart.notnull()].copy()
		for offset in offsets:
			lower = getLwrWin(cbed.copy(), offset, args.wsize)
			upper = getUprWin(cbed.copy(), offset, args.wsize)
			d = pd.concat([lower,upper])
			d[d.chromStart.notnull()].copy()

		## At this ppoint you could output the files to a series of BED files for using with VCFs to obtain data on the diversity around functional elements
		## What I'm going to do is to take this basic structure and use it to analyse SLiM output

main()
