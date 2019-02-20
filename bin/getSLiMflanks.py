## A script that is a Python rewrite of Dan's get_bed_flanks.R script.

import argparse, math, re, subprocess
import numpy as np
import pandas as pd
from interval import Interval, IntervalSet

class SLiM:
	def __init__(self, configFile):
		regexString = "[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?"
		lines = [i.strip().strip(';') for i in open(configFile)]
		elementsRaw = [i.split('(')[-1].strip(')').split(' ') for i in lines if i.startswith('initializeGenomicElement(')]
		self.elements = [i[0].split(',')[0] for i in elementsRaw]
		self.elementStarts = [int(i[1].split(',')[0]) for i in elementsRaw]
		self.elementEnds = [int(i[2].split(',')[0]) for i in elementsRaw]
#		self.mutRate = float([i for i in lines if i.startswith('initializeMutationRate(')][0].strip().strip(';').strip(')').split('(')[1])
		self.mutRate = float(re.findall(regexString, [i for i in lines if i.startswith('initializeMutationRate(')][0])[0])
		self.recom = [i for i in lines if i.startswith('initializeRecombinationRate(')]

	#	self.recRate = re.findall(regexString, [i for i in lines if i.startswith('initializeRecombinationRate(')])
		
	def SLiMelementBed(self, element):
		myBed = []
		if element in ['g1', 'g2', 'g3']: # I'll add strandedness to this in a bit
			for i,j,k in zip(self.elements, self.elementStarts, self.elementEnds):
				if i != element: continue
				myBed.append( ['chr1', j, k, '+','chr1:'+str(j)+'-'+str(k)])
			return pd.DataFrame(myBed, columns = ["chrom","chromStart","chromEnd","strand", "name"])
		else:
			return pd.DataFrame([[0,0,0,0,0]], columns = ["chrom","chromStart","chromEnd","strand", "name"])


	def SLiMaltElementBed(self, element):
		myBed = []
		if element in ['g1', 'g2', 'g3']: # I'll add strandedness to this in a bit
			for i,j,k in zip(self.elements, self.elementStarts, self.elementEnds):
				if i == element: continue
				myBed.append( ['chr1', j, k, '+','chr1:'+str(j)+'-'+str(k)])
			return pd.DataFrame(myBed, columns = ["chrom","chromStart","chromEnd","strand", "name"])
		else:
			return pd.DataFrame([[0,0,0,0,0]], columns = ["chrom","chromStart","chromEnd","strand", "name"])


	def SLiMBed(self):
		myBed = []
		for i,j,k in zip(self.elements, self.elementStarts, self.elementEnds):
			myBed.append( ['chr1', j, k, '+'])
		return pd.DataFrame(myBed, columns = ["chrom","chromStart","chromEnd","strand"])

	def recomMap(self):
		## This function is very reliant on the way that I've set up my SLiM simulations, so you'll have to double check that is works in al circumstances
		rawDat = [ [ i,j ] for i, j in zip( map(float,self.recom[0].split('(')[2].split(')')[0].split(',')), map(int, self.recom[0].split('(')[3].split(')')[0].split(','))) ]
		return pd.DataFrame(rawDat, columns = ['rate', 'endpoint'])
		
		
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
	x['window'] = x['direction'] + ':' + str(offset) + '-' + str(offset + length)

	return x
	
def makeNewRows(win, elements):
	print list(win), list(ele)
	A = IntervalSet([Interval(win['chromStart'],win['chromEnd'])])
	B = IntervalSet([Interval(ele['chromStart'],ele['chromEnd']) for ele in elements])
	for i in A - B :
		temp = win.copy()
		temp['chromStart'] = i.lower_bound
		temp['chromEnd'] = i.upper_bound
		
	
def masker(windows, elements):
	masked = []
	for index, w_row in windows.iterrows():	
		overlappers = []
		for index, e_row in elements.iterrows():
			
			if overlap( w_row['chromStart'], w_row['chromEnd'], e_row['chromStart'], e_row['chromEnd']):
				newRows = overlappers.append(e_row)




def overlap(start1, end1, start2, end2):
	"""Does the range (start1, end1) overlap with (start2, end2)?"""
	return end1 >= start2 and end2 >= start1

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
	x['window'] = x['direction'] + ':' + str(offset) + '-' + str(offset + length)
	return x
	
def gimmeFlanks(bed, offset, cend, minLen, wsize):
	n = bed.shape[0]
	## Let's get LOOOW - the lower threshold for each element
	bed['lwr_lim'] =  bed['chromStart'] - (bed['chromStart'] - np.insert(np.array(bed['chromEnd'][[i for i in range(0,(n-1))]]), 0, 0, axis=0))/2
## Let's get HIGH - the lower threshold for each element
	bed['upr_lim'] = bed['chromEnd'] + (np.append(np.array(bed['chromStart'][[i for i in range(1,(n))]]), cend) - bed['chromEnd'])/2
	## Make all entries below the size threshold NaNs
	bed[bed['chromStart'] - bed['lwr_lim'] < minLen] = np.nan
	bed[bed['upr_lim'] - bed['chromEnd'] < minLen] = np.nan
	
	## Now let's drop all invalid entries from the dataframe
	bed = bed[bed.chromStart.notnull()].copy()

	lower = getLwrWin(bed.copy(), offset, wsize)
	upper = getUprWin(bed.copy(), offset, wsize)
	d = pd.concat([lower,upper])
	d[d.chromStart.notnull()].copy()
	return d

	
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
	parser.add_argument("--element", 
			required = True,
			dest = "element",
			type = str, 
			help = "Which element are you analysing?")


	args = parser.parse_args()
		
	offsets = range(0, args.maxLen + args.wsize, args.wsize)

	slim = SLiM(args.input)
	cbed = slim.SLiMelementBed(args.element)
	fullBed = slim.SLiMaltElementBed(args.element)
	fullBed.to_csv('tempData.func.bed',sep = '\t',index = False, header = False)
	#print fullBed

	cend = max(slim.elementEnds) # Get this from the SLiM file 
	
	n = cbed.shape[0] # the number of genomic elements corresponding to this particular chromosome  

	if n == 1 and cbed['chrom'][0] == 0:
		print 'There are none of the focal element type in your simulation'
		return ## If this condition is met, it means that there are no entries in the BED file for that chromosome

	temp = []
	for offset in offsets:
		tempBed = gimmeFlanks(cbed, offset, cend, args.minLen, args.wsize)
		temp.append(tempBed)
	df = pd.concat(temp).reset_index()
	df = df[pd.notnull(df['chromStart'])].sort_values(by=['window']).copy()
	df = df.drop(['index'], axis=1)

	df['chromStart']=pd.to_numeric(df.chromStart, downcast='integer')
	df['chromEnd']=pd.to_numeric(df.chromStart, downcast='integer')
	df.to_csv('tempData.bed',sep = '\t',index = False, header = False)
	names_to_keep = list(df)
	
	command = subprocess.Popen(['bedtools', 'subtract', '-a', 'tempData.bed','-b', 'tempData.func.bed'],stdout=subprocess.PIPE)
	output = command.communicate()[0].split('\n')
### WRITE THE BED TO FILE AND THEN RELOAD INTO MEMORY HERE
	zz = pd.DataFrame([j.strip().split() for  j in output], columns = names_to_keep)
	zz['length'] = zz['chromEnd'] - zz['chromStart']
	print zz
	

	tempBed = gimmeFlanks(cbed, 0, cend, args.minLen, args.wsize)
	temp.append(tempBed)
	df = pd.concat(temp).reset_index()
	df = df[pd.notnull(df['chromStart'])].sort_values(by=['window']).copy()


main()
