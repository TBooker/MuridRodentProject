import argparse, math, re, subprocess, vcf
import numpy as np
import pandas as pd

## This script is written assuming Python 2.7

## This script relies on a number of modules that you may not have installed on
## the system you are using. The pandas and numpy library often don't come 
## as standard with most python distros (even though they should!). 
## Also, you need PyVCF, a very handy package.

## Finally, BEDtools is called from within the script, so you'll need to
## have that installed and callable from any given directory
## I'm running BEDtools version v2.27.1

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
		
		
		
def SFS_from_frequencies(frequencies, length,N):
	SFS = [0]*(N+1)
	for i in frequencies:
		if i > N:
			print "SFS_from_frequencies: Error in your frequencies vector: One of the values is greater than the number of individuals\nThe offending value is: " + str(i) +" and the sample is "+str(N)
			return
		SFS[i] += 1
	SFS[0] = length - len(frequencies)
	if sum(SFS) < length:
		print"SFS_from_frequencies: Error in your frequencies vector: Fewer items in the SFS than the length of the region"
		return
	if sum(SFS) > length:
		print"SFS_from_frequencies: Error in your frequencies vector: More items in the SFS than the length of the region"
		return
	return SFS
		
		
		
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
	parser = argparse.ArgumentParser(description="PYTHON 3 !!!! This script takes a simulation config files and a processed VCF and returns a list of site frequnecy spectra for a given class of functional elements. The script currently assumes that the VCF is in haploid format!")
	parser.add_argument("--input","-i", 
			required = True,
			dest = "input",
			type =str, 
			help = "The config file of the simulation that you want to analyse")
	parser.add_argument("--vcf","-v", 
			required = True,
			dest = "vcf",
			type =str, 
			help = "The VCF of mutations in the simulation that you want to analyse. It has to be BGZIPPED and TABIXED")
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
			help = "The minimum length an analysis window can be, the default is 2",
			default = 0)
	parser.add_argument("--element", 
			required = True,
			dest = "element",
			type = str, 
			help = "Which element are you analysing?")
	parser.add_argument("--output","-o", 
			required = True,
			dest = "output",
			type =str, 
			help = "Give a prefix for your output files. This script creates a number of temporary files, so use something that won't cause problems if you are running this script in parallel")
	parser.add_argument("--dont_clean", 
			required = False,
			dest = "dont_clean",
			action = 'store_true', 
			help = "Give this flag if you don't want the script to remove the temporary files. You might want to do this for testing etc.")


	args = parser.parse_args()
		
# Prepare the list of analysis windows that you'll be using.
	offsets = range(0, args.maxLen + args.wsize, args.wsize)

# Read the SLiM config file in. Remember that Eidos is a coding language 
# so there are a very large number of ways that a person could code up a 
# simulation. The SLiM class that I've put together depends VERY strongly
# on the way that the function in SLiMfunctions.py prepares them

	slim = SLiM(args.input)
# Make a bed file of the focal functional elements in the simulation
	cbed = slim.SLiMelementBed(args.element)
# Make a bed file of all the other functional elements, used for masking
	fullBed = slim.SLiMaltElementBed(args.element)

# How long is the simuated chromosome?
	cend = max(slim.elementEnds) # Get this from the SLiM file 

# How many of the focal elements are there in the simulation?
	n = cbed.shape[0] # the number of genomic elements corresponding to this particular chromosome  

	if n == 1 and cbed['chrom'][0] == 0:
		print 'There are none of the focal element type in your simulation'
		return ## If this condition is met, it means that there are no entries in the BED file for that chromosome

# For each offset value (i.e. moving away from 
# the edge of a focal element), get the analysis windows and store them 
# as a BEDframe
	temp = []
	for offset in offsets:
		tempBed = gimmeFlanks(cbed, offset, cend, args.minLen, args.wsize)
		temp.append(tempBed)
	df = pd.concat(temp).reset_index()
	df = df[pd.notnull(df['chromStart'])].sort_values(by=['window']).copy()
	df = df.drop(['index'], axis=1)
	
## I'm now going to write two dataframes to BED files. The first is 
## the list of all functional elements in the simulation
## The second is the analysis windows that came out of this simulation
## I then subtract the functional elements from the analysis windows 
## using BEDtools in order to prevent them from being analysed 
## I came to this decision aftr trying a bunch of different subtraction
## methods, but I cannot beat bedtools for speed.

	fullBed.to_csv( args.output+'.func.bed',sep = '\t',index = False, header = False)
	
	df['chromStart']=pd.to_numeric(df.chromStart, downcast='integer')
	df['chromEnd']=pd.to_numeric(df.chromEnd, downcast='integer')
	df.to_csv(args.output + '.analysisWindows.bed',sep = '\t',index = False, header = False)
	
# store the header of the file that will be subtracted
	names_to_keep = list(df)
	
# prepare the bedtools command that will be performed using SubProcesss
	bedtools_command = ['bedtools', 'subtract', '-a', args.output + '.analysisWindows.bed','-b',  args.output+'.func.bed']
	print bedtools_command
	tempBedfileName = args.output + '.analysisWindows.subtracted.bed'	
	with open(tempBedfileName, "w") as outfile:
		subprocess.call(bedtools_command, stdout= outfile)

## And now read back in the bedfile of analysis windows 		
	subtractedBed = pd.read_csv(tempBedfileName, names = names_to_keep, sep = '\t')
	
## Now load up the VCF file
	vcf_reader = vcf.Reader(open(args.vcf, 'r'))
	numAlleles = 0
	
	finalOutput = []
## Let's iterate over all the analysis windows and extract the relevant bits from the VCF to collate the SFS
	for index, w_row in subtractedBed.iterrows():	
		if  w_row.chromEnd -  w_row.chromStart < 2:continue
		windowWidth = w_row.chromEnd -  w_row.chromStart
		alleleFreqs = []
		for i in vcf_reader.fetch(w_row.chrom.split('r')[1], w_row.chromStart, w_row.chromEnd):
			if numAlleles == 0 :
				numAlleles = i.num_called
			alleleFreqs.append(i.num_hom_alt)
		## Now let's collate the SFS from the alleles present in the current analysis window
		sfs = SFS_from_frequencies(alleleFreqs, windowWidth,numAlleles)
		finalOutput.append( [w_row.window, '1:'+str(w_row.chromStart)+'-'+str(w_row.chromEnd)] + sfs)
# Write the final data to a gzipped file
	pd.DataFrame(finalOutput).to_csv(args.output + '.sfs.txt.gz', header= False,index = False, compression='gzip')


## Now we need to clean up and scrub away all those pesky intermediate files
	tempFiles = [ args.output+'.func.bed', args.output + '.analysisWindows.bed', args.output + '.analysisWindows.subtracted.bed']
	if not args.dont_clean:
		for  i in tempFiles:
			clean_command = ['rm', i]
			subprocess.call(clean_command)

# list of files created by this script:
#
#'tempData.bed'
#'myfile'
#'tempData.func.bed'
#args.output+'.func.bed'
		
	return

main()
