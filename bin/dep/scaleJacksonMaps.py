import argparse
import pandas as pd
import numpy as np

##  Script to convert values from the Jackson Lab converter to cM/Mb values I can use in my simulations

def positionSampler(simLength, file = '/home/booker/work/MouseAdap/testResources/mm10.fa.fai'):
	invalid = True
	while invalid:
		chromLengths = []
		chroms = []
		chromDict = {}
		for line in open(file).readlines():
			x = line.strip().split()
			chroms.append( x[0] )
			chromLengths.append( float(x[1]) )
			chromDict[ x[0] ] = float(x[1]) 
		chosenChrom = np.random.choice(chroms, p =np.array(chromLengths)/sum(chromLengths))
		chosenPosition =  int(np.random.random() * (chromDict[chosenChrom] - 3e6)) + 3e6
		if chosenPosition + simLength < chromDict[chosenChrom]:
			invalid = False
		else:
			print 'reject', chosenChrom,  int(chosenPosition), int(chosenPosition + simLength) 
	return chosenChrom,  int(chosenPosition), int(chosenPosition + simLength)

def sampleLDmap(region):
	print region
	recFile ='/home/booker/work/MouseAdap/recombinationMaps/castaneus/'+region[0]+'.rho.condensed.txt.gz'
	rawMap = pd.read_csv(recFile, names=['chrom','lab1','lab2','pos','posRep','rho','cumRho'], sep = '\t', header = None)	
	print rawMap['chrom']
	
def main():
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("--input","-i", 
			required = True,
			dest = "input",
			type =str, 
			help = "The recombination maps for mouse chromosomes")
	parser.add_argument("--LD", 
			required = False,
			dest = "LD",
			action = 'store_true', 
			help = "Use this flag if the maps you are analysing were generated using LDhelmet")

	args = parser.parse_args()

#	if args.LD:
#		rawMap = pd.read_csv(args.input, names=['chrom','lab1','lab2','pos','posRep','rho','cumRho'], sep = '\t', header = False)	
#	[chr14,recom_rate,LDhelmet,7413308,7413308.,1.6541e-02,73000.511087]
#	else:
#		rawMap = pd.read_csv(args.input)	
		
	np.random.seed(123)
	region = positionSampler(100000)
	sampleLDmap(region)

main()
