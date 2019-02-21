import numpy as np
import pysam, argparse

## Here the script grabs a position from the mouse genome at random, but with respect to the fact that
## different chromosomes have different lengths
## The file that the script assumes is a faidx format file. That is, an index of a genome assemble fasta file
## Each line has a chromosome and it's length in base pairs (note that in mmX assemblies,
## there are 3e6 Ns representing unassembled heterochromatin)

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


## Now, we read in a file that contains the DFE information for each of the genomic elements featured in your BED file
## This parser assumes that there are only deleterious mutations and that all the DFE are gamma distributions
## The parse also assumes that there are 1/3 neutral sites in the CDS class. However, this is a variable and you can change it if you wish

class genomicElement:
	def __init__(self, line):
		self.name = line[0]
		if line[1] == 'g': # it's a gamma DFE
			self.shape = float(line[3])
			self.mean = float(line[2])
	def SlimSyntaxMutType(self, N, index):
		return 'initializeMutationType("m' + str(index) + '", 0.5, "g", ' + \
		 str(self.mean/N) + ', ' + str(self.shape) +');'
	def SlimSyntaxGenElType(self, index):
			if self.name == 'CDS':
				return 'initializeGenomicElementType("g' + str(index) + '", c(m0, m' + \
		 str(index) + '), c(0.333, 0.664));'

			else:
				return 'initializeGenomicElementType("g' + str(index) + '", m' + \
		 str(index) + ', 1.0);'

def DFEparser(dfeFile, N, CDSneutral = 0.333):
	mutTypes = [] ## Slim mutation types
	genElTypes = [] ## genomic element types
	nameDict = {}
	DFEs = []
	for line in open(dfeFile).readlines():
		DFEs.append( genomicElement( line.strip().split() ) )	

	for i in range(len(DFEs)):
		gX = DFEs[i]
		nameDict[gX.name] = 'g'+str(i+1)
		mutTypes.append( gX.SlimSyntaxMutType(N, i+1) )
		genElTypes.append( gX.SlimSyntaxGenElType(i+1) )
	return mutTypes, genElTypes, nameDict


## This function grabs all the functional elements in the sampled region and adds the intergenic stuff 
## in between the features. It then returns a list of annotations that are going to be used for the 
## simulations in SLiM format
##

def genomicArchitechture(sampledRegion, nameDict, simLength, file = '/home/booker/work/MouseAdap/Annotations/combinedElements.s.bed.gz', treeSeq = False):
	go = pysam.Tabixfile(file)
	elements = []
	for i in go.fetch(sampledRegion[0],int(sampledRegion[1]),int(sampledRegion[2])):
		x = i.strip().split()
		x[1] = int(x[1])
		x[2] = int(x[2])
		if x[1] < int(sampledRegion[1]):
			continue
		elif x[2] > sampledRegion[2]:
			continue
		elements.append([nameDict[x[3]], x[1] -sampledRegion[1], x[2]-sampledRegion[1]])
	items = []
	for el in range(len(elements)):
		if el == 0:
			if not treeSeq:
				items.append(['g0',0, elements[el][1]-1])
			items.append(elements[el])
		elif el < len(elements)-1:
			if not treeSeq:
				items.append(['g0',items[-1][2]+1, elements[el][1]-1])
			items.append(elements[el])
		elif el == len(elements)-1:
			if not treeSeq:
				items.append(['g0',items[-1][2]+1, elements[el][1]-1])
			items.append(elements[el])
			if not treeSeq:
				items.append(['g0',elements[el][2]+1, simLength - 1])
	genomicElements = []
	for k in items:
		genomicElements.append( 'initializeGenomicElement(' + str(k[0]) + ', ' + str(k[1]) +', ' + str(k[2])+ ');' )
	return genomicElements

def parseDemography(demography):
	print demography
	sizes = []
	times = []
	for line in open(demography):
		if line.startswith('N1'):continue
		elif line.startswith('N'):
			sizes.append( float(line.strip().split()[1]) )
		elif line.startswith('t'):
			times.append( float(line.strip().split()[1]) )
	for i,j in zip(sizes, times):
		yield [i,j]


def main():
	parser = argparse.ArgumentParser(description="Generate a simulation config file for SLiM")
	parser.add_argument("--N", 
			required = True,
			dest = "N",
			type = int, 
			help = "The population size to be simulated")
	parser.add_argument("--theta", 
			required = True,
			dest = "theta",
			type = float, 
			help = "The population effective mutation rate")
	parser.add_argument("--length", 
			required = True,
			dest = "length",
			type = int, 
			help = "The length of the simulated chromosome")
	parser.add_argument("--rho", 
			required = False,
			dest = "rho",
			type = float, 
			help = "The population effective recombination rate. this is optional, if you specify it it will set a uniform recombination rate")
	parser.add_argument("--DFE", 
			required = False,
			dest = "DFE",
			type = str, 
			help = "Give the file that contains the DFE information")
	parser.add_argument("--treeSeq", 
			required = False,
			dest = "treeSeq",
			action = 'store_true', 
			help = "Use this flag if you want to make a config file that will be in the TreeSeq format")
	parser.add_argument("--generations", 
			required = True,
			dest = "generations",
			type = int, 
			help = "the multiple of the popuation size that you want to simualte for. If specifying a non-equlbruum demography, the number of generation will be the number you give + N generations ")
	parser.add_argument("-o", 
			required = True,
			dest = "output",
			type =str, 
			help = "name the output file")
	parser.add_argument("--demography", 
			required = False,
			dest = "demography",
			type = str, 
			help = "Give a file containing the demographic history of the population",
			default = 'Null')
	parser.add_argument("--seed", 
			required = False,
			dest = "seed",
			type = int, 
			help = "Give a random seed for reproducible runs",
			default = 0)
	args = parser.parse_args()

## Now, using a BED file of all the functional elements in the genome, we generate the simulation architechture
	simLength = args.length
	theta = args.theta
	rho = args.rho
	N = args.N
	if args.seed == 0:
		pass
	else:
		np.random.seed(args.seed)

	generations = args.N * args.generations
	sampledRegion = positionSampler(simLength)
	mutationTypes, genomicElementTypes, genomicElementKey = DFEparser('/home/booker/work/MouseAdap/testResources/DFE.txt', N)
	if args.treeSeq:
		formattedGenomicElements = genomicArchitechture(sampledRegion, genomicElementKey, simLength,treeSeq =True)
	else:
		formattedGenomicElements = genomicArchitechture(sampledRegion, genomicElementKey, simLength,treeSeq =False)

	output = []
	output.append( 'initialize() {')
	output.append(  'initializeTreeSeq();')

	output.append( 'initializeMutationRate(' + str(theta / (4.*args.N)) + ');')
	output.append( 'initializeMutationType("m0", 0.5, "f", 0.0);')
	for i in mutationTypes:
		output.append( i )
		
	if not args.treeSeq:
		output.append( 'initializeGenomicElementType("g0", m0, 1.0);' )

	for i in genomicElementTypes:
		output.append( i )
	for i in formattedGenomicElements:
		output.append( i )
	output.append( 'initializeRecombinationRate(' + str(args.rho / (4*N)) + ');' )
	output.append( '}' )
	output.append( '1 { sim.addSubpop("p1", ' + str(N) + ' ); }' )

	if args.demography == 'Null':

		if args.treeSeq:
			output.append( str(generations) + ' late() {sim.treeSeqOutput("'+args.output+'.trees");}' )
		else:
			output.append( str(generations) + ' late() {p1.outputSample(10);}')
		config = open(args.output+'.txt', 'w')

		config.write('\n'.join(output))
		config.close()

	else:
		output.append( str(generations) + ' late() {sim.treeSeqOutput("'+args.output+'.eq.trees");}' )

		generationsAdded = 0
		for time_point in parseDemography(args.demography):
			output.append( str(generationsAdded + generations) + ' { p1.setSubpopulationSize(' + str(int(N*time_point[0])) + '); }' )
			generationsAdded += int(time_point[1]*N)

		finalGeneration =  int(time_point[1]*N) + (generationsAdded + generations) - int(time_point[1]*N)
		if args.treeSeq:
			output.append( str(finalGeneration) + ' late() {sim.treeSeqOutput("'+args.output+'.neq.trees");}' )
		else:
			output.append( str(finalGeneration) + ' late() {p1.outputSample(10);}')
		config = open(args.output+'.txt', 'w')

		config.write('\n'.join(output))
		config.close()
####

#	print 'initialize() {'
#	print 'initializeTreeSeq();'
#
#	print 'initializeMutationRate(' + str(theta / (4.*args.N)) + ');'
#	print 'initializeMutationType("m0", 0.5, "f", 0.0);'
#	for i in mutationTypes:
#		print i
#		
#	if not args.treeSeq:
#		print 'initializeGenomicElementType("g0", m0, 1.0);'
#
#	for i in genomicElementTypes:
#		print i
#	for i in formattedGenomicElements:
#		print i
#	print 'initializeRecombinationRate(' + str(args.rho / (4*N)) + ');'
#	print '}'
#	print '1 { sim.addSubpop("p1", ' + str(N) + ' ); }'
#	if args.treeSeq:
#		print '20000 late() {sim.treeSeqOutput("'+args.output+'.trees");}'
#
#	else:
#		print '20000 late() {p1.outputSample(10);}'


main()
