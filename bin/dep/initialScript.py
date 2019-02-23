import numpy as np
import pysam

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

## This function parses the DFE file. Take a look at the example file to see how I format those files
## I haven't yet added flexibility to do discrete DFEs or anything other than Gamma DFEs

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


## Do your elements bookend one another?

def bookend(e1,e2):
	if e1[2]  == e2[1] -1:

		return True
	else:
		return False

## This function grabs all the functional elements in the sampled region and adds the intergenic stuff 
## in between the features. It then returns a list of annotations that are going to be used for the 
## simulations in SLiM format
##

def genomicArchitechture(sampledRegion, nameDict, simLength, file = '/home/booker/work/MouseAdap/Annotations/combinedElements.s.bed.gz'):
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
			items.append(['g0',0, elements[el][1]-1])
			items.append(elements[el])
		elif el < len(elements)-1:
			if bookend(elements[el-1],elements[el]):
				items.append(elements[el])
			else:
				items.append(['g0',items[-1][2]+1, elements[el][1]-1])
				items.append(elements[el])
		elif el == len(elements)-1:
			if bookend(elements[el-1],elements[el]):
				items.append(elements[el])
				items.append(['g0',elements[el][2]+1, simLength - 1])
			else:
				items.append(['g0',items[-1][2]+1, elements[el][1]-1])
				items.append(elements[el])
				items.append(['g0',elements[el][2]+1, simLength - 1])
	genomicElements = []
	for k in items:
		genomicElements.append( 'initializeGenomicElement(' + str(k[0]) + ', ' + str(k[1]) +', ' + str(k[2])+ ');' )
	return genomicElements





## Now, using a BED file of all the functional elements in the genome, we generate the simulation architechture
simLength = 100000
theta = 0.01
rho = 0.01
N = 10000
#np.random.seed(123)

sampledRegion = positionSampler(simLength)
mutationTypes, genomicElementTypes, genomicElementKey = DFEparser('/home/booker/work/MouseAdap/testResources/DFE.txt', N)
formattedGenomicElements = genomicArchitechture(sampledRegion, genomicElementKey, simLength)

print 'initialize() {'
print 'initializeMutationRate(' + str(theta / (4*N)) + ');'
print 'initializeMutationType("m0", 0.5, "f", 0.0);'

for i in mutationTypes:
	print i
print 'initializeGenomicElementType("g0", m0, 1.0);'
for i in genomicElementTypes:
	print i
for i in formattedGenomicElements:
	print i
print 'initializeRecombinationRate(' + str(rho / (4*N)) + ');'
print '}'
print '1 { sim.addSubpop("p1", ' + str(N) + ' ); }'
print str(N*10) +' late() {p1.outputSample(10);}'