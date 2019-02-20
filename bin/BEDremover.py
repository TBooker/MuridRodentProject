## This script takes a BED file that and subtracts a single base pair from the start of each element
## This is done because BED files store elements with an extra number at the beginning of each element 
## The output of this can be concatenated with other files and sorted to make the coombined elements file needed to make simulationconfig files
import argparse




def main():
	parser = argparse.ArgumentParser(description="Subtract he first base of each genomic element in a BED file")
	parser.add_argument("-i", 
			required = True,
			dest = "input",
			type = str, 
			help = "The input BED file")
	parser.add_argument("--name", 
			required = True,
			dest = "name",
			type = str, 
			help = "Name the element one of CDS, CNE or UTR")
	parser.add_argument("-o", 
			required = True,
			dest = "output",
			type = str, 
			help = "the name of the output file")
	parser.add_argument("--dont_change", 
			dest = "dont_change",
			action = 'store_true', 
			help = "use this flag if you don't want to change the positions, but do want to add the CDS/UTR/CNE label",
			default = False)
	args = parser.parse_args()
	output = open(args.output, 'w')

	for i in open(args.input):
		x = i.strip().split()
		if not x[0].startswith('c'):
			x[0] = 'chr'+x[0] 
		if not args.dont_change:
			x[1] = str(int(x[1])+1)
		
		output.write( '\t'.join( x[:3] +[args.name] ) +'\n')
	output.close()
main()