import sys, argparse, gzip

def main():
	parser = argparse.ArgumentParser(description="PYTHON 3 !!!! This script takes a .trees file from SLiM 3.x and add neutral mutations to the simulated populations. ")
	parser.add_argument("--input","-i", 
			required = True,
			dest = "input",
			type =str, 
			help = "name the input file, assumes that the file is Gzipped")
	parser.add_argument("--output","-o", 
			required = True,
			dest = "output",
			type =str, 
			help = "name the output file")
	parser.add_argument("--element", 
			required = True,
			dest = "element",
			type =str, 
			help = "Give the code for the element")
	args = parser.parse_args()

# BED format:
# chr1	3100870	3100938	CNE

# GTF format:
# chr1	some_text	CNE	3100871	3100938

## REMEMBER THAT BED FILES ARE 0-Based, or half open, or whatever it's called
	output = open(args.output, 'w')

	for i in gzip.open(args.input):
		x = i.split()

		chrom = x[0]
		if not chrom.startswith('c'):
			chrom = 'chr'+chrom
		start = int(x[3]) - 1
		end = int(x[4])
		element = args.element
		line = '\t'.join([chrom, str(start), str(end), element]) + '\n'
		output.write(line)
	output.close()

main()
