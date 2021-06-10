import sys, argparse, gzip

def main():
	parser = argparse.ArgumentParser(description="PYTHON 3 !!!! This script takes a .trees file from SLiM 3.x and add neutral mutations to the simulated populations. ")
	parser.add_argument("--input","-i", 
			required = True,
			dest = "input",
			type =str, 
			help = "name the input file, assumes that the file is Gzipped")
	args = parser.parse_args()

# BED format:
# chr1	3100870	3100938	CNE

# GTF format:
# chr1	some_text	CNE	3100871	3100938

## REMEMBER THAT BED FILES ARE 0-Based, or half open, or whatever it's called

	first = 10
	
	for i in open(args.input):
		
		x = i.split()
		chrom = x[0]
		start = int(x[1])
		end = int(x[2])
		classy = x[3]
		if first == 10:
			print '!'
			previousLine = [chrom,start, end, classy]
			first +=1
			continue
		else:
			if start - previousLine[1] < 0:
				print previousLine
				print i.strip()
				print start - previousLine[1]
				print
			previousLine = [chrom,start, end, classy]

main()
