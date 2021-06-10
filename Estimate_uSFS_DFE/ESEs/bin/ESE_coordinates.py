from Bio import SeqIO
import re, sys
#A script for obtaining the coordinates of ESEs within protein coding regions of the genome

if len(sys.argv) < 2:
	print 'Usage: \npython ESE_coordinates chr1 \n You must specify which chromosome you are analysing'
	sys.exit()
	
if sys.argv[1].startswith('chr'):
	chromosome = sys.argv[1]
else:
	chromosome = 'chr' + sys.argv[1]
	
	
# specify the paths to each of the data files you'll use
CDS_file = "/home/booker/work/MouseAdap/ESEs/dataFiles/AllCDS.gtf"
mm10_path = "/home/booker/work/MouseAdap/ESEs/dataFiles/mm10_split/"
ESEs_file = "/home/booker/work/MouseAdap/ESEs/dataFiles/merged_no_Ke_ESEs.txt"

# Make a list of all the ESE motifs
print 'loading ESEs'

ESEs = [i.strip() for i in open(ESEs_file).readlines() if not i.startswith('m')]
## The if not statement at the end removes the header from the fileq

print 'loading genome reference file dictionary'

mm10_dict = SeqIO.to_dict(SeqIO.parse(mm10_path + chromosome + ".fasta", "fasta"))

output = open('ESEcoords.' + chromosome + '.bed','w')

print 'iterating over GTF file:', CDS_file

count = 0
for i in open(CDS_file):
	count +=1
	x = i.strip().split()
	chrom = 'chr' + x[0]
	if chrom != chromosome: continue
	start = int(x[3])-1
	end = int(x[4])-1
	CDSsequence = str(mm10_dict[chrom].seq[start:end]).upper()
	for ESE in ESEs:
		temp = [(m.start(0) + start - 1, m.end(0)+start) for m in re.finditer(ESE, CDSsequence)]
		for match in temp:
			output.write(chrom + '\t' + str(match[0]) + '\t' + str(match[1]) + '\t' + str(match[1] -match[0]) + '\n')
	if count % 1000 == 0:
		print 'processed',count,'records'

output.close()
