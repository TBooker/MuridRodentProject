from Bio import SeqIO
#A quick script for splitting a FASTA into individual sequences

mm10_file = "/home/booker/work/MouseAdap/ESEs/dataFiles/mm10.fa"

print 'loading genome reference file dictionary'
mm10_dict = SeqIO.to_dict(SeqIO.parse(mm10_file, "fasta"))

for i in mm10_dict.keys():
	print i
	out = open(i+'.fasta', 'w')
	out.write('>' + i + '\n')
	out.write(str(mm10_dict[i].seq)  +'\n')
	out.close()
