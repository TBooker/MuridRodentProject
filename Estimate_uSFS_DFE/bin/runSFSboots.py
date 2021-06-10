import argparse, subprocess, random, os 
import numpy as np
from collections import Counter
from tqdm import tqdm

def commandLine(path_to_exec, config, input, seed, output_sfs, output_anc):
	return ' '.join([path_to_exec, config, input, seed, output_sfs, output_anc])

def main():

	parser = argparse.ArgumentParser(description="Bootstraps the SFS data and generates")

	parser.add_argument("--input",
				required = True,
				dest = "input",
				type =str,
				help = "the name of the file containing the data")

	parser.add_argument("--outputPrefix",
				required = True,
				dest = "outputPrefix",
				type =str,
				help = "A prefix for the output files, you'll get two per run")
	parser.add_argument("--output",
				required = True,
				dest = "output",
				type =str,
				help = "A prefix for the output files, you'll get two per run")

# 	parser.add_argument("--path_to_exec",
# 				required = True,
# 				dest = "path",
# 				type =str,
# 				help = "the full path to Peter's program")
# 
	parser.add_argument("--boots",
				required = True,
				dest = "boots",
				type = int,
				help = "The number of bootstrap windows you want to analyse")

 	parser.add_argument("--config",
  				required = True,
 				dest = "config",
  				type = str,
 				help = "The config file that you want to use for PK's program")

	args = parser.parse_args()
	
	bigDict = {}
	firstLine = True
	current = []	
	tempName = args.outputPrefix+'.tempDATA.out'
	temp = open(tempName, 'w')
	for i in open(args.input):
		if i.startswith('chr'):
			if firstLine:
				current = []
				firstLine = False
			else:
				bigDict[name] = current
				current = []
				name = i.strip()
			name = i.strip()
			continue
		else:
			current.append(i)
			temp.write(i)	
	temp.close()
	
	path_to_exec = '~/software/est-sfs-release-2.03/est-sfs'
	config_file = 'config-kimura.txt'
	seed_file = 'seedfile.txt'
	output_sfs = args.outputPrefix + '.sfs'
	output_p_anc = args.outputPrefix + '.anc.txt'
	
	names = np.array(bigDict.keys())
	numNames = len(names)

	Vanilla = [path_to_exec, config_file, tempName, seed_file, output_sfs, output_p_anc]

#	subprocess.call(Vanilla)

	os.system(commandLine(path_to_exec, config_file, tempName, seed_file, output_sfs, output_p_anc))
	
	finalFile = open(args.output, 'w') ## This will be the file that all the results will be written to
	
	[finalFile.write('vanilla,'+i) for i in open(output_sfs, 'r')]
	
	os.system('rm '+output_sfs)
	os.system('rm '+args.outputPrefix+'.tempDATA.out')
	
	for b in range(args.boots):
		print b+1
		sample = np.random.choice(names, numNames)
		temp = open(tempName, 'w')
		for j in sample:
			temp.write(''.join( bigDict[j] ) )
		temp.close()
		os.system(commandLine(path_to_exec, config_file, tempName, seed_file, output_sfs, output_p_anc))
		[finalFile.write(str(b+1)+','+i) for i in open(output_sfs, 'r')]
		os.system('rm '+output_sfs)
		os.system('rm '+args.outputPrefix+'.tempDATA.out')
		
	finalFile.close()
	
		
		
main()

