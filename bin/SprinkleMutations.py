import msprime, pyslim, argparse, random
import numpy as np
import re

def remove_mutations(ts, start, end, proportion):
	'''
	This function will return a new tree sequence the same as the input,
	but after removing each non-SLiM mutation within regions specified in lists
	start and end with probability `proportion`, independently. So then, if we
	want to add neutral mutations with rate 1.0e-8 within the regions and 0.7e-8
	outside the regions, we could do
	ts = pyslim.load("my.trees")
	first_mut_ts = msprime.mutate(ts, rate=1e-8)
	mut_ts = remove_mutations(first_mut_ts, start, end, 0.3)

	:param float proportion: The proportion of mutations to remove.
	'''
	new_tables = ts.dump_tables()
	new_tables.mutations.clear()
	mutation_map = [-1 for _ in range(ts.num_mutations)]
	for j, mut in enumerate(ts.mutations()):
	   keep_mutation = True
	for i in range(len(start)):
		left = start[i]
		right = end[i]
		assert(left < right)
		if i > 0:
			 assert(end[i - 1] <= left)
		if mut.position >= left and mut.position < right and len(mut.metadata) == 0:
			 keep_mutation = (random.uniform(0, 1) > proportion)
		if keep_mutation:
			mutation_map[j] = new_tables.mutations.num_rows
		if mut.parent < 0:
			 new_parent = -1
		else:
			new_parent = mutation_map[mut.parent]
			new_tables.mutations.add_row(site = mut.site, node = mut.node,
			derived_state = mut.derived_state,
			parent = new_parent,
			metadata = mut.metadata)
	return new_tables.tree_sequence()


class SLiM:
	def __init__(self, configFile):
		regexString = "[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?"
		lines = [i.strip().strip(';') for i in open(configFile)]
		elementsRaw = [i.split('(')[-1].strip(')').split(' ') for i in lines if i.startswith('initializeGenomicElement(')]
		self.elementStarts = [int(i[1].split(',')[0]) for i in elementsRaw]
		self.elementEnds = [int(i[2].split(',')[0]) for i in elementsRaw]
#		self.mutRate = float([i for i in lines if i.startswith('initializeMutationRate(')][0].strip().strip(';').strip(')').split('(')[1])
		self.mutRate = float(re.findall(regexString, [i for i in lines if i.startswith('initializeMutationRate(')][0])[0])
		print(self.mutRate)


def main():
	parser = argparse.ArgumentParser(description="PYTHON 3 !!!! This script takes a .trees file from SLiM 3.x and add neutral mutations to the simulated populations. ")
	parser.add_argument("--input","-i", 
			required = True,
			dest = "input",
			type =str, 
			help = "name the simulation config file")
	parser.add_argument("--trees","-t", 
			required = True,
			dest = "trees",
			type =str, 
			help = "name the trees file")
#	parser.add_argument("--N", 
#			required = True,
#			dest = "N",
#			type = int, 
#			help = "The population size that was simulated")
	parser.add_argument("--sample", 
			required = True,
			dest = "sample",
			type = int, 
			help = "The number of haploid copies to sample from the population")
#	parser.add_argument("--theta", 
#			required = True,
#			dest = "theta",
#			type = float, 
#			help = "The population effective mutation rate")
	parser.add_argument("-o", 
			required = True,
			dest = "output",
			type =str, 
			help = "name the output file")
	args = parser.parse_args()

	config = SLiM(args.input)

	ts = pyslim.load(args.trees).simplify()

#	mut_rate = args.theta/(4.*args.N)

	sprinkled = msprime.mutate(ts, rate= config.mutRate, keep=True)

	mutated = remove_mutations(sprinkled, config.elementStarts, config.elementEnds, 0.0)


#	Here we strip out the neutral mutations that were added on top of elements that were already simulated

	haploidSample = np.random.choice(mutated.samples(),args.sample, replace = False)

	ts2 = sprinkled.simplify(samples = haploidSample)

	with open(args.output, "w") as vcf_file:
		ts2.write_vcf(vcf_file, 1)

	for ind in pyslim.extract_mutation_metadata(ts2.tables):
		print(ind)

	return
#	for tree in ts2.trees():
#		print("-" * 20)
#		print("tree {}: interval = {}".format(tree.index, tree.interval))
#		print(tree.draw(format="unicode"))


main()
