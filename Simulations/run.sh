
#The following function will generate a config file for a simulation of 1000 individuals, with a population-scaled mutation rate (4Neu) of 0.0083, with unform recombination rate of 0.009, DFEs as laid out in the testResources/DFE.txt file, for 10N = 10000 generations based on a 500 Kbp chunk of the mouse genome. 

python bin/SlimFunctions.py --N 1000 \
							--theta 0.0083 \
							--rho 0.009 \
							--DFE ../testResources/DFE.txt \
							--generations 10 \
							-o 1.config \
							--length 500000 \
							--treeSeq
							
# The output config file is stored as 1.config.txt

# Making lots of similar configs files is easy peasy, you can run the above in a loop and use the loop variable in the output name, or use GNU parallel which is what I tend to do. For example, to make 100 replicate simulation configs, you could run something like:

parallel "python bin/SlimFunctions.py --N 1000 --theta 0.0083 --rho 0.009 --DFE ../testResources/DFE.txt --generations 10 -o mySims/{}.config --length 500000 --treeSeq" ::: $(seq 1 100)

# Note that the output is stored in a directory now

# Running the simulations is also easy peasy, and so is parallelisation
parallel "/path/to/slim.3/slim  mySims/{}.txt" ::: $(seq 1 100)

# If you're running the treeSeq mode, which I can't recommend highly enough, you can sprinkle neutral mutations onto the stored coalescent simulation using the SprinkleMutations.py script

# For a single replicate, it would look something like:
python3 bin/SprinkleMutations.py --input mySims/1.txt \
								--trees mySims/1.trees \
								--sample 30 \
								-o mySims/1.vcf

# Note that the script gzips the resulting VCF file

# Again, it is very easy to parallelise this process:
parallel "python3 bin/SprinkleMutations.py --input mySims/{}.txt --trees mySims/{}.trees --sample 30 -o mySims/{}.vcf" ::: $(seq 1 100)

# In the MBE paper we analysed the patterns of diversity around functional elements. The next script that takes a config file and the corresponding VCF and returns a file containing the site frequency spectrum in windows around a genomic element type in your simulation. It makes use of tabix and bedtools, so you'll need to have both of those installed on your system in order to run the script.


python bin/getSLiMflanks.py --input configs/1.txt \
							--vcf mySims/1.vcf.gz \
							--wsize 100 \
							--maxLen 10000 \
							--element g2 \
							-o mySims/1

# This call will extract the SFS in windows of 100bp, extending to 10,000 bp in the flanks of all g2 elements in the simulation. In this case, g2 correpsonds to CNEs.

# You're probably getting sick of this, but yes it is very easy to parallelise this as well... 
parallel "python bin/getSLiMflanks.py --input mySims/{}.txt --vcf mySims/{}.vcf.gz --wsize 100 --maxLen 10000 --element g2 -o {}" ::: $(seq 1 100)

# A final step is to summarise the data coming from multiple replicate simualtions. This is done using as follows:

python bin/summariseSFSfromDir.py --input mySims/ --output CNEs --boots 100

# This script then reads in all the SFS files, and collates the SFSs by distance to functional element. Summary stats like nucleotide diversity (pi) and Tajima's D are calculated from the collated SFSs. Bootstraps are perfoemd (100 replicates in this case) and 95% intervals are saved from the resulting data 



# Disclaimer, I'm lucky enought to have access to servers that don't use schedulers like sungrid or the like, so I've not had to worry about writing submission scripts. If you were running simulations on a cluster or something, I would recommend making the configs locally, copying them over to the cluster, then running a submission script.
