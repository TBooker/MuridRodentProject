## Simulating the effects of selection at linked sites in simulations incorporating realistic functional architechture and recombination rate variation

Here are scripts that can be used to generate SLiM 3.x config files incorporating recombination rate variation, arbitrary distributions of fitness effects and functional element architechture.

Each simulation replicate samples the genome and 

###The things you need are:

⋅⋅* Annotations of the functional elements that you want to model

I made a combined BED file that had the CDS and UTRs of protein-coding genes (making no distinction between 3' and 5' UTRs) as well as conserved non-coding elements (CNEs) obtine dby Halligan et al (2013). The simulations sample this file, and create simulation architechture based on the elements within a chosen window.
 
To see how I generated a list of the functional elements in the mouse genome, see the file ``run.sh`` in the Annotations/ directory.

⋅⋅* Estimates of the distribution of fitness effects for each of those elements

For each of the functional element types I obtained estimates of the full DFE (both +ve and -ve fitness effects) by analysis of the site frequency spectrum. To see how that was done, check out the DFE/ directory

⋅⋅* Recombination rate map(s)


___

## Notes

I've recently incorporated Tree Sequence recording stuff from SLiM 3.X into my config files. Using this option speeds up the simulation like the Dickens. For illustration, during my PhD the simulations I was performing in SLiM 1.x, which went into the MBE paper, took about 40 minutes each. Performing similar simulations (500Kbp chromosomes in 1000 individuals for 10,000 generations) in SLiM 3.x with TreeSeq, takes less than a minute (!!!) and post-processing takes very little time too.  

In the MBE paper (i.e. Booker and Keightley 2018), we used variant calls, functional annotations and LD-based recombination rate estimates that were generated using the mm9 version of the mouse genome. We are currently re-doing a lot of similar analyses usng mm10-based datasets, so get in touch if you want access to any of that stuff, but note that it is a work in progress. Other than the simulations, all of the raw data we used for the MBE paper was published in other studies, but I am more than happy to provide directions to that data.

I've tried not to hard code things into my scripts, but my laziness may have gotten the better of me in some places. I'm happy to help anyone who may want to use these scripts, so get in touch if you have any questions.
