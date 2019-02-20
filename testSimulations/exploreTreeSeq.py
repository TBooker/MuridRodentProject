import msprime, pyslim
import numpy as np

def ld_matrix_example(ts):
    ld_calc = msprime.LdCalculator(ts)
    A = ld_calc.r2_matrix()
    # Now plot this matrix.
    x = A.shape[0] / pyplot.rcParams['figure.dpi']
    x = max(x, pyplot.rcParams['figure.figsize'][0])
    fig, ax = pyplot.subplots(figsize=(x, x))
    fig.tight_layout(pad=0)
    im = ax.imshow(A, interpolation="none", vmin=0, vmax=1, cmap="Blues")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in 'top', 'bottom', 'left', 'right':
        ax.spines[s].set_visible(False)
    pyplot.gcf().colorbar(im)
    pyplot.savefig("ld.svg")


ts = pyslim.load("./recipe_16.4.trees").simplify()
mutated = msprime.mutate(ts, rate=1e-7, random_seed=1, keep=True)
NewData = np.random.choice(mutated.samples(),20, replace = False)
print(NewData)
ts2 = mutated.simplify(samples = NewData)

with open("output.vcf", "w") as vcf_file:
	ts2.write_vcf(vcf_file, 1)

#for variant in mutated.variants():
#	print(variant.site.id, variant.site.position,
#	variant.alleles, variant.genotypes, sep="\t")
#import msprime
#import matplotlib.pyplot as pyplot

#ld_matrix_example(mutated)

for tree in ts2.trees():
	print("-" * 20)
	print("tree {}: interval = {}".format(tree.index, tree.interval))
	print(tree.draw(format="unicode"))
