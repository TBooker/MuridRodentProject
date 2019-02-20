# Keywords: Python, tree-sequence recording, tree sequence recording

# This is a Python recipe; note that it runs the SLiM model internally, below

import subprocess, msprime, pyslim
import matplotlib.pyplot as plt
import numpy as np

# Run the SLiM model and load the resulting .trees
#subprocess.check_output(["slim", "-m", "-s", "0", "./recipe_16.4.slim"])
ts = pyslim.load("./recipe_16.4.trees").simplify()

## Let's add mutations to the tree
