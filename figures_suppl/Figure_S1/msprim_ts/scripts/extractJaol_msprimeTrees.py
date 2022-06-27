infile = str(snakemake.input)
samplingtimes = int(snakemake.params.nsam)
outfile = str(snakemake.output)

# script
import tskit, numpy as np

################################################################################
from heapq import nlargest
import random

def sample_from_iterable(iterable, samplesize):
    return (x for _, x in nlargest(samplesize, ((random.random(), x) for x in iterable)))

################################################################################

def ri(): return random.randint(0,1)

def getOnlyDiINDS(ts):
    return np.array([node.id for node in ts.nodes() if node.is_sample()]).flatten('C')

def getPairTs(ts):

    # get two random individuals
    i, j = sample_from_iterable(getOnlyDiINDS(ts), 2)

    # simplify tree by individuals
    return ts.simplify([i, j])

def getJaol(ts, samplingtimes):

    t_age, t_length = [], []

    # one loop per sampling
    for i in range(samplingtimes):

        # get a simplified pair tree
        pairTS = getPairTs(ts).simplify()

        # loop over all trees to get age and length for each one
        # remind: 1 tree for 1 pair of haploids, there is no recombination
        # within a single tree

        age, length = [], []

        for tree in pairTS.trees():
            length.append(tree.interval[1]-tree.interval[0])
            age.append(tree.time(tree.root))

        t_age.extend(age), t_length.extend(length)

    return np.array([t_age, t_length]).transpose()

################################################################################
# execution

# load tree; simplify!
ts = tskit.load(infile).simplify()

# get jaol (joint age on length)
jaol = getJaol(ts, samplingtimes)

# write table into file
np.savetxt(
	outfile,
	jaol,
	delimiter="\t",
	header="age\tlength"
	)
