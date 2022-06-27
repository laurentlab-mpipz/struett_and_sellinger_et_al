import sys
import datetime
import warnings
import itertools
import pickle
import math
import numpy as np
import pandas as pd
import tskit
import json
import pyfuncs  # from file


# read sample names
sample_names = np.loadtxt(snakemake.input.sample_list, dtype=str)


# read tree sequence
treeseq_athal = tskit.load(snakemake.input.athal_treeseq)


# find the node ids for the sample of the population
population_sample = []
for individual in treeseq_athal.individuals():
    if str(json.loads(individual.metadata)["id"]) in sample_names:
        population_sample.extend(individual.nodes)


# sample treeseq to provided samples
treeseq_athal_population = treeseq_athal.simplify(samples=population_sample)
del treeseq_athal
