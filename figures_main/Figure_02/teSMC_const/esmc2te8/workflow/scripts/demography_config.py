
"""
provide demography as needed for the input of the generate ts data pipeline
needed: 2 lists, 1. Ne and time, 2. sigma and time, usually t_0 = 0. The lists
do not have to have the same time segmentation; we provide the demography as
piecewise constant,
e. g.
 - [
        [Ne_0, t_0],
        [Ne_1, t_1],
        ...
        [Ne_k, t_k]
    ]
"""

import numpy as np
import msprime # i am using this to cast the zigzag of thibaut

Ne = 100_000
sigma_0 = 0.95
sigma_1 = 0
t_sigma = [0.01, 1, 2.5]
r = 3.6e-8
mu = 4e-9
# special demographies
zigzag_amplitude = [5]
bottleneck_time = [0.01, 0.5, 2.5]
bottleneck_length = 0.05
bottleneck_strength = [3, 10]

# adjust the parameters to your satisfaction
time_segmentation_n = 100 # segment it over time in log Ne space
time_limit = 100 * Ne
def gen_log_space(limit, n):
    result = [1]
    if n>1:  # just a check to avoid ZeroDivisionError
        ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    while len(result)<n:
        next_value = result[-1]*ratio
        if next_value - result[-1] >= 1:
            # safe zone. next_value will be a different integer
            result.append(next_value)
        else:
            # problem! same integer. we need to find next_value by artificially
            # incrementing previous value
            result.append(result[-1]+1)
            # recalculate the ratio so that the remaining values will scale
            # correctly
            ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    # round, re-adjust to 0 indexing (i.e. minus 1) and return np.uint64 array
    return np.array(list(map(lambda x: round(x)-1, result)), dtype=np.uint64)
time_segmentation = gen_log_space(time_limit, time_segmentation_n)

piecewise_Ne = []
piecewise_sigma = []

# all constant; neg control
piecewise_Ne.append([[Ne, 0]])
piecewise_sigma.append([[sigma_0, 0]])

# const Ne, change in sigma
for ts in t_sigma:
    piecewise_Ne.append([[Ne, 0]])
    piecewise_sigma.append([[sigma_0, 0], [sigma_1, ts*Ne]])

# zigzag; here we use the time segmentation
def zigzag_function(time_segmentation=100, my_amplitude, Ne_0):
    """
    zig zag function with amplitude
    returns Ne for all the time segments, we use the msprime debugger to do so
    """
    # for each time segment find start and end Ne
    Ne = 1e4 # is the scaling parameter for the growth rate. The growth rate should generally be independent of the actual Ne
    assert my_amplitude == 5, "the growth rates are hardcoded for amplitude 5"
    demographic_events = [
                msprime.PopulationParametersChange(
                    time=20, growth_rate=6437.7516497364/(4*Ne)
                ),
                msprime.PopulationParametersChange(
                    time=30, growth_rate=-378.691273513906/(4*Ne)
                ),
                msprime.PopulationParametersChange(
                    time=200, growth_rate=-643.77516497364/(4*Ne)
                ),
                msprime.PopulationParametersChange(
                    time=300, growth_rate=37.8691273513906/(4*Ne)
                ),
                msprime.PopulationParametersChange(
                    time=2000, growth_rate=64.377516497364/(4*Ne)
                ),
                msprime.PopulationParametersChange(
                    time=3000, growth_rate=-3.78691273513906/(4*Ne)
                ),
                msprime.PopulationParametersChange(
                    time=20000, growth_rate=-6.4377516497364/(4*Ne)
                ),
                msprime.PopulationParametersChange(
                    time=30000, growth_rate=0.378691273513906/(4*Ne)
                ),
                msprime.PopulationParametersChange(
                    time=200000, growth_rate=0.64377516497364/(4*Ne)
                ),
                msprime.PopulationParametersChange(
                    time=300000, growth_rate=-0.0378691273513906/(4*Ne)
                ),
                msprime.PopulationParametersChange(
                    time=2000000, growth_rate=-0.064377516497364/(4*Ne)
                ),
                msprime.PopulationParametersChange(
                    time=3000000, growth_rate=0.00378691273513906/(4*Ne)
                ),
                msprime.PopulationParametersChange(
                    time=20000000, growth_rate=0, initial_size=Ne
                ),
        ]
    pop_conf=[msprime.PopulationConfiguration(initial_size=Ne_0,sample_size=2)]
    debugger = msprime.DemographyDebugger(Ne=Ne_0,
        population_configurations=pop_conf, migration_matrix=None,
        demographic_events=demographic_events, model='hudson')
    pop_size_traj = debugger.population_size_trajectory(time_segmentation)
    # add the timing
    pop_size = []
    for my_Ne, my_time in zip(pop_size_traj, time_segmentation):
        pop_size.append([my_Ne[0], my_time])

    # remove steps if no change in pop size
    pop_size_cleared = [pop_size[0]]
    for p,t in pop_size:
        if p != pop_size_cleared[-1][0]:
            pop_size_cleared.append([p, t])

    # debug print
    #print("_" * 80); print(pop_size_cleared); print("=" * 80, end="\n\n")

    return pop_size_cleared # alias pop size trajectory at times of time segmentation

for ts in t_sigma:
    for za in zigzag_amplitude:
        piecewise_Ne.append(zigzag_function(time_segmentation, za, Ne))
        piecewise_sigma.append([
            [sigma_0, 0], [sigma_1, ts*Ne]
            ])

# bottleneck
for ts in t_sigma:
    for tb in bottleneck_time:
        for sb in bottleneck_strength:
            piecewise_Ne.append([
                [Ne, 0], [Ne/sb, tb*Ne], [Ne, (tb+bottleneck_length*Ne)]
                ])
            piecewise_sigma.append([
                [sigma_0, 0], [sigma_1, ts*Ne]
                ])

# clean print for config file
def pretty_print_for_config_file(Ne, sigma):
    # round the times and num of individual to integers
    my_func = lambda x: int(round(x,0))
    Ne = [[[my_func(k) for k in j] for j in i] for i in Ne]
    # round the times to integers
    sigma = [[[my_func(k) if kx == 1 else k for kx, k in enumerate(j)] for j
        in i] for i in sigma]

    print("_"*80)
    
    print("population_sizes_backward_in_time:")
    for a in Ne:
        print(f" - {str(a)}")

    print("selfing_rates_backward_in_time:")
    for b in sigma:
        print(f" - {str(b)}")
    
    print("="*80, end="\n\n")

assert len(piecewise_Ne) == len(piecewise_sigma), "piecewise function for all demographies"
pretty_print_for_config_file(piecewise_Ne, piecewise_sigma)























