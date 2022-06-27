
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

here: we read the config file of the ABC and create the lists as needed for
the teSMC pipeline
"""

import yaml
import sys

path_to_config_of_abc = "resources/ABC_config/config.yaml"

with open(path_to_config_of_abc, "r") as file:
    config = yaml.load(file, Loader=yaml.FullLoader)

N_r = config["pod_population_size_recent"]
N_a = config["pod_population_size_ancient"]
t_N = config["pod_population_size_change_time"]
s_r = config["pod_selfing_rate_recent"]
s_a = config["pod_selfing_rate_ancient"]
t_s = config["pod_selfing_rate_change_time"]

# transform the model to piecewise constant per time
piecewise_Ne = []
piecewise_sigma = []
for i, (nr, na, tN, sr, sa, ts) in enumerate(zip(N_r, N_a, t_N, s_r, s_a, t_s)):
    # filter for wanted conditions
    if not (
        nr == 40_000 and
        ts <= 300_000
        ):
        continue

    piecewise_Ne.append([[nr, 0], [na, tN]])
    piecewise_sigma.append([[sr, 0], [sa, ts]])

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



