
"""
example how to run varying selfing rates / or recombination rates; we use
Norcborg's approximation whithout the refinment on biallelic loci
"""

import tskit
import msprime
from mspts import *
import random

def main(Ne,nsim,self,sc,data):
    class params:
        data=data+1

        if data == 1:
           L=5*1e6

        if data == 1:
           i=1

        r =5*1e-8
        nsam = 10
        mu =1e-8

        if sc == 1:
           pop_size_fn = [[Ne, 0]]

        if sc == 2:
           pop_size_fn = [[5*Ne, 0],[Ne,0.125*4*Ne],[5*Ne,0.5*4*Ne]]

        if sc == 3:
           pop_size_fn = [[5*Ne, 0],[Ne,2*Ne]]

        if sc == 4:
           pop_size_fn = [[Ne, 0],[5*Ne,2*Ne]]

        if self == 1:
           sigma_fn = [[0.95, 0],[0, 0.5*Ne]] 

        if self == 2:
           sigma_fn = [[0.95, 0],[0, Ne]] 

        if self == 3:
           sigma_fn = [[0.95, 0],[0, 2*Ne]] 

        if self == 4:
           sigma_fn = [[0.95, 0],[0, 4*Ne]] 

    class output:
        ts = f"example.ts"

    # parameter preparation
        pop_sizes, pop_size_times = pop_size_over_time(params.pop_size_fn)
        recombination_rates, recombination_times = (r_over_time_from_sigma_over_time(params.r, params.sigma_fn))
        rescaled_pop_sizes, rescaled_pop_size_times = rescale_pop_size_by_sigma(params.pop_size_fn, params.sigma_fn)

    # simulation parameters
        simulation_parameters = {
            "sample_size": params.nsam,
            "pop_size_over_time": rescaled_pop_sizes,
            "pop_size_times": rescaled_pop_size_times,
            "length": params.L,
            "recombination_rate_over_time": recombination_rates,
            "recombination_rate_times": recombination_times,
            "population_configurations": None,
            "demographic_events": [
        	msprime.SimulationModelChange(1000, "hudson"),
        	msprime.SimulationModelChange(1e6, "smc_prime")
            ],
            "model": "dtwf"# starting model
        }

	# simulate trees
        for x in range(nsim):
            for chr in range(params.i):
                ts = simulate_change_in_recombination(simulation_parameters, False)
                ts=ts.simplify()
                ts = msprime.mutate(ts, rate=params.mu, random_seed=random.randint(1,2**32-1), model=None, keep=True)
                f = open("Figure_Setfan"+str(params.mu)+str(params.r)+"x"+str(x)+"chr"+str(chr)+"self"+str(self)+"sc"+str(sc)+"data"+str(params.data)+".txt","w")
                f.write("num sites"+str(ts.get_num_sites())+'\n')
                for variant in ts.variants():
                    f.write(str(variant.site.position)+str(variant.genotypes)+'\n')
                f.close()
data_v=[0]
for data in data_v:
    self_v=[1,2,3,4]
    for self in self_v:
        sc_v=[2]
        for sc in sc_v:
            main(10000,10,self,sc,data)
