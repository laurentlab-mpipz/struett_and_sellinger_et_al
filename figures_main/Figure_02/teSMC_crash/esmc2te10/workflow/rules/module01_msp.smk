
"""
simulate tree_sequences using msprime
"""

def demography2int(wildcards):
    return int(wildcards.demography)

rule simulate_tree_sequence_using_msprime:
    output:
        ts = temp("results/tree_sequences/demographic_scenario_{demography}/"
            "rep_{replication}.nsam_{sample_size}.mu_{mutation_rate}."
            "r_{recombination_rate}.length_{chromosome_length}.ts_{chr_number}")
    params:
        pop_size_fn = lambda wildcards: config[
            "population_sizes_backward_in_time"][int(float(
                wildcards.demography))],
        sigma_fn = lambda wildcards: config["selfing_rates_backward_in_time"][
            int(float(wildcards.demography))],
        nsam = lambda wildcards: int(float(wildcards.sample_size)),
        mu = lambda wildcards: float(wildcards.mutation_rate),
        r = lambda wildcards: float(wildcards.recombination_rate),
        L = lambda wildcards: int(float(wildcards.chromosome_length))
    threads: 1
    run:
        # parameter preparation
        pop_sizes, pop_size_times = pop_size_over_time(params.pop_size_fn)
        recombination_rates, recombination_times = (
            r_over_time_from_sigma_over_time(params.r, params.sigma_fn))
        rescaled_pop_sizes, rescaled_pop_size_times = rescale_pop_size_by_sigma(
            params.pop_size_fn, params.sigma_fn)

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
        ts = simulate_change_in_recombination(simulation_parameters, False)
        print("Number of trees: {}".format(ts.num_trees))

        # mutate trees
        print(f"mu {params.mu}")
        mutated = msprime.mutate(ts, rate=params.mu, random_seed=random.randint(
            1,2**32-1), model=None, keep=True)
        print("Number of mutations: {}".format(mutated.num_mutations))

        mutated.dump(output.ts)
                
