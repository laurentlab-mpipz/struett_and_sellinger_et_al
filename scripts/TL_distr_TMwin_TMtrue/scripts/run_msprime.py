import argparse
import msprime
import tskit
import random

parser = argparse.ArgumentParser(description=
    "Run msprime with a change in recombination_rate")
parser.add_argument("-o", "--output_file", type=str, required=True,
    help="path to output file (trees file)")
parser.add_argument("-sigma", "--selfing_rate", type=float, required=True,
    help="rate of selfing to change to")
parser.add_argument("-t_sigma", "--t_change_in_selfing_rate", type=float, required=True,
    help="timing of transition from outcrossing to selfing")
parser.add_argument("-Ne", "--population_size", type=int, required=True,
    help="constant population size (equal to effective pop size under outcrossing)")
parser.add_argument("-r", "--recombination_rate", type=float, required=True,
    help="recombination_rate (effective rate under outcrossing)")
parser.add_argument("-n", "--sample_size", type=int, required=True,
    help="sample size of simulated trees file")
parser.add_argument("-L", "--chromosome_length", required=True,
    help="length of simulated chromosome (uniform recombination rate)")

args = parser.parse_args()

# parameters from cmd line
outfile = args.output_file
selfing_rate = args.selfing_rate
t_sigma = args.t_change_in_selfing_rate
population_size = args.population_size
recombination_rate = args.recombination_rate
nsam = args.sample_size
L = int(float(args.chromosome_length))

################################################################################
# functions
def getNandRfromSigma(Ne_phase_2, recombination_rate_phase_2, sigma):
    Fis = sigma / ( 2 - sigma )
    Ne_phase_1 = Ne_phase_2 * 1 / ( 1 + Fis )
    recombination_rate_phase_1 = recombination_rate_phase_2 * ( 1 - Fis )
    return Ne_phase_1, recombination_rate_phase_1

def simulate_change_in_recombination(simulation_parameters: dict, verbose=False):

    if verbose: print("")

    # prune demographic events
    events_phase_1 = [event for event in
        simulation_parameters["demographic_events"] if
        event.time < simulation_parameters["t_recombination_change"]]
    events_phase_2 = [event for event in
        simulation_parameters["demographic_events"] if
        event.time >= simulation_parameters["t_recombination_change"]]

    # set events for phase 1
    simulation_parameters["demographic_events"] = events_phase_1

    if verbose: print("\tDemographic events configured for phase 1")

    # Phase 1 simulation
    ts_phase_1 = msprime.simulate(
        sample_size=simulation_parameters["sample_size"],
        Ne=simulation_parameters["Ne_phase_1"],
        length=simulation_parameters["length"],
        recombination_rate=simulation_parameters["recombination_rate_phase_1"],
        population_configurations=
            simulation_parameters["population_configurations"],
        demographic_events=simulation_parameters["demographic_events"],
        random_seed=simulation_parameters["random_seed_phase_1"],
        end_time=simulation_parameters["t_recombination_change"],
        model=simulation_parameters["model"]
    )

    if verbose: print("\tPhase 1 finished")

    # set events for phase 1
    simulation_parameters["demographic_events"] = events_phase_2

    if verbose: print("\tDemographic events configured for phase 2")

    # it occured a few times, that the end_time of the trees was slightly
    # different from the new start_time (e.g. 0.00000000006 difference)
    # here I test, that they will not differ more than 1 generation
    assert max(ts_phase_1.tables.nodes.time) - simulation_parameters["t_recombination_change"] < 1, "Problems with the transition {}, {}".format(
        simulation_parameters["t_recombination_change"],
        max(ts_phase_1.tables.nodes.time))

    # Phase 2
    ts_phase_2 = msprime.simulate(
        Ne=simulation_parameters["Ne_phase_2"],
        recombination_rate=simulation_parameters["recombination_rate_phase_2"],
        demographic_events=simulation_parameters["demographic_events"],
        random_seed=simulation_parameters["random_seed_phase_2"],
        from_ts=ts_phase_1,
#        start_time=simulation_parameters["t_recombination_change"]
    )

    if verbose: print("\tPhase 2 finished")
    if verbose: print("")

    if verbose: print("num trees\n\tphase 1: {}\n\tphase 2: {}\n".format(
        ts_phase_1.num_trees, ts_phase_2.num_trees))

    return ts_phase_2.simplify()

################################################################################
# script

# set seed manually
set_seed = False
manual_seed = 1

# set seed
seed = manual_seed if set_seed else random.randint(1,2**32-1)
random.seed(seed)

# translate selfing into the two phase model
sigma = selfing_rate
Ne = population_size
r = recombination_rate

Ne_phase_1, recombination_rate_phase_1 = getNandRfromSigma(Ne, r, sigma)

# define simulation parameters in dictionary
simulation_parameters = {
    "sample_size": nsam,
    "Ne_phase_1": Ne_phase_1,
    "Ne_phase_2": Ne,
    "length": L,
    "recombination_rate_phase_1": recombination_rate_phase_1,
    "recombination_rate_phase_2": r,
    "population_configurations": None,
    "demographic_events": [
        msprime.SimulationModelChange(1000, "hudson"),
        msprime.SimulationModelChange(1e6, "smc_prime")
    ],
    "random_seed_phase_1": random.randint(1,2**32-1),
    "random_seed_phase_2": random.randint(1,2**32-1),
    "t_recombination_change": int(t_sigma*Ne),
    "model": "dtwf"
}

# run the simulation
ts = simulate_change_in_recombination(simulation_parameters)

# test output
print("Number of trees: {}".format(ts.num_trees))

# save trees file
ts.dump(outfile)
