
# helper functions for running mspts

import msprime
import tskit
import numpy as np
import sys
import warnings

def simulate_change_in_recombination(simulation_parameters, verbose=True):
    # collection of demographic events
    events = []

    events.extend(simulation_parameters["demographic_events"])
    if verbose: print("\tadded explicit demographic events")

    # add pop size changes
    for t, p in zip(simulation_parameters["pop_size_times"],
        simulation_parameters["pop_size_over_time"]):
        events.append(msprime.PopulationParametersChange(t, initial_size=p,
            growth_rate=None, population=None, population_id=None))
    if verbose: print("\tadded pop size changes as demographic events")

    # sort events by time
    events.sort(key=lambda x: x.time)
    if verbose: print("\tsorted demographic events")

    if verbose:
        print("_" * 80)
        print("all events")
        for e in events:
            print(e.time, end="\t")
            print(e)
        print("=" * 80, end="\t\t")
        print()

    # loop through all phases
    ts = None
    for i, (t, r) in enumerate(zip(simulation_parameters[
        "recombination_rate_times"], simulation_parameters[
        "recombination_rate_over_time"])):
        # define end_time of this phase
        try:
            end_time = simulation_parameters["recombination_rate_times"][i+1]
        except:
            end_time = np.inf

        phase_events = [e for e in events if end_time > e.time >= t]

        if verbose:
            print("_" * 80)
            print("phase events")
            for e in phase_events:
                print(e.time, end="\t")
                print(e)
            print("=" * 80, end="\t\t")
            print()

        if verbose:
            if not ts is None:
                print(f"max root time: {ts.max_root_time} ", end="")
                print(f"compared to wanted t = {t}")

        # assert floating time issues in phase events
        if not ts is None and ts.max_root_time != t:
            if verbose: print("+"*80)
            assert ts.max_root_time > t, "unexpected error root time error"
            if not np.isclose(ts.max_root_time, t):
                warings.warn(" root time and phase time differs", Warning)
            for e_index, e in enumerate(phase_events):
                if verbose: print("phase event:", e)
                if np.isclose(e.time, ts.max_root_time):
                    e.time = ts.max_root_time
                    phase_events[e_index] = e
                elif e.time < ts.max_root_time:
                    e.time = ts.max_root_time
                    phase_events[e_index] = e
                    warnings.warn(" demographie might be affected")
                if verbose: print("phase event :", e)

        if verbose:
            if not ts is None:
                print(f"corrected max root time: {ts.max_root_time} ", end="")
                print(f"compared to wanted t = {phase_events[0].time}")
                print("+"*80)

        if i == 0:
            ts = msprime.simulate(
                sample_size=simulation_parameters["sample_size"],
                Ne=simulation_parameters["pop_size_over_time"][0],
                length=simulation_parameters["length"],
                recombination_rate=r,
                demographic_events=phase_events,
                end_time=end_time,
                model=simulation_parameters["model"],
                from_ts=ts
                )
        else:
            ts = msprime.simulate(
                length=simulation_parameters["length"],
                recombination_rate=r,
                demographic_events=phase_events,
                end_time=end_time,
                from_ts=ts
                )
        if verbose: print(f"{i} {ts}")

    return ts

def pop_size_over_time(pop_size_fn):
    """get pop_sizes, pop_size_times from config file
    """
    pop_sizes = []
    pop_times = []
    for p, t in pop_size_fn:
        pop_sizes.append(p)
        pop_times.append(t)

    return pop_sizes, pop_times

def r_over_time_from_sigma_over_time(r, sigma_fn):
    """get rec rates, rec rate times from config file
    """
    sigmas = []
    rec_times = []
    for s, t in sigma_fn:
        sigmas.append(s)
        rec_times.append(t)

    rec_rates = []
    for s in sigmas:
        fis = s/(2-s)
        rec_rates.append((1-fis)*r)

    return rec_rates, rec_times

def rescale_pop_size_by_sigma(pop_size_fn, sigma_fn):
    pop_sizes = []
    pop_size_times = []
    for p, t in pop_size_fn:
        pop_sizes.append(p)
        pop_size_times.append(t)

    sigmas = []
    sigma_times = []
    for s, t in sigma_fn:
        sigmas.append(s)
        sigma_times.append(t)

    # create rescaled pop sizes
    pop_sizes_rescaled = []
    pop_size_times_rescaled = sorted(list(set(pop_size_times + sigma_times)))
    for t in pop_size_times_rescaled:
        current_pop_size_unscaled = pop_sizes[np.where(np.array(pop_size_times)
            <= t)[0].max()]
        current_sigma = sigmas[np.where(np.array(sigma_times) <= t)[0].max()]

        pop_sizes_rescaled.append(current_pop_size_unscaled * (1 - 0.5 *
            current_sigma))

    return pop_sizes_rescaled, pop_size_times_rescaled
