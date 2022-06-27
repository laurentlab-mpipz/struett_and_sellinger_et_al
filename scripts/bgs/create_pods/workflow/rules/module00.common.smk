import os
import numpy as np
import pandas as pd


def all_output(wildcards, verbose=False, more_stats=True):
    """list of output files in final target (overall) rule

    Args:
        wildcards: wildcards of input rule

    Returns:
        list of files to be created in this snake
    """
    num_pod_simulations = int(float(config["pod_number_of_simulations_per"]))
    pod_param_setup = [
        config["pod_population_size_recent"],
        config["pod_population_size_ancient"],
        config["pod_population_size_change_time"],
        config["pod_selfing_rate_recent"],
        config["pod_selfing_rate_ancient"],
        config["pod_selfing_rate_change_time"],
    ]
    nums_pods = [len(i) for i in pod_param_setup]
    assert len(set(nums_pods)) == 1, "definitions of pod parameters vary in size"
    id_pod = list(range(num_pod_simulations))
    pod_iloc = list(range(len(config["pod_selfing_rate_change_time"])))

    if verbose:
        print("_" * 80)
        print("\twildcards")
        print(wildcards)
        print("=" * 80, end="\n" * 2)

    output_list = []

    # burnin files
    #    output_list.extend(expand("results/burnin/pod_{pod_no}.burnin.ts",
    #        pod_no=range(num_pod_simulations)))

    # transition files
    output_list.extend(
        expand(
            "results/bgs/pod_{pod_no}/tsigma_{tsigma}_locusno_{locus}.ts",
            locus=range(config["num_independent_regions"]),
            tsigma=tsigma_slim_asList,
            pod_no=range(num_pod_simulations),
        )
    )

    # collected stats
    output_list.extend(
        expand(
            "results/stats/params_and_sumstats_pod_{pod_iloc}.table.gzip",
            pod_iloc=pod_iloc,
        )
    )

    # measured sumstats over time
    if (
        more_stats
    ):  # the pipeline is not created to produce these files; I recommend to run this part only when all the work is done
        print("You set to calculate more of the sumstats. Please revise!")
        output_list.extend(
            expand("results/more_stats/tsigma_{tsigma}.gzip", tsigma=grep_all_tsigma())
        )
        output_list.append(
            "results/more_stats/more_stats.collected.csv"
        )  # collect also the differen tsigma
        output_list.append(
            "results/more_stats_fpop/more_stats.collected.csv"
        )  # if we do not downsample

    # mhs for teSMC
    output_list.extend(
        expand("results/mhs/pod_{pod_no}_tsigma_{tsigma}.mhs",
            pod_no=range(num_pod_simulations),
            tsigma=grep_all_tsigma())
        )

    return output_list


def sumstat_on_pyslim_input_func(wc):
    """input for rule sumstat_on_pyslim

    Provide the different loci from different burn-ins. That is very crucial
    to mimic different chromosomes. We use a simple algorithm to reshuffle the
    loci from different burn-ins. For pod_no of output, we take the locus 1 of
    the same pod_no as input. For locus 2 we go one pod_no up and so forth. We
    just overflow the number of pods back to zero, once we reached the maximal
    number of pods.
    """
    file_dependecies = []

    # output from previous rule
    # "results/bgs/pod_{pod_no}/tsigma_{tsigma}_locusno_{locus}.ts"

    # all possible pod_no and locus wildcard values
    pods = list(range(num_pod_simulations))
    loci = list(range(config["num_independent_regions"]))

    for p in [int(float(wc.pod_no))]:
        this_p = p
        for l in loci:
            file_dependecies.append(
                f"results/bgs/pod_{this_p}/tsigma_{wc.tsigma}_locusno_{l}.ts"
            )
            this_p = (this_p + 1) % len(pods)

    return file_dependecies


def random_sample_from_treeseq(my_ts, sample_size):
    """Get a random subsample from the provided tree sequence

    Args:
        my_ts: a tree sequence
        sample_size: the number of samples to choose randomly from

    Returns:
        simplified tree sequence
    """
    sample_nodes = [i for i in my_ts.samples()]
    chosen_samples = np.random.choice(sample_nodes, sample_size, replace=False)
    simpel_ts = my_ts.simplify(samples=chosen_samples)
    return simpel_ts


def sample_1perInd(my_ts, seed=None):
    """Simplify to one haplotype per individual

    Intended to use on pyslim treeseq, we sample only 1 single
    haplotype per individual.

    Args:
        my_ts: tree sequence from pyslim

    Returns:
        simplified tree sequence with a single haplotype per individual
    """
    np.random.seed(seed)

    list_of_1hap_samples = []
    for i in my_ts.individuals():
        list_of_1hap_samples.append(np.random.choice(i.nodes))
    list_of_1hap_samples = np.array(list_of_1hap_samples)

    my_simple_ts = my_ts.simplify(samples=list_of_1hap_samples)

    return my_simple_ts


def tsigma_from_pod_iloc(wc):
    """get tsigma by pod index

    We use the config to find the index of a pod based on the given tsigma. We
    do so, to keep the file nameing consistent with the previous pipeline, which
    should allow us to run model choice as well as the performance analysis in
    parameter estimation

    Args:
        wc: wildcards object from snakemake rule

    Returns:
        The tsigma based on the snakemake config object by the index.
    """
    return config["pod_selfing_rate_change_time"][int(float(wc.pod_iloc))]


def grep_all_tsigma():
    """find all tsigma to expand over

    This is a hard coded function. The only purpose is to return all possible
    tsigma for the created tree sequences to loop over. This is done by finding
    all possible values based on a file string pattern.
    """
    tsigma_list = []

    rootdir = "results/bgs/"

    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            complete_path = os.path.join(subdir, file)
            if "tsigma" in complete_path:
                this_tsigma = int(complete_path.split("tsigma_")[1].split("_")[0])
                tsigma_list.append(this_tsigma)
    tsigma_list = sorted(list(set(tsigma_list)))

    return tsigma_list


def more_stats(tree_seq, params, log):
    """Calculate vector of stats from tree sequence.

    Args:
        tree_seq: pyslim or tskit tree sequence

    Returns:
        named vector of summary statistics; type(pd.DataFrame)
    """
    # recaptitate; only for the use of the summarizing function; no mutating
    recap = tree_seq.recapitate(
        recombination_rate=params.recombrate, Ne=params.pop_size_outcrossing
    )

    ts_1hapInd = sample_1perInd(recap)  # sample a single haplotype per individual

    # subsample
    ts_1perInd_nsam = random_sample_from_treeseq(ts_1hapInd, sample_size=params.nsam)

    # some tree stats
    with open(log.std, "w") as loge:
        print("_" * 80, file=loge)
        print(ts_1perInd_nsam, file=loge)
        print("=" * 80, end="\n" * 2, file=loge)

    stats = np.array(
        [
            ts_1perInd_nsam.Tajimas_D(mode="site"),
            ts_1perInd_nsam.diversity(mode="site"),
            ts_1perInd_nsam.segregating_sites(mode="site", span_normalise=False),
        ]
    )

    colnames = ["Tajimas_D", "diversity", "segregating_sites"]

    # multi-dim stats; 1d-sfs
    stats_sfs = ts_1perInd_nsam.allele_frequency_spectrum(
        span_normalise=False, polarised=True
    )
    stats = np.append(stats, stats_sfs)  # append sfs to stats list
    _ = [
        colnames.append("sfs_" + str(i)) for i, _ in enumerate(stats_sfs)
    ]  # append colnames for sfs

    df_stats = pd.DataFrame([stats], columns=colnames)

    return df_stats


def more_stats_fpop(tree_seq, params, log):
    """Calculate vector of stats from tree sequence.

    Different from the more_stats function, here we calculate all the stats on
    the whole population (1 haplotype per individual). We do not downsample to
    a given pop size or something.

    Args:
        tree_seq: pyslim or tskit tree sequence

    Returns:
        named vector of summary statistics; type(pd.DataFrame)
    """
    # recaptitate; only for the use of the summarizing function; no mutating
    recap = tree_seq.recapitate(
        recombination_rate=params.recombrate, Ne=params.pop_size_outcrossing
    )

    ts_1hapInd = sample_1perInd(recap)  # sample a single haplotype per individual

    # some tree stats
    with open(log.std, "w") as loge:
        print("_" * 80, file=loge)
        print(ts_1hapInd, file=loge)
        print("=" * 80, end="\n" * 2, file=loge)

    stats = np.array(
        [
            ts_1hapInd.Tajimas_D(mode="site"),
            ts_1hapInd.diversity(mode="site"),
            ts_1hapInd.segregating_sites(mode="site", span_normalise=False),
        ]
    )

    colnames = ["Tajimas_D", "diversity", "segregating_sites"]

    # multi-dim stats; 1d-sfs
    stats_sfs = ts_1hapInd.allele_frequency_spectrum(
        span_normalise=False, polarised=True
    )
    stats = np.append(stats, stats_sfs)  # append sfs to stats list
    _ = [
        colnames.append("sfs_" + str(i)) for i, _ in enumerate(stats_sfs)
    ]  # append colnames for sfs

    df_stats = pd.DataFrame([stats], columns=colnames)

    return df_stats

