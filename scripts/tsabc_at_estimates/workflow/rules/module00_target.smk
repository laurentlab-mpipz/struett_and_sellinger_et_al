"""
Function for target rule to provide the list of target files
"""


import sys
from os import walk


def target_files(wilcards, verbose=True, return_expanded_wildcards_only=False):
    """
    create a list of target files
    """
    # expansion list of wildcards
    SSCOMP = list(range(1, len(config["sumstat_combination_to_use"]) + 1))
    PLS = config["pls_components_for_loclinear"]
    OBS = get_observation_index(path="resources/", pattern="stats.|.table.gzip")


    if return_expanded_wildcards_only:
        return SSCOMP, PLS, OBS
    else:
        output_list = []

        # transformed stats tables for abc sims and observations
        if False:
            output_list.extend(
                expand(
                    "results/performance/pls_transformed/statset_{sscomp}/abc_sim.transformed.txt",
                    sscomp=SSCOMP,
                )
            )
            output_list.extend(
                expand(
                    "results/performance/pls_transformed/statset_{sscomp}/params_and_sumstats.{obs}.table.txt",
                    sscomp=SSCOMP,
                    obs=OBS,
                )
            )

        # make estimates
        if False:
            output_list.extend(
                expand("results/posteriors/statset_{sscomp}/pop_{obs}.pls_{pls}.RDS",
                    sscomp=SSCOMP,
                    obs=OBS,
                    pls=PLS
                )
            )

        # make average posterior plots and main result table
        if True:
            output_list.extend(
                expand("results/figures/statset_{sscomp}/sumPlot_pop_{obs}.pls_{pls}.pdf",
                    sscomp=SSCOMP,
                    obs=OBS,
                    pls=PLS
                )
            )
            output_list.append("results/figures/aggregated_posteriors.pdf")
            output_list.append("results/main_table.txt")

        # print final list
        if verbose:
            print("_" * 80, file=sys.stderr)
            print("  You are asking for following output files:", file=sys.stderr)
            for file_index, file in enumerate(output_list, start=1):
                print(f"    {file_index}) ", file, sep="", file=sys.stderr)
            print("=" * 80, end="\n" * 2, file=sys.stderr)

        return output_list


def get_observation_index(path="resources/", pattern="stats.|.table.gzip"):
    """Find the index of the files providing the observed statistics

    The files must be of the given pattern in the given path.

    Args:
        path: String, path to the folder containting the stats files
        pattern: String, must be of following format: "pre|post". Pre is the
            string preceding the index, post is the string postceding the index.

    Returns:
        A list with the indexes.
    """
    # list the files in the given path
    filenames = next(walk(path), (None, None, []))[2]  # [] if no file

    # prepare the pattern of the wanted stats files
    pre, post = pattern.split("|")

    # split until getting the index
    indices = [
        file.split(pre)[1].split(post)[0]
        for file in filenames
        if file.startswith(pre) and file.endswith(post)
    ]

    return indices
