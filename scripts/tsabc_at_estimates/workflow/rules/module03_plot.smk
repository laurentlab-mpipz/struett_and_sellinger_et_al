"""
plot the posteriors in a useful way
"""


# get expanded wildcards from module00_target.smk
SSCOMP, PLS, OBS =  target_files(None, verbose=False,
    return_expanded_wildcards_only=True)


rule abc_summarizing_plots:
    output:
        posterior_plot = "results/figures/statset_{sscomp}/sumPlot_pop_{obs}.pls_{pls}.pdf"
    input:
        posteriors = "results/posteriors/statset_{sscomp}/pop_{obs}.pls_{pls}.RDS"
    script: "../scripts/abc.plot.posterior.R"


rule plot_aggregated_estimates:
    output:
        agg_sumplot = "results/figures/aggregated_posteriors.pdf",
        table = "results/main_table.txt"
    input:
        agg_post = expand("results/posteriors/statset_{sscomp}/pop_{obs}.pls_{pls}.RDS",
            sscomp=SSCOMP,
            obs=OBS,
            pls=PLS)
    script: "../scripts/abc.plot.aggregated.posteriors.R"
