"""
Thist part of the pipeline shall make the estimates with the r-abc using the
previously prepared files from module 01
"""


rule abc_estimates:
    output:
        posteriors = "results/posteriors/statset_{sscomp}/pop_{obs}.pls_{pls}.RDS",
        posterior_plots = directory("results/posteriors/statset_{sscomp}/pop_{obs}.pls_{pls}")
    input:
        abc = "results/performance/pls_transformed/statset_{sscomp}/abc_sim.transformed.txt",
        obs = "results/performance/pls_transformed/statset_{sscomp}/params_and_sumstats.{obs}.table.txt"
    params:
        script = "workflow/scripts/abc.loclinear.R"
    shadow: "shallow"
    script: "../scripts/abc.loclinear.R"


