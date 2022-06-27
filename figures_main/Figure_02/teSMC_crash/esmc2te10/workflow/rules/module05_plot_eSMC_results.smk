
rule plot_eSMC_results:
    output:
        pdf = "results/figures/esmc_full_results.pdf"
    input:
        rds =  "results/esmc/table.rds"
    script: "../scripts/plot_esmc_results.R"

rule plot_mse_eSMC:
    output:
        pdf = "results/figures/mse_full.pdf"
    input:
        df = "results/mse/table.csv"
    script: "../scripts/plot_esmc_mse_results.R"