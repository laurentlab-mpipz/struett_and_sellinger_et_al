
"""
estimate piecewise constant selfing rates
"""

# localrules: plot_true_vs_estimated_teSMC, plot_aggregated_estimations_teSMC
localrules: aggregate_inference_models

rule teSMC_inferenceModel:
    output:
        full_data = "results/esmc/demographic_scenario_{demography}/"
            "rep_{replication}.nsam_{sample_size}.mu_{mutation_rate}."
            "r_{recombination_rate}.length_{chromosome_length}.{inference_model}.RDS"
    input:
        mhs = "results/mhs/demographic_scenario_{demography}/"
            "rep_{replication}.nsam_{sample_size}.mu_{mutation_rate}."
            "r_{recombination_rate}.length_{chromosome_length}.mhs",
        esmc_package = "resources/esmc/eSMC2_4.0.3.tar.gz"
    log:
        esmc_log = "results/esmc/demographic_scenario_{demography}/"
            "rep_{replication}.nsam_{sample_size}.mu_{mutation_rate}."
            "r_{recombination_rate}.length_{chromosome_length}.{inference_model}.log"
    params:
        inference_model = "{inference_model}",
        number_of_hidden_states = 40,
        recent_selfing_rate = 0.99,
        ancient_selfing_rate = 0.1,
    resources:
        swap_gb = 14
    threads: 1
    priority: 50
    shadow: "shallow"
    group: "esmc2"
    script: "../scripts/model_esmc.R"

rule aggregate_inference_models:
    output:
        csv = "results/esmc/demographic_scenario_{demography}/"
            "rep_{replication}.nsam_{sample_size}.mu_{mutation_rate}."
            "r_{recombination_rate}.length_{chromosome_length}.csv"
    input:
        full = expand("results/esmc/demographic_scenario_{{demography}}/"
            "rep_{{replication}}.nsam_{{sample_size}}.mu_{{mutation_rate}}."
            "r_{{recombination_rate}}.length_{{chromosome_length}}.{inference_model}.RDS",
            inference_model = config["teSMC"]["inference_models"])
    script: "../scripts/model_esmc_agg.R"

rule aggregate_all_inferences:
    output:
        rds =  "results/esmc/table.rds" # rds to save disk space 
    input:
        expand("results/esmc/demographic_scenario_{demography}/"
            "rep_{replication}.nsam_{sample_size}.mu_{mutation_rate}."
            "r_{recombination_rate}.length_{chromosome_length}.csv",
            demography=range(len(config["population_sizes_backward_in_time"])),
            replication=range(config["technical_replicates"]),
            sample_size=config["sample_size"],
            mutation_rate=config["mutation_rate"],
            recombination_rate=config["recombination_rate"],
            chromosome_length=config["chromosome_length"]
        )
    script: "../scripts/full_esmc_agg.R"

''' old script, what I used when arriving to Muc
rule estimate_with_teSMC:
    output:
        data = "results/esmc/demographic_scenario_{demography}/"
            "rep_{replication}.nsam_{sample_size}.mu_{mutation_rate}."
            "r_{recombination_rate}.length_{chromosome_length}.RDS"
    input:
        mhs = "results/mhs/demographic_scenario_{demography}/"
            "rep_{replication}.nsam_{sample_size}.mu_{mutation_rate}."
            "r_{recombination_rate}.length_{chromosome_length}.mhs"
    log:
        esmc_log = "results/esmc/demographic_scenario_{demography}/"
            "rep_{replication}.nsam_{sample_size}.mu_{mutation_rate}."
            "r_{recombination_rate}.length_{chromosome_length}.log"
    params:
        esmc_package = "resources/esmc/source_teSMC_sigmaEstim_bw.R"
    threads: 1
    script: "../scripts/esmc.R"

rule plot_true_vs_estimated_teSMC:
    output:
        pdf = "results/esmc/demographic_scenario_{demography}/"
            "rep_{replication}.nsam_{sample_size}.mu_{mutation_rate}."
            "r_{recombination_rate}.length_{chromosome_length}.pdf"
    input:
        data = "results/esmc/demographic_scenario_{demography}/"
            "rep_{replication}.nsam_{sample_size}.mu_{mutation_rate}."
            "r_{recombination_rate}.length_{chromosome_length}.RDS"
    params:
        path_to_config = path_to_config
    priority: 50
    threads: 1
    script: "../scripts/plot_true_vs_estimate.R"

rule plot_aggregated_estimations_teSMC:
    output:
        pdf = "results/esmc/agg_plots/demographic_scenario_{demography}."
            "aggregated.nsam_{sample_size}.mu_{mutation_rate}."
            "r_{recombination_rate}.length_{chromosome_length}.pdf"
    input:
        data = expand("results/esmc/demographic_scenario_{{demography}}/"
            "rep_{replication}.nsam_{{sample_size}}.mu_{{mutation_rate}}."
            "r_{{recombination_rate}}.length_{{chromosome_length}}.RDS",
            replication=range(config["technical_replicates"]))
    params:
        path_to_config = path_to_config
    priority: 60
    threads: 1
    script: "../scripts/plot_aggregated_true_vs_estimate.R"
'''
