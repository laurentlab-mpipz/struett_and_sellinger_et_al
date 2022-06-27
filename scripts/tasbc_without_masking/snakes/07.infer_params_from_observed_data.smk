
"""
snake to infer the parameters of the model using the summary statistics from in
inferred tree sequence (stage 06) using pls and everyithing as in the pipeline
 e.g. when was the transition to selfing in CEU of A. thaliana
"""

localrules: obs_pls_flag

# expansion over pls
obs_PLS = [int(i) for i in config['pls_components_for_loclinear'] ]

# load helper functions
hf = print_utils(False)

rule obs_performance_analysis_of_neuralnetwork_regression:
    output: "data/inference/obs_set_{obs_iloc}/par_inference_{sscomp}_method_neuralnet.table.gzip"
    input:
        abc_stats = "data/table/params_and_sumstats.table.gzip",
        pod_stats = "data/observed/stats.{obs_iloc}.table.gzip"
    params:
        script = "scripts/performance_analysis_of_neuralnetwork_regression.R",
        sumstats_to_use = hf.get_string_from_yaml_configuration_file_for_Rscript(
            config['sumstat_combination_to_use']),
        num_simulations_to_tolerate = config['par_inf_number_of_tolerated_simulations'],
        path_to_config_file = "config/config.yaml"
    threads: 1
    shell:
        r"""
        Rscript --vanilla {params.script} {output} \
            {input.abc_stats} {input.pod_stats} \
            {wildcards.sscomp} \
            "{params.sumstats_to_use}" \
            {params.num_simulations_to_tolerate} \
            {params.path_to_config_file}
        """

###loclin#### this part is somewhat fragile to changing directories and similar
# this part profits a lot from already defined pls as in the performance analysis

rule provide_obs_table_as_txt:
    output: "data/table/txt/params_and_sumstats_obs_{obs_iloc}.table.txt"
    input: table = "data/observed/stats.{obs_iloc}.table.gzip"
    group: "pls_transform_stats"
    run:
        # use pandas to write the data frame as txt
        print("prepare.. writing table as txt")
        pd.read_pickle(input.table).to_csv(str(output), sep="\t", na_rep="NA", index=False, compression=None)
        print("done.. writing table as txt")

rule transform_sumstats_of_obs_sims:
    output:
        pod_sim_transformed = "data/inference/pls_transformed/{sscomp}/params_and_sumstats_obs_{obs_iloc}_sim.transformed.txt"
    input:
        pls = "data/performance/pls_transformed/{sscomp}/Routput_{sscomp}",
        pod_sim_raw = "data/table/txt/params_and_sumstats_obs_{obs_iloc}.table.txt"
    params:
        script = "scripts/transformer",
        script_subset = "scripts/subset_table_and_separate_params.R"
    shadow: "shallow"
    group: "pls_transform_stats"
    shell:
        r"""
        ## step 1: subset the table to transform by the columns identified for transformation and separate the params
        # last arg is the prefix for the output: subset.table[.param/.sumstat]
        Rscript --vanilla {params.script_subset} {input.pls} {input.pod_sim_raw} subset.table
        echo "_______________"
        echo "finished step 1"
        echo "==============="

        ## step 2: transform
        # output file prefix here is simply output.transformed
        {params.script} {input.pls} subset.table.sumstat {output.pod_sim_transformed} boxcox
        echo "_______________"
        echo "finished step 2"
        echo "==============="
        """

rule parameter_inference_of_loclinear_regression:
    output: temp("data/inference/obs_{obs_iloc}/pls_{pls}/perf_sscomp_{sscomp}_method_loclinear.single.gzip")
    input:
        transformed_abc_stats = "data/performance/pls_transformed/{sscomp}/abc_sim.transformed.txt",
        transformed_pod_stats = "data/inference/pls_transformed/{sscomp}/params_and_sumstats_obs_{obs_iloc}_sim.transformed.txt",
        pls_flag = "data/inference/pls_flag/{obs_iloc}_{sscomp}_{pls}.flag"
    params:
        script = "scripts/parameter_inference_using_loclinear_regression.R",
        num_simulations_to_tolerate = config['par_inf_number_of_tolerated_simulations'],
        path_to_config_file = "config/config.yaml" # the name of the prior definition has to match the definition in the r-function
    threads: 1
    shell:
        r"""
        # please leave threads 1, otherwise writing between cores seems to slow down the process
        export OPENBLAS_NUM_THREADS={threads} OMP_NUM_THREADS={threads} MKL_NUM_THREADS={threads}

        # do the point estimation and perfromance for all requested pls
        Rscript --vanilla {params.script} {output} \
            {input.transformed_abc_stats} {input.transformed_pod_stats} \
            {params.num_simulations_to_tolerate} \
            {params.path_to_config_file} \
            {wildcards.pls}
        """

rule obs_pls_flag:
    output: touch("data/inference/pls_flag/{obs_iloc}_{sscomp}_{pls}.flag")
    shell: r"""echo "provide {wildcards} flag file" """

rule obs_aggregate_performance_analysis_of_loclinear_regression_by_pls:
    output: "data/inference/obs_set_{obs_iloc}/par_inference_{sscomp}_method_loclinear.table.gzip"
    input:
        expand("data/inference/obs_{{obs_iloc}}/pls_{pls}/perf_sscomp_{{sscomp}}_method_loclinear.single.gzip",
            pls = PLS)
    run:
        # concatenate the pickles
        df = pd.concat([pd.read_pickle(f) for f in input])
        df.to_pickle(str(output))






