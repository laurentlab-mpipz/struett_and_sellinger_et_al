"""
Thist part of the pipeline shall prepare the tables into the form that the abc
can directly be applied

"""


hf = print_utils(False)


rule get_pls_dim_reduced_sumstats_for_abc_sims:
    output:
        rmse_plot = "results/performance/pls_transformed/statset_{sscomp}/RMSE_{sscomp}.pdf",
        transformer = "results/performance/pls_transformed/statset_{sscomp}/Routput_{sscomp}"
    input:
        "resources/params_and_sumstats.table.1to135909.gzip"
    params:
        num_pls_max = config['number_max_pls_for_transformation'],
        sumstats_to_use = hf.get_string_from_yaml_configuration_file_for_Rscript(
            config['sumstat_combination_to_use']),
        script = "workflow/scripts/find_pls_by_wanted_sumstats.R",
    threads: 1
    shell:
        r"""
        # please leave threads 1, otherwise writing between cores seems to slow down the process
        export OPENBLAS_NUM_THREADS={threads} OMP_NUM_THREADS={threads} MKL_NUM_THREADS={threads}
        Rscript --vanilla {params.script} \
            {output.transformer} \
            {output.rmse_plot} \
            {input} \
            {params.num_pls_max} \
            {wildcards.sscomp} \
            "{params.sumstats_to_use}"
        """


rule provide_abc_table_as_txt:
    output: "results/table/txt/params_and_sumstats.table.txt"
    input:
    	table = "resources/params_and_sumstats.table.1to135909.gzip"
    run:
        # use pandas to write the data frame as txt
        print("prepare.. writing table as txt")
        pd.read_pickle(input.table, compression="gzip").to_csv(
            str(output), sep="\t", na_rep="NA", index=False)
        print("done.. writing table as txt")


rule provide_pod_table_as_txt:
    output: "results/table/txt/params_and_sumstats.{obs}.table.txt"
    input:
    	table = "resources/stats.{obs}.table.gzip"
    group: "pls_transform_stats"
    run:
        # use pandas to write the data frame as txt
        print("prepare.. writing table as txt")
        pd.read_pickle(input.table).to_csv(str(output), sep="\t", na_rep="NA", index=False)
        print("done.. writing table as txt")


rule transform_sumstats_of_abc_sims:
    output:
        abc_sim_transformed = "results/performance/pls_transformed/statset_{sscomp}/abc_sim.transformed.txt"
    input:
        pls = "results/performance/pls_transformed/statset_{sscomp}/Routput_{sscomp}",
        abc_sim_raw = "results/table/txt/params_and_sumstats.table.txt"
    params:
        script = "workflow/scripts/transformer",
        script_subset = "workflow/scripts/subset_table_and_separate_params.py"
    shadow: "shallow"
    threads: 1
    shell:
        r"""
        ## step 1: subset the table to transform by the columns identified for transformation and separate the params
        # last arg is the prefix for the output: subset.table[.param/.sumstat]
        python {params.script_subset} {input.pls} {input.abc_sim_raw} subset.table
        echo "_______________"
        echo "finished step 1"
        echo "==============="

        ## step 2: transform
        # output file prefix here is simply output.transformed
        {params.script} {input.pls} subset.table.sumstat {output.abc_sim_transformed} boxcox
        echo "_______________"
        echo "finished step 2"
        echo "==============="
        """


rule transform_sumstats_of_pod_sims:
    output:
        pod_sim_transformed = "results/performance/pls_transformed/statset_{sscomp}/params_and_sumstats.{obs}.table.txt"
    input:
        pls = "results/performance/pls_transformed/statset_{sscomp}/Routput_{sscomp}",
        pod_sim_raw = "results/table/txt/params_and_sumstats.{obs}.table.txt"
    params:
        script = "workflow/scripts/transformer",
        script_subset = "workflow/scripts/subset_table_and_separate_params.py"
    shadow: "shallow"
    group: "pls_transform_stats"
    shell:
        r"""
        ## step 1: subset the table to transform by the columns identified for transformation and separate the params
        # last arg is the prefix for the output: subset.table[.param/.sumstat]
        #Rscript --vanilla {params.script_subset} {input.pls} {input.pod_sim_raw} subset.table
        python {params.script_subset} {input.pls} {input.pod_sim_raw} subset.table
        head -n 6 subset.table* | cut -f 1-10
        echo "..."
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
