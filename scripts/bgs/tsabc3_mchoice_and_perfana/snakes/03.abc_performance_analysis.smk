
"""
snake to analyse the performance of the abc using an abc parameter/sumstat-
table and pods
"""

from tqdm import tqdm

localrules: pls_flag

# expansion over pls
PLS = [int(i) for i in config['pls_components_for_loclinear'] ]

# load helper functions
hf = print_utils(False)

rule performance_analysis_of_neuralnetwork_regression:
    output: "data/performance/pod_{pod_iloc}/perf_sscomp_{sscomp}_method_neuralnet.table.gzip"
    input:
        abc_stats = "data/table/params_and_sumstats.table.gzip",
        pod_stats = "data/table/params_and_sumstats_pod_{pod_iloc}.table.gzip"
    params:
        script = "scripts/performance_analysis_of_neuralnetwork_regression.R",
        sumstats_to_use = hf.get_string_from_yaml_configuration_file_for_Rscript(
            config['sumstat_combination_to_use']),
        num_simulations_to_tolerate = config['number_of_tolerated_simulations'],
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

# for each summary statistic composition create the maximum number of
# components, but no more than the number given in the config file (usually 50)
rule get_pls_dim_reduced_sumstats_for_abc_sims:
    output:
        rmse_plot = "data/performance/pls_transformed/{sscomp}/RMSE_{sscomp}.pdf",
        transformer = "data/performance/pls_transformed/{sscomp}/Routput_{sscomp}"
    input:
        "data/table/params_and_sumstats.table.gzip"
    params:
        num_pls_max = config['number_max_pls_for_transformation'],
        sumstats_to_use = hf.get_string_from_yaml_configuration_file_for_Rscript(
            config['sumstat_combination_to_use']),
        script = "scripts/find_pls_by_wanted_sumstats.R",
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
    output: "data/table/txt/params_and_sumstats.table.txt"
    input: table = "data/table/params_and_sumstats.table.gzip"
    run:
        # use pandas to write the data frame as txt
        print("prepare.. writing table as txt")
        # pd.read_pickle(input.table, compression="gzip").to_csv(str(output), sep="\t", na_rep="NA", index=False)
        pd.read_pickle(input.table, compression=None).to_csv(
            str(output), sep="\t", na_rep="NA", index=False)
        print("done.. writing table as txt")

rule provide_pod_table_as_txt:
    output: "data/table/txt/params_and_sumstats_pod_{pod_iloc}.table.txt"
    input: table = "data/table/params_and_sumstats_pod_{pod_iloc}.table.gzip"
    group: "pls_transform_stats"
    run:
        # use pandas to write the data frame as txt
        print("prepare.. writing table as txt")
        pd.read_pickle(input.table).to_csv(str(output), sep="\t", na_rep="NA", index=False)
        print("done.. writing table as txt")


rule transform_sumstats_of_abc_sims:
    output:
        abc_sim_transformed = "data/performance/pls_transformed/{sscomp}/abc_sim.transformed.txt"
    input:
        pls = "data/performance/pls_transformed/{sscomp}/Routput_{sscomp}",
        abc_sim_raw = "data/table/txt/params_and_sumstats.table.txt"
    params:
        script = "scripts/transformer",
        script_subset = "scripts/subset_table_and_separate_params.py"
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
        pod_sim_transformed = "data/performance/pls_transformed/{sscomp}/params_and_sumstats_pod_{pod_iloc}_sim.transformed.txt"
    input:
        pls = "data/performance/pls_transformed/{sscomp}/Routput_{sscomp}",
        pod_sim_raw = "data/table/txt/params_and_sumstats_pod_{pod_iloc}.table.txt"
    params:
        script = "scripts/transformer",
        script_subset = "scripts/subset_table_and_separate_params.py"
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

rule performance_analysis_of_loclinear_regression:
    output:
        perf_ana = temp("data/performance/pod_{pod_iloc}/pls_{pls}/perf_sscomp_{sscomp}_method_loclinear.single.gzip"),
        point_estims = "data/point_estims/pod_{pod_iloc}/pls_{pls}/perf_sscomp_{sscomp}_method_loclinear.single.gzip",
        quants = "data/point_estims/pod_{pod_iloc}/pls_{pls}/quant_sscomp_{sscomp}_method_loclinear.single.gzip"
    input:
        transformed_abc_stats = "data/performance/pls_transformed/{sscomp}/abc_sim.transformed.txt",
        transformed_pod_stats = "data/performance/pls_transformed/{sscomp}/params_and_sumstats_pod_{pod_iloc}_sim.transformed.txt",
        pls_flag = "data/performance/pls_flag/{pod_iloc}_{sscomp}_{pls}.flag"
    params:
        script = "scripts/performance_analysis_using_loclinear_regression.R",
        num_simulations_to_tolerate = config['number_of_tolerated_simulations'],
        path_to_config_file = "config/config.yaml" # the name of the prior definition has to match the definition in the r-function
    threads: 1
    shell:
        r"""
        # please leave threads 1, otherwise writing between cores seems to slow down the process
        export OPENBLAS_NUM_THREADS={threads} OMP_NUM_THREADS={threads} MKL_NUM_THREADS={threads}

        # do the point estimation and perfromance for all requested pls
        Rscript --vanilla {params.script} {output.perf_ana} \
            {output.point_estims} \
            {output.quants} \
            {input.transformed_abc_stats} {input.transformed_pod_stats} \
            {params.num_simulations_to_tolerate} \
            {params.path_to_config_file} \
            {wildcards.pls}
        """

rule aggregate_quants:
    output:
        outfile = "data/table/quants.pod.performance.table.gzip"
    input:
        expand("data/point_estims/pod_{pod_iloc}/pls_{pls}/quant_sscomp_{sscomp}_method_loclinear.single.gzip",
            pod_iloc=pod_iloc,
            pls=PLS,
            sscomp=sscomp)
    run:
        # read in all tables
        df_list = []
        for infile in tqdm.tqdm(input, total=len(input), desc=" quant"):
            df = pd.read_pickle(infile)

            # add wildcards; only sscomp is left
            sscomp, _, method_containing_str = infile.split("_")[-3:]
            method = method_containing_str.split(".")[0]

            df["sscomp"] = sscomp
            df["method"] = method
            df["file"] = infile

            df_list.append(df)
            del df

        # this is not clean
        # we want to only use non empty lines
        df_list = [df for df in df_list if df.shape[0] > 1]

        df = pd.concat(df_list, axis=0, ignore_index=True)
        del df_list

        df.to_pickle(output.outfile, compression="gzip")


rule pls_flag:
    output: touch("data/performance/pls_flag/{pod_iloc}_{sscomp}_{pls}.flag")
    shell: r"""echo "provide {wildcards.pls} flag file" """

rule aggregate_performance_analysis_of_loclinear_regression_by_pls:
    output: "data/performance/pod_{pod_iloc}/perf_sscomp_{sscomp}_method_loclinear.table.gzip"
    input:
        expand("data/performance/pod_{{pod_iloc}}/pls_{pls}/perf_sscomp_{{sscomp}}_method_loclinear.single.gzip",
            pls = PLS)
    run:
        # concatenate the pickles
        df = pd.concat([pd.read_pickle(f) for f in input])
        df.to_pickle(str(output))
