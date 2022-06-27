
"""
produce a new table to summarize the model choice results
for each of the pod sets, for each pod under each given condition
calculate a bayes factor with and without transformation of parameters

1) using transformed data with the pls definition from model 1 (transition to selfing)
2) using Box-Cox transformed data (but no pls transformation)
"""

ruleorder: model_choice_without_pls > model_choice_after_pls_transformation

# expansion over pls; 0 means no pls transformation
PLS = [0] + [int(i) for i in config['pls_components_for_loclinear'] ]

# load helper functions
hf = print_utils(False)

rule model_choice_without_pls:
    output: "data/model_choice/pod_{pod_iloc}/sscomp_{sscomp}/model_choice_pls_0.table.gzip"
    input:
        abc_stats = "data/table/params_and_sumstats.table.gzip",
        alt_stats = "data/table/params_and_sumstats.alternate_model.table.gzip",
        pod_stats = "data/table/params_and_sumstats_pod_{pod_iloc}.table.gzip"
    params:
        script = "scripts/model_choice_without_pls_transformation.R",
        sumstats_to_use = hf.get_string_from_yaml_configuration_file_for_Rscript(
            config['sumstat_combination_to_use']),
        num_simulations_to_tolerate = config['model_choice_tolerance'],
        pls_comp = 0, # if zero: w/o pls transformation, only boxcox
    threads: 1
    shell:
        r"""
        export OPENBLAS_NUM_THREADS={threads} OMP_NUM_THREADS={threads} MKL_NUM_THREADS={threads}

        Rscript --vanilla {params.script} {output} \
            {input.abc_stats} {input.alt_stats} {input.pod_stats} \
            {params.num_simulations_to_tolerate} \
            {wildcards.sscomp} \
            "{params.sumstats_to_use}" \
            {params.pls_comp}
        """

rule model_choice_after_pls_transformation:
    output: "data/model_choice/pod_{pod_iloc}/sscomp_{sscomp}/model_choice_pls_{pls}.table.gzip"
    input:
        transformed_abc_stats = "data/performance/pls_transformed/{sscomp}/abc_sim.transformed.txt",
        transformed_alt_stats = "data/model_choice/pls_transformed/{sscomp}/alt_sim.transformed.txt",
        transformed_pod_stats = "data/performance/pls_transformed/{sscomp}/params_and_sumstats_pod_{pod_iloc}_sim.transformed.txt",
    params:
        script = "scripts/model_choice_after_pls_transformation.R",
        sumstats_to_use = hf.get_string_from_yaml_configuration_file_for_Rscript(
            config['sumstat_combination_to_use']),
        num_simulations_to_tolerate = config['model_choice_tolerance'],
        pls_flag = "{pls}"
    threads: 1
    shell:
        r"""
        export OPENBLAS_NUM_THREADS={threads} OMP_NUM_THREADS={threads} MKL_NUM_THREADS={threads}

        Rscript --vanilla {params.script} {output} \
            {input.transformed_abc_stats} {input.transformed_alt_stats} {input.transformed_pod_stats} \
            {params.num_simulations_to_tolerate} \
            {wildcards.sscomp} \
            "{params.sumstats_to_use}" \
            {wildcards.pls}
        """ 

rule pls_transform_alt_model:
    output:
        transformed_alt_stats = "data/model_choice/pls_transformed/{sscomp}/alt_sim.transformed.txt"
    input:
        pls_transformer = "data/performance/pls_transformed/{sscomp}/Routput_{sscomp}",
        alt_stats = "data/table/params_and_sumstats.alternate_model.table.gzip"
    params:
        script = "scripts/transformer",
        script_subset = "scripts/subset_table_and_separate_params.py"
    shadow: "shallow"
    threads: 1
    shell:
        r"""
        ## step 1: subset the table to transform by the columns identified for transformation and separate the params
        # last arg is the prefix for the output: subset.table[.param/.sumstat]
        #Rscript --vanilla {params.script_subset} {input.pls_transformer} {input.alt_stats} subset.table
        python {params.script_subset} {input.pls_transformer} {input.alt_stats} subset.table
        head -n 6 subset.table* | cut -f 1-10
        echo "..."
        echo "_______________"
        echo "finished step 1"
        echo "==============="

        ## step 2: transform
        # output file prefix here is simply output.transformed
        {params.script} {input.pls_transformer} subset.table.sumstat {output.transformed_alt_stats} boxcox
        echo "_______________"
        echo "finished step 2"
        echo "==============="
        """

rule aggregate_model_choice_by_pls:
    output: "data/table/model_choice_pod_aggregation.table.gzip"
    input:
        expand("data/model_choice/pod_{pod_iloc}/sscomp_{sscomp}/model_choice_pls_{pls}.table.gzip",
            pls=PLS, pod_iloc=pod_iloc, sscomp=sscomp)
    run:
        # concatenate the pickles
        def get_pod_sscomp_pls_from_path(path_to_file):
            a = ["".join(x) for _, x in itertools.groupby(path_to_file, key=str.isdigit)]
            for ix,x in enumerate(a):
                if "pod" in x:
                    pod_iloc = a[ix+1]
                elif "sscomp" in x:
                    sscomp = a[ix+1]
                elif "pls" in x:
                    pls = a[ix+1]
            return int(pod_iloc), int(sscomp), int(pls)

        def read_table_and_add_wildcards_to_columns(path_to_file):
            df = pd.read_pickle(path_to_file)
            pod_iloc, sscomp, pls = get_pod_sscomp_pls_from_path(path_to_file)
            df['pod'] = pod_iloc
            df['sscomp'] = sscomp
            df['pls'] = pls
            df['path'] = path_to_file
            return df

        df = pd.concat([read_table_and_add_wildcards_to_columns(f) for f in tqdm.tqdm(input)])
        df.to_pickle(str(output))
