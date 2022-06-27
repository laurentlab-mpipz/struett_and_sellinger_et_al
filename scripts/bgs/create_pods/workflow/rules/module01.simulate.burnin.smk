
# expanding vars
tsigma_slim_asList = sorted([int(float(i)) for i in list(set(
    config["pod_selfing_rate_change_time"]))])
tsigma_slim = "c(" + ", ".join([str(int(i)) for i in tsigma_slim_asList]) + ")"

rule simulate_bgs_burnin:
    output:
        treeseq = "results/burnin/pod_{pod_no}.burnin.ts"  # within the slim script we define .burnin.ts as ending
    params:
        slim_script = "workflow/scripts/bgs_adapted_create_burnin.slim",
        mutation_rate = float(config['mutation_rate']),
        recombination_rate = float(config['recombination_rate']),
        locus_length = int(float(config["chromosome_length"])),
        pop_size = int(float(set(config["pod_population_size_recent"]).pop())),
        outfileprefix = lambda w, output: os.path.splitext(output[0])[0],
        model = "EX",
        sigma = float(set(config["pod_selfing_rate_recent"]).pop()),
        tsigma = tsigma_slim
    log:
        std = "results/log/burnin/pod_{pod_no}/burnin.log",
        err = "results/log/burnin/pod_{pod_no}/burnin.err"
    conda: "config/env.yml"
    resources:
        mem_mb = 20000
    threads: 1
    priority: 50
    shell:
        r"""
        # -m outputs memory usage and peak

        working_directory=\"`pwd -P`\"

        cat {params.slim_script} | slim  \
            -m \
            -d "mu='{params.mutation_rate}'" \
            -d "r='{params.recombination_rate}'" \
            -d "Lpre='{params.locus_length}'" \
            -d "Nepre='{params.pop_size}'" \
            -d "outfileprefix='{params.outfileprefix}'" \
            -d "model='{params.model}'" \
            -d "sigma='{params.sigma}'" \
            -d "tsigma={params.tsigma}" \
            -d "logfile='{log.std}'" \
            -d "working_directory=$working_directory" \
            2>{log.err} 

        wait
        """

rule simulate_bgs_transitioning_to_selfing:
    output:
        expand("results/bgs/pod_{{pod_no}}/tsigma_{tsigma}_locusno_{{locus}}.ts",
            tsigma=tsigma_slim_asList)
    input:
        burnin_treeseq = "results/burnin/pod_{pod_no}.burnin.ts"
    log:
        std = "results/log/bgs/pod_{pod_no}/burnin_locusno_{locus}.log",
        err = "results/log/bgs/pod_{pod_no}/burnin_locusno_{locus}.err"
    params:
        slim_script = "workflow/scripts/bgs_adapted_from_burnin.slim",
        mutation_rate = float(config['mutation_rate']),
        recombination_rate = float(config['recombination_rate']),
        locus_length = int(float(config["chromosome_length"])),
        pop_size = int(float(set(config["pod_population_size_recent"]).pop())),
        model = "EX",
        sigma = float(set(config["pod_selfing_rate_recent"]).pop()),
        tsigma = tsigma_slim,
        outfileprefix = lambda w, output: output[0].split("tsigma")[0] + "tsigma_",
        outfilesuffix = lambda w, output: "_locusno" + output[0].split("locusno")[1]
    conda: "config/env.yml"
    threads: 1
    shell:
        r"""
        # -m outputs memory usage and peak

        working_directory=\"`pwd -P`\"

        cat {params.slim_script} | slim  \
            -m \
            -d "mu='{params.mutation_rate}'" \
            -d "r='{params.recombination_rate}'" \
            -d "Lpre='{params.locus_length}'" \
            -d "Nepre='{params.pop_size}'" \
            -d "burninfile='{input}'" \
            -d "model='{params.model}'" \
            -d "sigma='{params.sigma}'" \
            -d "tsigma={params.tsigma}" \
            -d "logfile='{log.std}'" \
            -d "working_directory=$working_directory" \
            -d "outfileprefix='{params.outfileprefix}'" \
            -d "outfilesuffix='{params.outfilesuffix}'" \
            2>{log.err} 

        """
