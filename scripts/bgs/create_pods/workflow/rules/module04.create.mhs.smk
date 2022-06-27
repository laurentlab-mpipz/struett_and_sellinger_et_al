
rule mhs_from_tss:
    output:
        mhs = "results/mhs/pod_{pod_no}_tsigma_{tsigma}.mhs"
    input:
        tss = sumstat_on_pyslim_input_func
    log:
        std = "results/log/mhs/pod_{pod_no}_tsigma_{tsigma}.log"
    params:
        recombrate = float(config["recombination_rate"]),
        pop_size_outcrossing = int(float(set(config["pod_population_size_recent"]).pop())),
        nsam = int(float(config["sample_size_simulation"])),
    conda: "config/env.yml"
    threads: 1
    notebook:
        "../notebooks/create_mhs.py.ipynb"
